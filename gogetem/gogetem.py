import argparse
import requests
import gzip
from pathlib import Path
from typing import Generator
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--go-terms", nargs="+", type=int, help="GO terms to search for"
    )

    args = parser.parse_args()

    return args


def main():

    args = arguments()


def query_prefix() -> str:

    prefixen = [
        "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>",
        "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>",
        "PREFIX up: <http://purl.uniprot.org/core/>",
        "PREFIX go: <http://purl.obolibrary.org/obo/GO_>",
    ]

    return "\n".join(prefixen)


def query_select(include_amino_acids: bool = False) -> str:

    aa_var = "?aa_sequence" if include_amino_acids else ""
    select = f"SELECT ?protein ?name ?link {aa_var}"
    return select


def query_match(go_terms: list[int], include_amino_acids: bool = False) -> str:

    go_terms_str = " ".join(["{"] + [f"go:{g}" for g in go_terms] + ["}"])

    get_aa_sequence = """
        ?protein up:sequence ?seq .
        ?seq rdf:value ?aa_sequence . 
    """

    aa_sequence_query = get_aa_sequence if include_amino_acids else ""

    where = [
        "WHERE {",
        f"values ?go_terms {go_terms_str}",
        "?protein a up:Protein ;",
        "    rdfs:label ?name ;",
        "    up:classifiedWith|(up:classifiedWith/rdfs:subClassOf) ?go_terms .",
        aa_sequence_query,
        "?protein rdfs:seeAlso ?link .",
        "?link up:database ?database .",
        "?database rdfs:label 'EMBL nucleotide sequence database' . ",
        "}",
    ]

    return "\n".join(where)


def query_build(go_terms: list[int], include_amino_acids: bool = False) -> str:

    prefixen = query_prefix()
    select_stmt = query_select(include_amino_acids)
    match_stmt = query_match(go_terms, include_amino_acids)

    query = "\n".join([prefixen, select_stmt, match_stmt])

    return query


def query_submit(query: str) -> list[dict[str, dict[str, str]]]:

    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql")
    sparql.setReturnFormat(JSON)
    sparql.setQuery(query)
    ret = sparql.queryAndConvert()

    return ret["results"]["bindings"]


def parse_results(results: list[dict[str, dict[str, str]]]) -> pd.DataFrame:

    flattened = []
    for record in results:
        simplified_record = {}
        for field, value in record.items():
            simplified_record[field] = value["value"]
        flattened.append(simplified_record)

    df = pd.DataFrame(flattened)
    df["ena_accession"] = [Path(p).name for p in df["link"]]

    return df


def ena_query_format(results_table: pd.DataFrame) -> Generator[str, None, None]:

    base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
    ena_queries: list[str] = []

    current_accessions: list[str] = []
    query_length = len(base_url)

    for accession in results_table["ena_accession"]:

        if query_length + len(accession) < 1000:
            current_accessions.append(accession)

        else:
            accession_list = ",".join(current_accessions)
            ena_queries.append(f"{base_url}{accession_list}")
            current_accessions = []
            query_length = len(base_url)
    else:
        accession_list = ",".join(current_accessions)
        ena_queries.append(f"{base_url}{accession_list}")

    return ena_queries


def ena_fetch(ena_query: str) -> str:

    response = requests.get(ena_query, params={"download": "true", "gzip": "true"})

    return gzip.decompress(response.content).decode()


if __name__ == "__main__":
    main()
