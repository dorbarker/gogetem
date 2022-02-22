import argparse
import requests
import time
import gzip
import hashlib
from pathlib import Path
from typing import Generator
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON

from gogetem import __version__


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--go-terms", nargs="+", type=str, help="GO terms to search for"
    )

    parser.add_argument(
        "--limit", type=int, default=0, help="Number of records to retrieve [no limit]"
    )

    parser.add_argument(
        "--include-amino-acids",
        action="store_true",
        help="Download amino acid sequences from UniProt",
    )

    parser.add_argument(
        "--download-path", type=Path, required=True, help="Results directory"
    )

    parser.add_argument(
        "-v", "--version", action="version", version=f"gogetem {__version__}"
    )

    args = parser.parse_args()

    # validate GO terms
    for go_term in args.go_terms:
        if not is_valid_go_term(go_term):
            invalid_go_term = f"{go_term} is not a valid numeric GO term."
            raise argparse.ArgumentTypeError(invalid_go_term)

    return args


def main():

    args = arguments()

    uniprot_query = query_build(args.go_terms, args.include_amino_acids, args.limit)
    uniprot_results = query_submit(uniprot_query)
    uniprot_results_table = parse_results(uniprot_results)

    ena_queries = ena_query_format(uniprot_results_table)
    ena_download_accessions(ena_queries, args.download_path)


def ena_retrieve():
    pass


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


def query_match(
    go_terms: list[str], include_amino_acids: bool = False, limit: int = 0
) -> str:

    go_terms_str = " ".join(["{"] + [f"go:{g}" for g in go_terms] + ["}"])

    get_aa_sequence = """
        ?protein up:sequence ?seq .
        ?seq rdf:value ?aa_sequence . 
    """

    aa_sequence_query = get_aa_sequence if include_amino_acids else ""
    limit_keyword = f"LIMIT {limit}" if limit else ""

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
        limit_keyword,
    ]

    return "\n".join(where)


def query_build(
    go_terms: list[int], include_amino_acids: bool = False, limit: int = 0
) -> str:

    prefixen = query_prefix()
    select_stmt = query_select(include_amino_acids)
    match_stmt = query_match(go_terms, include_amino_acids, limit)

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

    current_accessions: list[str] = []
    query_length = len(base_url)

    for accession in results_table["ena_accession"]:

        if query_length + len(accession) < 1000:
            current_accessions.append(accession)

        else:
            accession_list = ",".join(current_accessions)
            yield f"{base_url}{accession_list}"
            current_accessions = []
            query_length = len(base_url)
    else:
        accession_list = ",".join(current_accessions)
        yield f"{base_url}{accession_list}"


def ena_fetch(ena_query: str) -> str:
    # TODO: consider changing this a POST request
    response = requests.get(ena_query, params={"download": "true", "gzip": "true"})

    return gzip.decompress(response.content).decode()


def ena_download_accessions(
    ena_queries: Generator[str, None, None], download_path: Path
) -> list[str]:
    # uses hash of query url as the file name, so it can tell if a particular
    # file has been downloaded already or not
    #
    # returns a list of failed queries
    failed: list[str] = []

    download_path.mkdir(parents=True, exist_ok=True)

    for ena_query in ena_queries:

        fasta_name = hashlib.md5(ena_query.encode()).hexdigest()
        fasta_path = download_path.joinpath(f"{fasta_name}.fasta")

        if fasta_path.exists():  # skip files we've already made
            continue

        fasta = ena_fetch(ena_query)
        if fasta:
            fasta_path.write_text(fasta)
        else:
            # as of right now, a failed query returns an empty string
            failed.append(ena_query)

        time.sleep(1)  # DBAD

    return failed


def is_valid_go_term(go_term: str) -> bool:
    # Leading 0s in GO terms are significant,
    # so we can't just test if the whole thing is an int
    try:
        return all([int(x) for x in go_term])
    except ValueError:
        return False


if __name__ == "__main__":
    main()
