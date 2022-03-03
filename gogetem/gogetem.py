import argparse
import requests
import time
import gzip
import hashlib
from pathlib import Path
from collections.abc import Generator, Iterable
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON

from gogetem import __version__


def arguments():

    parser = argparse.ArgumentParser()

    download_location = parser.add_mutually_exclusive_group(required=True)

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

    download_location.add_argument(
        "--download-path",
        type=Path,
        help="Results directory. Will overwrite existing data found at that location.",
    )

    download_location.add_argument(
        "--resume",
        type=Path,
        help="If this path points to a previously-specified path given by --download-path,"
        " and uniprot_results.tsv has been generated under it, gogetem will attempt to"
        " resume downloading nucleotide data",
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

    results_path = args.download_path or args.resume

    if args.resume is None:  # if no prior results exist

        uniprot_query = query_build(args.go_terms, args.include_amino_acids, args.limit)

        uniprot_results = query_submit(uniprot_query)
        uniprot_results_table = parse_results(uniprot_results)
        uniprot_table_save(uniprot_results_table, args.download_path)

        args.download_path.joinpath("query.txt").write_text(uniprot_query)

    # load previously generated table if the --resume option is used
    else:
        uniprot_results_table = pd.read_csv(
            args.resume.joinpath("uniprot_results.tsv"), sep="\t"
        )

    ena_retrieve(uniprot_results_table, results_path)
    amino_acids_write(uniprot_results_table, results_path)


def ena_retrieve(uniprot_results_table: pd.DataFrame, download_path: Path):

    ena_queries = ena_query_format(uniprot_results_table)
    failed_queries = ena_download_accessions(ena_queries, download_path)

    while failed_queries:

        n_failed = len(failed_queries)
        current_time = time.strftime("%Y-%m-%d %H:%M:%S")

        print(f"{current_time} | Reattempting to download {n_failed} failed queries")

        time.sleep(30)
        failed_queries = ena_download_accessions(failed_queries, download_path)


def amino_acids_write(uniprot_results_table: pd.DataFrame, download_path: Path):

    aa_download_path = download_path.joinpath("amino_acid")
    aa_download_path.mkdir(parents=True, exist_ok=True)

    fasta_path = aa_download_path.joinpath("amino_acids.fasta")

    formatted_aa_fasta = amino_acids_format(uniprot_results_table)

    with fasta_path.open("w") as f:
        f.writelines(formatted_aa_fasta)


def amino_acids_format(
    uniprot_results_table: pd.DataFrame,
) -> Generator[str, None, None]:

    fasta_format = (
        amino_acid_format(row) for _, row in uniprot_results_table.iterrows()
    )

    return fasta_format


def amino_acid_format(uniprot_row: pd.Series) -> str:

    versioned_accession = uniprot_row["ena_accession"]
    versionless_accession, _ = versioned_accession.split(".")
    aa_sequence = uniprot_row["aa_sequence"]
    description = uniprot_row["name"]

    return f">ENA|{versionless_accession}|{versioned_accession}|{description}\n{aa_sequence}\n"


def uniprot_table_save(uniprot_results_table: pd.DataFrame, download_path: Path):

    download_path.mkdir(exist_ok=True)
    uniprot_results_table.to_csv(
        download_path.joinpath("uniprot_results.tsv"), sep="\t", index=False
    )


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
    go_terms: list[str], include_amino_acids: bool = False, limit: int = 0
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

    current_accessions: list[str] = []
    query_length = 0

    for accession in results_table["ena_accession"]:
        acc_len = len(accession)

        if query_length + acc_len < 1000:
            current_accessions.append(accession)
            query_length += acc_len

        else:
            accession_list = ",".join(current_accessions)
            yield accession_list

            current_accessions = [accession]
            query_length = acc_len
    else:
        accession_list = ",".join(current_accessions)
        yield accession_list


def ena_fetch(ena_query: str) -> str:
    # TODO: consider changing this a POST request
    base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
    query = f"{base_url}{ena_query}"
    response = requests.get(query, params={"download": "true", "gzip": "true"})

    return gzip.decompress(response.content).decode()


def ena_download_accessions(
    ena_queries: Iterable[str], download_path: Path
) -> list[str]:
    # uses hash of query url as the file name, so it can tell if a particular
    # file has been downloaded already or not
    #
    # returns a list of failed queries
    failed: list[str] = []

    nt_download_path = download_path.joinpath("nucleotide")
    nt_download_path.mkdir(parents=True, exist_ok=True)

    for ena_query in ena_queries:

        fasta_name = hashlib.md5(ena_query.encode()).hexdigest()
        fasta_path = nt_download_path.joinpath(f"{fasta_name}.fasta")

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
        [int(x) for x in go_term]
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main()
