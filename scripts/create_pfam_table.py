#!/usr/bin/env python3

"Extract Pfam annotations from hmmscan output file"

import click
import pandas as pd
from Bio import SearchIO


def create_pfam_table(input, go_terms, output):
    """
    Extracts the Pfam hits from a hmmscan file into a table format
    """

    # Load pfam2go mappings
    go_terms = pd.read_table(go_terms,
                             sep=" > ",
                             skiprows=6,
                             names=["Pfam", "GO"],
                             engine="python")
    pfam_dom = go_terms.Pfam.str.split(" ", expand=True)[1]
    go_ids = go_terms.GO.str.split("; ", expand=True)[1]
    pfam2go = pd.DataFrame({"Pfam": pfam_dom,
                            "GO": go_ids})

    # Load hmmscan results
    results = SearchIO.parse(input, "hmmer3-text")
    # Iterate over the records
    with open(output, "w") as fh:
        for res in results:
            # Save query name
            query = res.id.split("_")[2]
            # Get Pfam and GO hits
            if len(res.hits) > 0:
                pfam = res.hits[0].id
                go = ",".join(pfam2go[pfam2go.Pfam == pfam].GO)
            # Set negative results to "No hit"
                if len(go) == 0:
                    go = "No hit"
            else:
                pfam = "No hit"
                go = "No hit"
            # Write to output file
            fh.write(f"{query}\t{pfam}\t{go}\n")


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-i",
              "--input",
              help="Input search file")
@click.option("-g",
              "--go_terms",
              help="Pfam to GO terms mapping file")
@click.option("-o",
              "--output",
              help="Output table name")
# CLI main function
def cli(input, go_terms, output):
    """
    Command line interface
    """

    create_pfam_table(input, go_terms, output)


if __name__ == "__main__":
    cli()
