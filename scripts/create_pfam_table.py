#!/usr/bin/env python3

"Extract Pfam annotations from hmmscan output file"

import click
from Bio import SearchIO


def create_pfam_table(input, output):
    """
    Extracts the Pfam hits from a hmmscan file into a table format
    """

    # Load hmmscan results
    results = SearchIO.parse(input, "hmmer3-text")
    # Iterate over the records
    with open(output, "w") as fh:
        for res in results:
            # Save query name
            query = res.id.split("_")[2]
            # Write Pfam hit
            if len(res.hits) > 0:
                fh.write(f"{query}\t{res.hits[0].id}\n")
            else:
                fh.write(f"{query}\tNo hit\n")


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-i",
              "--input",
              help="Input search file")
@click.option("-o",
              "--output",
              help="Output table name")
# CLI main function
def cli(input, output):
    """
    Command line interface
    """

    create_pfam_table(input, output)


if __name__ == "__main__":
    cli()
