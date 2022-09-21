#!/usr/bin/env python3.8

"""
"""

import click
import gzip
import pandas as pd
import requests


def download_proteome(proteome, identifier):
    """
    """

    base_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Viruses/"
    endpoint = f"{proteome}/{proteome}_{identifier}.fasta.gz"
    url = base_url + endpoint
    # url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Viruses/{proteome}/{proteome}_{identifier}.fasta.gz"
    output = f"{proteome}.gz"

    r = requests.get(url, stream=True)
    with open(output, "wb") as fh:
        for chunk in r.raw.stream(1024, decode_content=False):
            fh.write(chunk)


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-p",
              "--proteomes_list",
              help="")
# @click.option("-i",
#               "--identifier",
#               help="")
def cli(proteomes_list):
    """
    """

    df = pd.read_csv(proteomes_list, sep="\t")
    proteomes = df["Proteome Id"]
    identifiers = df["Organism Id"]
    for p, i in zip(proteomes, identifiers):
        download_proteome(p, i)


if __name__ == '__main__':
    cli()
