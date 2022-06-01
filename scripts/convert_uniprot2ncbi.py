#!/usr/bin/env python3

import click
import pandas as pd


def convert_uniprot2ncbi(mapping, proteome, output):
    """
    Remaps the UniProt identifiers in a protein table
    into the NCBI accession number, based on a
    mapping table (e.g. Blast hits table)
    """

    # Load uniprot2ncbi mappings table
    columns = ["uniprot",
               "ncbi",
               "%_identity",
               "length",
               "mismatch",
               "gapopen",
               "qstart",
               "qend",
               "sstart",
               "send",
               "p-value",
               "bitscore"
               ]
    df = pd.read_csv(mapping,
                     sep="\t",
                     names=columns)
    # Extract UniProt accessions
    df.uniprot = df.uniprot.str.split("|", expand=True)[1]
    # Extract NCBI accessions
    df.ncbi = df.ncbi.str.split("_", expand=True)[2]
    # Create a UniProt -> NCBI dictionary
    dic = dict(zip(df.uniprot, df.ncbi))
    # Load proteome dataset
    data = pd.read_csv(proteome)
    # Replace UniProt with NCBI
    data = data.replace(dic)
    # Drop duplicates know to be a problem
    drop_list = ["ORF038",
                 "Vp91/p95",
                 "ORF108"]
    data = data[~data["Protein name"].isin(drop_list)]
    # Save data
    data.to_csv(f"{output}", index=False)


# Command line interface
@click.command()
@click.option("-m", "--mappings")
@click.option("-p", "--proteome")
@click.option("-o", "--output")
def cli(mappings, proteome, output):
    "Command line interface"
    convert_uniprot2ncbi(mappings, proteome, output)


if __name__ == "__main__":
    cli()
