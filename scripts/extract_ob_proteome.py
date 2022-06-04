#!/usr/bin/env python3

import click
import pandas as pd
from Bio import SeqIO


def extract_ob_proteome(proteome, protein_list, output):
    """
    Filters a multifasta proteome and returns only the
    proteins present on the list provided
    """

    # Load proteome dictionary
    seqs = SeqIO.to_dict(SeqIO.parse(proteome, "fasta"),
                         key_function=lambda rec: rec.id.split("_")[2])
    # Load protein list
    df = pd.read_csv(protein_list)
    with open(output, "w") as fh:
        for item in df.Accession:
            if item in seqs.keys():
                seq = seqs[item]
                fh.write(f">{seq.description}\n{seq.seq}\n")


# Command line interface
@click.command()
@click.option("-p", "--proteome",
              help="Multifasta proteome")
@click.option("-l", "--protein_list",
              help="Protein list to retain")
@click.option("-o", "--output",
              help="Output file")
def cli(proteome, protein_list, output):
    "Command line interface"
    extract_ob_proteome(proteome, protein_list, output)


if __name__ == "__main__":
    cli()
