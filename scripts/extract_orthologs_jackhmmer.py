#!/usr/bin/env python3

"Download Drosophila melanogaster extracellular domain batabase (FlyXCDB) table, published in the Journal of Molecular Biology"

import click
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO


def extract_sequences_jackhmmer(sequences, jackhmmer, output):
    """
    """

    seqs = SeqIO.to_dict((SeqIO.parse(sequences, "fasta")))
    search = SearchIO.read(jackhmmer, "hmmer3-tab")

    with open(output, "w") as fh:
        for i in search:
            data = seqs[i.id]
            header = data.description.split("[")[1][:-1]
            header = i.id + "_" + "_".join(header.split())
            sequence = data.seq
            fh.write(f">{header}\n{sequence}\n")


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-s",
              "--sequences",
              help="")
@click.option("-j",
              "--jackhmmer",
              help="")
@click.option("-o",
              "--output",
              help="")
# CLI main function
def cli(sequences, jackhmmer, output):
    """
    """

    extract_sequences_jackhmmer(sequences, jackhmmer, output)


if __name__ == '__main__':
    cli()
