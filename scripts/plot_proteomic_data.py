#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import holoviews as hv
import numpy as np
hv.extension("matplotlib")
import matplotlib.pyplot as plt
from holoviews import opts


def plot_proteomic_data(annotations, abundances):
    """
    Plot the functional groups, genomic and abundance distribution for proteomic data
    """

    bars = hv.Bars(annotations, hv.Dimension("Functional Annotation"), "# Pfam Domains")
    frequencies, edges = np.histogram(abundances, 10)
    hist = hv.Histogram((edges, frequencies))
    tem = hv.RGB.load_image('resources/sfgv_ob_tem.png')
    sem = hv.RGB.load_image('resources/sfgv_ob_sem.png')
    layout = hv.Layout(tem + sem + bars + hist).cols(2)
    layout.opts(
        opts.RGB(xaxis=None, yaxis=None),
        opts.Bars(color="gray"),
        opts.Histogram(facecolor="gray", xlabel="emPAI %VP39", ylabel="# OB Proteins"),
            )
    opts.defaults(opts.Layout(hspace=0.1, vspace=0.1))
    hv.save(layout, "layout.svg", backend="matplotlib", dpi=600)

# Command line interface
@click.command()
@click.option("-a", "--annotations",
              help="Pfam functional annotations")
@click.option("-p", "--protein_list",
              help="Protein list with abundance values")
@click.option("-f", "--proteome",
              help="Fasta file with protein positions")
@click.option("-o", "--output",
              help="Output file")
def cli(annotations, protein_list, proteome, output):
    "Command line interface"

    names = ["Protein",
             "InterProscan",
             "Description",
             "Annotation"]
    ann = pd.read_csv(annotations,
                      sep="\t",
                      names=names)
    g = ann.groupby("Annotation").size()
    ann = list(zip(g.index, g.values))

    df = pd.read_csv(protein_list)
    vp39 = float(df[df["Protein name"] == "Vp39"]["Average emPAI"])
    abund = df["Average emPAI"] * 100 / vp39
    abund = abund[abund <= 100]

    plot_proteomic_data(ann, abund)


if __name__ == "__main__":
    cli()


# def get_protein_locations(proteome):
#     """
#     Parse protein positions from a fasta file and returns them as a list
#     """

#     # Load proteome
#     seqs = SeqIO.parse(proteome, "fasta")
#     # Initialize coordinates list
#     locs = []
#     # Get location for each protein on the file
#     for seq in seqs:
#         loc = seq.description.split("[location=")[1].split("]")[0]
#         if loc.startswith("complement"):
#             locs.append(int(loc[11:-1].split("..")[1]))
#         else:
#             locs.append(int(loc.split("..")[0]))

#     return locs
