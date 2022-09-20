#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def plot_proteomic_data(annotations, abundances, output):
    """
    Plot the functional groups, genomic and abundance distribution for proteomic data
    """

    figure, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
    tem = plt.imread('../resources/sfgv_ob_tem.png')
    axs[0, 0].imshow(tem)
    axs[0, 0].axis("off")
    sem = plt.imread('../resources/sfgv_ob_sem.png')
    axs[0, 1].imshow(sem)
    axs[0, 1].axis("off")
    bars = sns.barplot(x="Annotation", y="# Pfam Domains",
                       data=annotations, ax=axs[1, 0],
                       palette="deep")
    hist = sns.histplot(abundances, ax=axs[1, 1], color="black", bins=10)
    plt.savefig(f"{output}.png")


# Command line interface
@click.command()
@click.option("-a", "--annotations",
              help="Pfam functional annotations")
@click.option("-p", "--protein_list",
              help="Protein list with abundance values")
@click.option("-o", "--output",
              help="Output file")
def cli(annotations, protein_list, output):
    "Command line interface"

    names = ["Protein",
             "InterProscan",
             "Description",
             "Annotation"]
    ann = pd.read_csv(annotations,
                      sep="\t",
                      names=names)
    g = ann.groupby("Annotation").size()
    # ann = list(zip(g.index, g.values))
    ann = pd.DataFrame(data={"Annotation": g.index, "# Pfam Domains": g.values})
    # ann = list(zip(g.index, g.values))

    df = pd.read_csv(protein_list)
    vp39 = float(df[df["Protein name"] == "Vp39"]["Average emPAI"])
    abund = df["Average emPAI"] * 100 / vp39
    abund = abund[abund <= 100]

    plot_proteomic_data(ann, abund, output)


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
