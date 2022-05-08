#!/usr/bin/env python3

"Aggregates individual SfGV proteomics datasets into a summary table"

import click
import pandas as pd


def create_proteins_table(dataset1, dataset2, dataset3, output):
    """
    Merges the three protein tables of each replicate into a summary table.
    """

    # Read the three datasets and parse protein names
    df1 = pd.read_csv(dataset1)
    df1["Protein"] = df1["Description"].str.split(" ", expand=True)[0]
    df2 = pd.read_csv(dataset2)
    df2["Protein"] = df2["Description"].str.split(" ", expand=True)[0]
    df3 = pd.read_csv(dataset3)
    df3["Protein"] = df3["Description"].str.split(" ", expand=True)[0]
    # Merge datasets
    summary = pd.concat([df1, df2, df3])
    # Get de-duplicated accession IDs
    accessions = summary["Accession"].drop_duplicates()
    # Get protein names
    protein_names = summary["Protein"].drop_duplicates().set_axis(accessions)
    # Get protein lengths
    aa = summary.groupby("Accession")["# AAs"].mean()
    # Get protein coverages
    coverage = summary.groupby("Accession")["Coverage"].mean()
    # Get protein emPAIs
    empai = summary.groupby("Accession")["emPAI"].mean()
    # Write summary table
    table = pd.DataFrame({"Protein name": protein_names,
                          "Protein Size (amino acids)": aa,
                          "Average Coverage (%)": coverage,
                          "Average emPAI": empai})
    table = table.round(decimals=2)
    table.to_csv(output)


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-d1",
              "--data1",
              help="First dataset")
@click.option("-d2",
              "--data2",
              help="Second dataset")
@click.option("-d3",
              "--data3",
              help="Third dataset")
@click.option("-o",
              "--output",
              help="Output table name")
# CLI main function
def cli(data1, data2, data3, output):
    """
    Command line interface
    """

    create_proteins_table(data1, data2, data3, output)


if __name__ == "__main__":
    cli()
