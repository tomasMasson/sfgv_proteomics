#!/usr/bin/env python3

"""
Parse individual JSON files from Hyphy FUBAR and aBSREL methods into a TSV row
"""

import click
import csv
import json
import pandas as pd
from pathlib import Path


def load_fubar(fubar):
    """
    Read FUBAR JSON file
    """

    # Load file
    with open(fubar, "r") as fh:
        data = json.load(fh)
    # Read column headers and contents
    headers = [h[0] for h in data["MLE"]["headers"]]
    # Intialize a DataFrame for further calculations
    content = data["MLE"]["content"]["0"]
    # Drop empty columns and assign headers
    df = pd.DataFrame(content).drop([6, 7], axis=1)
    df.columns = headers

    return df


def compute_fubar_feaures(fubar):
    """
    Calculates a set of features from a FUBAR JSON file
    """

    df = load_fubar(fubar)
    # Create a dN/dS (Omega) column
    df["omega"] = df["beta"] / df["alpha"]
    # Define the features to be saved from the DataFrame
    alpha = df["alpha"].mean()
    beta = df["beta"].mean()
    omega = df["omega"].mean()
    # Fractions of negatively and positively selected residues
    f_neg = df[df["Prob[alpha>beta]"] > 0.9].shape[0] / df.shape[0]
    n_pos = df[df["Prob[alpha<beta]"] > 0.9].shape[0]
    # Store features into a list before returning them
    features = [alpha, beta, omega, f_neg, n_pos]

    return features


def load_meme(meme):
    """
    Read MEME JSON file
    """

    # Load file
    with open(meme, "r") as fh:
        data = json.load(fh)
    # Read column headers and contents
    headers = [h[0] for h in data["MLE"]["headers"]]
    # Intialize a DataFrame for further calculations
    content = data["MLE"]["content"]["0"]
    # Drop empty columns and assign headers
    df = pd.DataFrame(content)
    df.columns = headers

    return df


def compute_meme_feaures(meme):
    """
    Calculates a set of features from a FUBAR JSON file
    """

    df = load_meme(meme)
    n_pos = df[df["p-value"] < 0.1].shape[0]

    return [n_pos]


# def merge_fubar_absrel(fubar, absrel):
def merge_fubar_meme(fubar, meme):
    """
    Merge into a single list the orthogroup name,
    FUBAR and aBSREL results
    """

    # Get OG from file's name
    og = Path(fubar).name.split(".")[0]
    ff = compute_fubar_feaures(fubar)
    fm = compute_meme_feaures(meme)
    # branches = parse_absrel_branches(absrel)

    return [og] + ff + fm
    # return [og] + features + [branches]


@click.command()
@click.option("-f",
              "--fubar",
              help="FUBAR results from Hyphy")
@click.option("-m",
              "--meme",
              help="MEME results from Hyphy")
@click.option("-o",
              "--output",
              help="Output file")
# def cli(fubar, absrel, output):
def cli(fubar, meme, output):
    """
    Parse individual JSON files from Hyphy FUBAR and aBSREL methods into a TSV row
    """

    # Get results
    results = merge_fubar_meme(fubar, meme)
    # Save results as a single row (tab-delimited)
    with open(output, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(results)


if __name__ == "__main__":
    cli()
