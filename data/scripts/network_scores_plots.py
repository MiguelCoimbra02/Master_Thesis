#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script loads a network.txt and generates multiple visualizations for the distribution of the aggregated scores.
Options include histogram, box plot, and density plot.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def network_histogram(network, score, outname):
    df = pd.read_table(network, sep=" ")
    fig, ax = plt.subplots()
    df.hist(column=score, ax=ax)
    plt.xlabel(score)
    plt.ylabel('Number of Edges')
    plt.title(f"{outname} - Histogram")
    fig.savefig(f"{outname}_histogram.png")

def network_boxplot(network, score, outname):
    df = pd.read_table(network, sep=" ")
    fig, ax = plt.subplots()
    df.boxplot(column=score, ax=ax)
    plt.xlabel('Network')
    plt.ylabel(score)
    plt.title(f"{outname} - Box Plot")
    fig.savefig(f"{outname}_boxplot.png")

def network_density_plot(network, score, outname):
    df = pd.read_table(network, sep=" ")
    fig, ax = plt.subplots()
    df[score].plot(kind='density', ax=ax)
    plt.xlabel(score)
    plt.ylabel('Density')
    plt.title(f"{outname} - Density Plot")
    fig.savefig(f"{outname}_density_plot.png")

def main():
    parser = argparse.ArgumentParser(description="Generate multiple visualizations (Histogram, Box Plot, Density Plot) for aggregated network scores.")
    parser.add_argument('--network', required=True, metavar='tab', help='Network file from which the scores will generate visualizations')
    parser.add_argument('--score', type=str, metavar="STR", default="irp_score", help='Score column to visualize (e.g., irp_score)')
    parser.add_argument('--out', type=str, metavar="STR", default="network_visualization", help='Output file name base')
    parser.add_argument('--visualizations', type=str, choices=['histogram', 'boxplot', 'density'], nargs='+', default=['histogram', 'boxplot', 'density'], 
                        help='List of visualizations to generate. Options: histogram, boxplot, density')

    args = parser.parse_args()

    # Generate selected visualizations
    if 'histogram' in args.visualizations:
        network_histogram(args.network, args.score, args.out)
    if 'boxplot' in args.visualizations:
        network_boxplot(args.network, args.score, args.out)
    if 'density' in args.visualizations:
        network_density_plot(args.network, args.score, args.out)
main()