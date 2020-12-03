#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel

version = "0.0.1"
##################################################
# Modules
##################################################

try:
    from yoda_powers.toolbox import existant_file, sort_human, welcome_args
except ImportError:
    print(f"Best module 'yoda_powers' not install, please check:'https://yoda-powers.readthedocs.io/en/latest/'")
    exit()
import sys
from pathlib import Path
import argparse
from datetime import datetime
from collections import defaultdict, OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

plt.style.use('ggplot')
mpl.use('Agg')


def build_parser():
    epilog_tools = """Documentation avail at: https://yoda-powers.readthedocs.io/en/latest/ \n\n"""
    description_tools = f"""
    More information:
        Script version: {version}
    """
    parser_mandatory = argparse.ArgumentParser(add_help=False)

    mandatory = parser_mandatory.add_argument_group('Input mandatory infos for running')
    mandatory.add_argument('-l', '--lib', metavar="path/to/file/csv", type=existant_file, required=True,
                           dest='csv_file', help='path to file with size library, use to normalize data')
    parser_other = argparse.ArgumentParser(
            parents=[parser_mandatory],
            add_help=False,
            prog=Path(__file__).name,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=description_tools,
            epilog=epilog_tools
    )

    optional = parser_other.add_argument_group('Input infos not mandatory')
    optional.add_argument('-v', '--version', action='version', version=version,
                          help=f'Use if you want to know which version of {Path(__file__).name} you are using')
    optional.add_argument('-h', '--help', action='help', help=f'show this help message and exit')
    optional.add_argument('-d', '--debug', action='store_true', help='enter verbose/debug mode')
    optional.add_argument('-p', '--plot', action='store_true', help='plot connections distribution')
    optional.add_argument('-c', '--chromosome', metavar="int", type=int, default=1000000, dest='scaff_min',
                          help='Minimum of scaffold size (default = 1000000)')
    optional.add_argument('-f', '--fragments', metavar="int", type=int, default=5000, dest='fragments_min',
                          help='Minimum of connection size (default = 5000)')

    return parser_other


@welcome_args(version, build_parser())
def main():
    prog_args = build_parser().parse_args()
    # GLOBAL VARIABLE
    sep = "\t"
    sample1, sample2, *_ = prog_args.csv_file.name.split(".pp1.fasta.")
    print(sample1, sample2)

    # open mummer dataframe
    header = [f"{sample2}-start1", f"{sample2}-stop1", f"{sample1}-start2", f"{sample1}-stop2", "len_1", "len_2", "t7",
              "t8", "t9", f"{sample2}_len", f"{sample1}_len",
              "t12", "t13", f"{sample2}-chrom", f"{sample1}-chrom"]
    df = pd.read_csv(prog_args.csv_file, sep=sep, names=header)
    # df[f"{sample1}-chrom"] = df[f"{sample1}-chrom"].replace("Scaffold_", f"{sample1}_")
    df[f"{sample1}-chrom"] = [item.replace("Scaffold_", f"{sample1}_") for item in df[f"{sample1}-chrom"]]
    df[f"{sample2}-chrom"] = [item.replace("Scaffold_", f"{sample2}_") for item in df[f"{sample2}-chrom"]]
    # df[f"{sample2}-chrom"] = df[f"{sample2}-chrom"].replace("Scaffold_", f"{sample2}_")
    if prog_args.debug: print(df)
    # exit()
    # see statistics of connections
    print(f"\n{'*' * 80}\nSTATS\n{df['len_1'].describe()}\n{'*' * 80}\n")

    # apply filter with scaffold length and connections size
    df = df[(df["len_1"] >= prog_args.fragments_min) & (df[f"{sample1}_len"] >= prog_args.scaff_min)]

    if prog_args.plot:
        font_size = 18
        chunk_size = 80
        fig, axes = plt.subplots(figsize=(45, 20), dpi=150, nrows=1, ncols=1)
        plt.subplots_adjust(top=0.95)
        fig.suptitle(f'Graph {sample1} VS {sample2}', fontsize=font_size + 6)

        df.hist(column='len_1', grid=True, xrot=90, ax=axes, figsize=None, bins=5000)

        plt.xticks(np.arange(1, max(df['len_1']), int(max(df['len_1']) / chunk_size)),
                   rotation=70)  # ajoute la grille an 80 partie
        # axes.set_xlim(1, 1000)
        plt.savefig(f"Plot_{sample1}-{sample2}.png", bbox_inches='tight', pad_inches=1)

    # write df to file filter
    with open(f"CIRCA_{sample1}-{sample2}_filter.csv", "w") as filter_file:
        df.to_csv(filter_file, index=False, sep=sep, header=True)

    # build scaffold size file for CIRCA
    df_ch1 = df[[f"{sample1}-chrom", f"{sample1}_len"]].drop_duplicates(ignore_index=True).sort_values(
            by=f"{sample1}_len",
            ascending=False)
    if prog_args.debug: print(df_ch1)
    df_ch2 = df[[f"{sample2}-chrom", f"{sample2}_len"]].drop_duplicates(ignore_index=True).sort_values(
            by=f"{sample2}_len",
            ascending=False)
    if prog_args.debug: print(df_ch2)

    # write to file
    with open(f"CIRCA_{sample1}-{sample2}_scaffold_len.csv", "w") as chrom_len_file:
        chrom_len_file.write(f"Scaffold{sep}size\n")
        df_ch2[df_ch2[f"{sample2}_len"] >= prog_args.scaff_min].to_csv(chrom_len_file, index=False, sep=sep,
                                                                       header=False)
        df_ch1[df_ch1[f"{sample1}_len"] >= prog_args.scaff_min].to_csv(chrom_len_file, index=False, sep=sep,
                                                                       header=False)


if __name__ == '__main__':
    main()
