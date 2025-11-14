"""
Module for managing the CLI layer on top of the API; also contains the CLI entry point under `main()`.
"""
from argparse import RawTextHelpFormatter, RawDescriptionHelpFormatter, ArgumentParser
from sys import stdout
from pathlib import Path
from typing import get_args
from pathogenx import RESOURCES
from pathogenx.utils import bold
from pathogenx.io import _GENOTYPE_FLAVOURS, _META_FLAVOURS, _DIST_FLAVOURS

# Constants ------------------------------------------------------------------------------------------------------------
_LOGO = (f"\033[1;35m========================|> PathoGenX |>========================\n"
         f"{'A Python library for Pathogen Genotype eXploration':^63}\033[0m")


# Functions ------------------------------------------------------------------------------------------------------------
def prevalence_parser(subparsers):
    name, desc = 'prevalence', 'Calculate prevalence in a dataset'
    parser = subparsers.add_parser(
        name, description=_LOGO, prog=f'{RESOURCES.package} {name}',
        formatter_class=RawTextHelpFormatter, help=desc,
        usage="%(prog)s <genotype> [metadata] [distance] [options]", add_help=False
    )
    inputs = parser.add_argument_group(bold('Inputs'), '')
    inputs.add_argument('genotypes', metavar='<genotypes>', help='Genotype file', type=Path)
    inputs.add_argument('metadata', metavar='<metadata>', help='Optional metadata file', nargs='?', type=Path)
    inputs.add_argument('distances', metavar='<distances>', help='Optional distance file', nargs='?', type=Path)
    inputs.add_argument('--genotype-flavour', help='Genotype file flavour (default: %(default)s)\n'
                                                   '(choices: %(choices)s)', metavar='',
                        choices=get_args(_GENOTYPE_FLAVOURS), default='pw-kleborate')
    inputs.add_argument('--metadata-flavour', help='Metadata file flavour (default: %(default)s)\n'
                                                   '(choices: %(choices)s)', metavar='',
                        choices=get_args(_META_FLAVOURS), default='pw-metadata')
    inputs.add_argument('--distance-flavour', help='Distance file flavour (default: %(default)s)\n'
                                                   '(choices: %(choices)s)', metavar='',
                        choices=get_args(_DIST_FLAVOURS), default='pw-dist')

    calc = parser.add_argument_group(bold('Calculator options'), '')
    calc.add_argument('--stratify-by', help='List of columns to stratify the analysis by', nargs='+', metavar='')
    calc.add_argument('--adjust-for', help='Optional list of columns for adjustment (e.g., Cluster)', nargs='*', metavar='')
    calc.add_argument('--n-distinct', help='Optional list of columns to calculate distinct counts for', nargs='*', metavar='')
    calc.add_argument('--denominator', help='Optional column to use as the primary grouping for denominators\n'
                                            'If None, the first column in stratify-by is used', metavar='')

    calc = parser.add_argument_group(bold('Clustering options'), '')
    calc.add_argument('--snp-distance', type=int, default=20, metavar='',
                      help='The maximum distance for two samples to be considered connected.\n'
                           'Only used when `method` is connected_components')

    opts = parser.add_argument_group(bold('Other options'), '')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')


# Main CLI Entry Point -------------------------------------------------------------------------------------------------
def main():
    parser = ArgumentParser(
        description=_LOGO,
        usage="%(prog)s <command>", add_help=False, prog=RESOURCES.package, formatter_class=RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(
        title=bold('Command'), dest='command', metavar='<command>', required=True, help=None,
    )
    prevalence_parser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), '')
    opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
    opts.add_argument('-h', '--help', help='Show this help message and exit', action='help')

    args = parser.parse_args()

    if args.command == 'prevalence':
        from pathogenx.io import GenotypeFile, MetaFile, DistFile
        from pathogenx.dataset import Dataset
        from pathogenx.calculators import PrevalenceCalculator
        metadata_file, distance_file = None, None
        genotype_file = GenotypeFile.from_flavour(args.genotypes, args.genotype_flavour)
        if args.metadata is not None:
            metadata_file = MetaFile.from_flavour(args.metadata, args.metadata_flavour)
        if args.distances is not None:
            distance_file = DistFile.from_flavour(args.distances, args.distance_flavour)

        dataset = Dataset.from_files(genotype_file, metadata_file, distance_file)
        if dataset.distances is not None:
            dataset.calculate_clusters(distance=args.snp_distance)
        calculator = PrevalenceCalculator(args.stratify_by, args.adjust_for, args.n_distinct, args.denominator)
        result = calculator.calculate(dataset)
        result.data.to_csv(stdout, sep='\t', index=False)
