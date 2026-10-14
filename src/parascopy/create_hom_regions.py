#!/usr/bin/env python3

import argparse
import sys
import gzip
import os


def _parse_region(region_str):
    """Parse region string like 'chr1:1000-2000' into (chrom, start, end)."""
    parts = region_str.split(':')
    chrom = parts[0]
    coord_parts = parts[1].split('-')
    start = int(coord_parts[0]) - 1
    end = int(coord_parts[1]) - 1
    return (chrom, start, end)


def _sort_regions(regions):
    """Sort regions by chromosome and position."""
    def chrom_sort_key(region):
        chrom = region[0]
        chrom_num_str = chrom.replace("chr", "").replace("X", "30").replace("Y", "31")
        try:
            chrom_num = int(chrom_num_str)
        except ValueError:
            chrom_num = 999
        return (chrom_num, region[1])

    return sorted(regions, key=chrom_sort_key)


def _get_model_files(input_paths):
    """Expand input paths to list of .gz files."""
    filenames = []
    for path in input_paths:
        if os.path.isdir(path):
            filenames.extend([os.path.join(path, f) for f in os.listdir(path) if f.endswith('.gz')])
        elif path.endswith('.gz'):
            filenames.append(path)
        else:
            with open(path) as inp:
                curr_dir = os.path.dirname(path)
                filenames.extend([os.path.join(curr_dir, line.strip()) for line in inp])
    return filenames


def create_hom_regions_bed(model_paths, output_file):
    """
    Extract homologous regions from model files and write to BED file.

    Args:
        model_paths: List of directories or files containing model files (.gz)
        output_file: Output BED file path
    """
    model_files = _get_model_files(model_paths)
    regions = []

    for filename in model_files:
        try:
            with gzip.open(filename, 'rt') as f:
                first_line = f.readline().strip()
                region_str = first_line.split()[1]
                region_tuple = _parse_region(region_str)
                regions.append(region_tuple)
        except (IndexError, ValueError) as e:
            sys.stderr.write(f'Warning: Cannot parse region from "{filename}": {e}\n')

    sorted_regions = _sort_regions(regions)

    with open(output_file, 'w') as out:
        for r in sorted_regions:
            out.write("\t".join((r[0], str(r[1]), str(r[2]))))
            out.write("\n")


def main(prog_name, args):
    """Command-line interface for create_hom_regions."""
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description='Create homologous regions BED file from parascopy model files.'
    )

    parser.add_argument(
        'model',
        nargs='+',
        metavar='<model>',
        help='Path(s) to model directory or model file(s). Can be a directory containing .gz files, '
             'a single .gz file, or a text file with list of .gz files (one per line).'
    )

    parser.add_argument(
        '-o', '--output',
        metavar='<file>',
        required=True,
        help='Output BED file'
    )

    parsed_args = parser.parse_args(args)

    create_hom_regions_bed(parsed_args.model, parsed_args.output)


if __name__ == '__main__':
    main('parascopy create-hom-regions', sys.argv[1:])
