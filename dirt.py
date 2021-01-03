#! /nv/hp10/adas30/bin/python

import csv
import os
import time
from os.path import join

import yaml
from math import ceil

import click
import imageio
import numpy as np

import thresholding
from options import DIRTOptions
from orientation import fix_orientation
from segmentation import segment
from traits import traits
from utils import print_cli_header


@click.command()
@click.argument('input_file')
@click.option('--output_directory', required=False, type=str, default='')
@click.option('--excised_roots', required=False, type=int, default=0)
@click.option('--marker_diameter', required=False, type=float, default=25.4)
@click.option('--stem_reconstruction', required=False, type=bool, default=False)
@click.option('--plot', required=False, type=bool, default=True)
def cli(input_file,
        output_directory,
        excised_roots,
        marker_diameter,
        stem_reconstruction,
        plot):
    print_cli_header()

    start = time.time()
    options = DIRTOptions(
        input_file=input_file,
        output_directory=output_directory,
        excised_roots=excised_roots,
        marker_diameter=marker_diameter,
        stem_reconstruction=stem_reconstruction,
        plot=plot)

    # fix orientation of the image in TIFF and JPG files
    fix_orientation(options.input_file, replace=True)
    output_prefix = join(options.output_directory, options.input_stem)
    image = imageio.imread(options.input_file, as_gray=True)
    if len(image) == 0:
        raise ValueError(f"Image is empty: {options.input_name}")

    # apply binary threshold
    masked = thresholding.binary_threshold(image.astype(np.uint8))

    # segment crown, marker, tag, and excised root(s)
    results = segment(options, masked)

    # extract traits
    results = traits(options, masked, results)

    print(f"All done in just {ceil((time.time() - start))} seconds! Writing output file")

    with open(f"{output_prefix}.results.yml", 'w') as file:
        yaml.dump(results, file, default_flow_style=False)

    # with open('output.csv', 'a' if os.path.isfile('output.csv') else 'w', newline='') as output_file:
    #     writer = csv.writer(output_file, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # TODO write header row
        # TODO write output rows


if __name__ == '__main__':
    cli()
