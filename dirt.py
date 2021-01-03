#! /nv/hp10/adas30/bin/python

import csv
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
from results import DIRTResults
from segmentation import segment
from traits import traits
from utils import print_cli_header


@click.command()
@click.argument('input_file')
@click.option('--output_directory', required=False, type=str, default='')
@click.option('--marker_diameter', required=False, type=float, default=25.4)
def cli(input_file, output_directory, marker_diameter):
    print_cli_header()
    start = time.time()
    options = DIRTOptions(
        input_file=input_file,
        output_directory=output_directory,
        marker_diameter=marker_diameter)

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

    # compute total runtime
    duration = ceil((time.time() - start))
    print(f"All done in just {duration} seconds! Writing output file")
    results = {**results, **DIRTResults(duration=duration)}

    # write YAML output file
    with open(f"{output_prefix}.results.yml", 'w') as file:
        yaml.dump(results, file, default_flow_style=False)

    # write CSV output file
    with open(f"{output_prefix}.results.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(list(results.keys()))
        writer.writerow(list(results.values()))


if __name__ == '__main__':
    cli()
