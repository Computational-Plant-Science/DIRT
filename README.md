# DIRT 

![CI](https://github.com/Computational-Plant-Science/DIRT/workflows/CI/badge.svg) [![Coverage Status](https://coveralls.io/repos/github/Computational-Plant-Science/DIRT/badge.svg?branch=plantit)](https://coveralls.io/github/Computational-Plant-Science/DIRT?branch=master)

An automatic root phenotyping tool.

## Requirements

- [Python3](https://www.python.org)
- Docker

### Python dependencies

The software depends on:
- [`graph-tool`](http://graph-tool.skewed.de)
- [`mahotas`](http://luispedro.org/software/mahotas)
- [`numpy`](http://sourceforge.net/projects/numpy/)
- [`scipy`](http://www.scipy.org/SciPy)
- [`matplotlib`](https://matplotlib.org/)
- [`scikit-image`](https://scikit-image.org/)
- [`imageio`](https://imageio.github.io/)
- [`click`](https://click.palletsprojects.com/en/7.x/)
- [`pandas`](https://pandas.pydata.org/)

### Binary dependencies

- [tesseract](https://code.google.com/p/tesseract-ocr/)
- [zbar](http://zbar.sourceforge.net)

## Installation

Clone the repo with `git clone `.

## Usage

A good way to get a feel for DIRT is to run the test cases:

```bash
docker run -it -v $PWD:/opt/dev -w /opt/dev computationalplantscience/dirt pytest -s
```

At its simplest, DIRT can be invoked with:

```bash
python dirt.py <image file>
```

For instance, `dirt.py root1.png` or `dirt.py root2.jpg`. Supported filetypes are JPEG, PNG, and TIFF.

You can also specify an `--output_directory` to write files to.

Note that it is not possible to analyze only an excised root when a root crown is in the image. However, it is possible to analyze compute images containing only excised roots.

## References

This software uses, modifies, or references:
- http://www.daniweb.com/software-development/python/threads/31449/k-means-clustering
- Adaptive thresholding code from [`scikit-image`](http://scikit-image.org)
- Orientation correction code from [`python-image-orientation-patch`](https://github.com/kylefox/python-image-orientation-patch)

## Attribution

Please cite the DIRT publication if you use DIRT for your project.:

Bucksch et al., 2014 "Image-based high-throughput field phenotyping of crop roots", Plant Physiology

## Author

Alexander Bucksch
Department of Plant Biology
Warnell School of Forestry and Natural Resources
Institute of Bioinformatics
University of Georgia
[Email](mailto:bucksch@uga.edu)
[Website](http://www.computational-plant-science.org)

## Changelog

### January 2021

CLI refactored and tests added anticipating installation on PlantIT.

### 21 June 2019

Some bug fixes on the avg. root density. There was a problem with very young and sparse root system. The formula changed and is now normed to the max. width instead of the max. width of the line.
The bug was found by Peng Wang at the University of Nebraska.

### 11 January 2016

Minor bug fixes in Preprocessing.py to allow smaller circle markers and fix a possible missdetection of the experiment tag as the circle. 
Thanks to Linda Zamariola (U Bologna) for finding this issue.

### 4 November 2015

Minor bug fixes in the excised root calculations. Thanks to Alexandre Grondin (U Nebraska) for discovering and validating the fixes.

### 14 January 2015

- storage of trait values is changed from a list data structure to a dictionary to allow trait selection controlled by the file traits.csv
- added support for trait selection to reduce computation time. See example file traits.csv (1 - trait is computed, 0 - trait is not computed)
- removed unused tip-diameter switch on the command line
- add stem reconstruction switch on the command line to turn the experimental stem reconstruction on/off
- output file now uses the codes in the trait.csv file and only contains selected traits
- removed several unused variables and minor bugs fixed
- added command line option to turn storage of numpy arrays on/off. These files can be used to plot the individual root statistics and can be found in the "Plots" folders.
- new (experimental, not validated) traits added due to community requests: projected root area, width and depth of the skeleton (medial axis), top and bottom angle for monocots, segmentation of adventious and basal roots for legumes to retrieve taproot and hypocotyl diameter and adventious and basal root counts.
- added computational statistics such as computation time and graph size to help balancing grid installations
- added an option to have an output file with all possible traits that contains empty cells for not computed traits in the output.csv file. This was a developer request to enable faster ingestion into data bases
