#!/usr/bin/env python
import click
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import rasterio
import fiona
import os
import csv


@click.command()
@click.argument('raster1', metavar='<raster 1>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('raster2', metavar='<raster 2>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('lut', metavar='<lut>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('output', metavar='<output raster>', nargs=1, type=click.Path(resolve_path=True))
@click.option('--delim', metavar='', type=str, default=',', show_default=True, help='Text delimiter')
def create_strata(raster1, raster2, lut, output, delim):
    """ Script to create a stratification raster from two maps
    using a lookup table"""

    ds1 = rasterio.open(raster1)
    ds2 = rasterio.open(raster2)

    lut_file = np.loadtxt(lut, delimiter=delim, skiprows=1, dtype='uint8')

    # TODO add code from Chris to calculate dtype automatically

    # Init output
    y1 = ds1.read(1)
    y2 = ds2.read(1)
    profile = ds1.profile
    map = np.zeros_like(y1)

    # Match codes and create map. UPDATE LUT TO INCLUDE NODATA TRANSITIONS!
    for row in lut_file:
        m1 = np.in1d(y1, row[0]).reshape(y1.shape)
        m2 = np.in1d(y2, row[1]).reshape(y2.shape)
        m_1_2 = np.logical_and(m1, m2)

        #  Add to map
        map[m_1_2] = row[2]

    # Save using the profile taken from the input file directly. If needed, options
    # could be specified manually or by updating the profile itself. Map.astype
    # required to match the types and avoid writing larger files than needed

    with rasterio.open(output, mode='w', **profile) as dst:
        dst.write(map.astype('uint8'), 1)

#
if __name__ == '__main__':
    create_strata()