#!/usr/bin/env python
import click
import numpy as np
import rasterio
from rasterio.plot import show_hist
import matplotlib.pyplot as plt


@click.command()
@click.argument('input', metavar='<input raster>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('varname', metavar='<variable name>', default='Value',  type=click.STRING)
def raster_histogram(input, varname):
    """Create raster histograms each band in a raster and save them to file"""
    
    # Read the data
    raster = rasterio.open(input)

    for i in raster.indexes:
        fig, ax = plt.subplots(1)
        # TODO: Make this calculation automatic based on SD or something
        ax.set_xlim([-0.5, 0.5])
        show_hist(raster.read(i, masked=True), bins=50, lw=2, masked=True, alpha=0.6,
                  title="Histogram - band{}".format(i), ax=ax, facecolor='blue')
        plt.xlabel(varname)
        fig.savefig(varname + "_hist_band{}".format(i), dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    raster_histogram()
