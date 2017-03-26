#!/usr/bin/env python
import click
import numpy as np
import rasterio
from rasterio.plot import show_hist
import matplotlib.pyplot as plt

@click.command()
@click.argument('input', metavar='<input raster>', nargs=1, type=click.Path(exists=True, resolve_path=True))
def raster_histogram(input):
    """Create raster histograms each band in a raster and save them to file"""
    
    # Read the data
    raster = rasterio.open(input)

    for i in raster.indexes:
        fig, ax = plt.subplots(1)
        show_hist(raster.read(i, masked=True), bins=50, lw=2, masked=True, alpha=0.3,
                  title="Histogram - band{}".format(i), ax=ax, color='blue')
        fig.savefig("hist_band{}".format(i), dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    raster_histogram()
