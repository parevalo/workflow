#!/usr/bin/env python
import click
import numpy as np
import rasterio

@clicl.command()
@click.argument('input', metavar='<input raster>', nargs=1, type=click.Path(exists=True, resolve_path=True)
@click.argument('output', metavar='<output csv>', nargs=1, type=click.Path(resolve_path=True)
def count_pixels(input, output):
    """Script to count the number of pixels per class in a given raster file. 
   Outputs a csv with the class numbers and pixel count.
   Current version only checks the first band."""
    
    # Read the data
    raster = rasterio.open(input)
    data = raster.read(1)
    
    # Count unique values 
    flatraster = data.flatten()
    stats = np.bincount(flatraster)
    max = flatraster.max()
    min = flatraster.min()
    outarray = np.zeros((max+1, 2))
    outarray[:, 0] = np.arange(max+1)
    outarray[:, 1] = stats


    # Save csv with classes and pixels per class
    np.savetxt(output, outarray, delimiter=',', header="stratum,pixels", fmt='%01d')

if __name__ == '__main__':
    count_pixels()

