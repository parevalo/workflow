#!/usr/bin/env python
import click
import numpy as np
import rasterio
import csv


@click.command()
@click.argument('raster1', metavar='<raster 1>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('raster2', metavar='<raster 2>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('csvfile', metavar='<lut>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('output', metavar='<output raster>', nargs=1, type=click.Path(resolve_path=True))
@click.option('--header', metavar='', type=bool, default=True, show_default=True, help='Does CSV have header?')
@click.option('--delim', metavar='', type=str, default=',', show_default=True, help='CSV Text delimiter')
@click.option('--quote', metavar='', type=str, default='"', show_default=True, help='CSV quote character')
@click.option('--format', metavar='', type=str, default='GTiff', show_default=True, help='Output file format')
def create_strata(raster1, raster2, csvfile, header, output, delim, quote, format):
    """ Script to create a stratification raster from two maps
    using a lookup table"""

    ds1 = rasterio.open(raster1)
    ds2 = rasterio.open(raster2)

    with open(csvfile) as f:
        csvreader = csv.reader(f, delimiter=delim, quotechar=quote)
        if header:
            next(csvreader)
        lut_file = {(int(row[0]), int(row[1])): int(row[2]) for row in csvreader}

    # Calculate output dtype based on LUT values
    out_dt = 'byte'
    if min(lut_file.values()) < 0:
        # Must be signed int
        if max(np.abs(lut_file.values())) < 2 ** 15:
            out_dt = 'int8'
        elif max(np.abs(lut_file.values())) < 2 ** 31:
            out_dt = 'int32'
        elif max(np.abs(lut_file.values())) < 2 ** 63:
            out_dt = 'int64'
        else:
            click.echo('Required output data type is unknown')
            click.Abort()

    else:
        # Can be unsigned
        if max(lut_file.values()) < 2 ** 8:
            out_dt = 'uint8'
        elif max(lut_file.values()) < 2 ** 16:
            out_dt = 'uint16'
        elif max(lut_file.values()) < 2 ** 32:
            out_dt = 'uint32'
        elif max(lut_file.values()) < 2 ** 64:
            out_dt = 'uint64'
        else:
            click.echo('Required output data type is unknown')
            click.Abort()

    # Init output
    y1 = ds1.read(1)
    y2 = ds2.read(1)
    profile = ds1.profile
    profile.update(driver=format, dtype=out_dt)
    map = np.zeros_like(y1, dtype=out_dt)

    # Match codes and create map.
    for k, v in lut_file.items():
        m1 = np.in1d(y1, k[0]).reshape(y1.shape)
        m2 = np.in1d(y2, k[1]).reshape(y2.shape)
        m_1_2 = np.logical_and(m1, m2)

        #  Add to map
        map[m_1_2] = v

    # Save using the profile taken from the input file directly. If needed, options
    # could be specified manually or by updating the profile itself. Map.astype
    # required to match the types and avoid writing larger files than needed

    with rasterio.open(output, mode='w',  **profile) as dst:
        dst.write(map.astype(out_dt), 1)


if __name__ == '__main__':
    create_strata()
