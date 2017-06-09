#!/usr/bin/env python
import click
import os
import fiona
import rasterio
import shapely
from shapely.geometry import mapping, shape
from fiona import collection

@click.command()
@click.argument('sample_shp', metavar='<input sample shapefile>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('east_shp', metavar='<input east zone>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('west_shp', metavar='<input west zone>', nargs=1, type=click.Path(exists=True, resolve_path=True))
@click.argument('east_sample_fname', metavar='<output east sample filename>', nargs=1, type=click.STRING)
@click.argument('west_sample_fname', metavar='<output west sample filename>', nargs=1, type=click.STRING)
def intersect_and_split(sample_shp, east_shp, west_shp, east_sample_fname, west_sample_fname):
    """Intersect and split samples in the east or west zones and create new separate files"""

    # Open shapefiles
    east = fiona.open(east_shp, "r")
    west = fiona.open(west_shp, "r")
    samples = fiona.open(sample_shp, "r")

    # Set output schema
    schema = {'geometry': 'Polygon', 'properties': {'ID': 'int:10'}}
    crs = {'init': 'epsg:32618'}

    # Iterate over points to check if they intersect with east or west zone, and write shapes
    with collection(east_sample_fname, "w", "ESRI Shapefile", schema, crs=crs) as output_east:
        for point in samples:
           if shape(east[0]['geometry']).intersects(shape(point['geometry'])):
                output_east.write({
                    'properties': {
                        'ID': point['properties']['ID']
                    },
                    'geometry': point['geometry']
                })

    with collection(west_sample_fname, "w", "ESRI Shapefile", schema, crs=crs) as output_west:
        for point in samples:
           if shape(west[0]['geometry']).intersects(shape(point['geometry'])):
                output_west.write({
                    'properties': {
                         'ID': point['properties']['ID']
                    },
                    'geometry': point['geometry']
                })

    # Print matched records per zone and total
    click.echo("Number of samples, east zone: {}".format(len(output_east)), err=True)
    click.echo("Number of samples, west zone: {}".format(len(output_west)), err=True)
    click.echo("Total number of samples: {}".format(len(output_east) + len(output_west)), err=True)

if __name__ == '__main__':
    intersect_and_split()

