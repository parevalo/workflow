#!/usr/bin/env python
import click
import fiona
from shapely.geometry import mapping, shape
from fiona import collection
from functools import partial
import pyproj
import sys


@click.command()
@click.argument('sample_shp', metavar='<input sample shapefile>', nargs=1,
                type=click.Path(exists=True, resolve_path=True))
@click.argument('footprints_shp', metavar='<landsat footprints>', nargs=1,
                type=click.Path(exists=True, resolve_path=True))
@click.argument('output_shp', metavar='<output shapefile>', nargs=1, type=click.STRING)
def extract_pathrow(sample_shp, footprints_shp, output_shp):
    """Intersect samples with path-rows and extract their info.
    Assumes both are in the same coordinate system
    If samples intersect the boundary of two scenes, duplicate samples will be created"""

    # Open shapefile and print CRS
    samples = fiona.open(sample_shp, "r")
    click.echo("Samples' coordinate system is: {}".format(samples.crs['init']))
    footprints = fiona.open(footprints_shp, "r")
    click.echo("Footprints' coordinate system is: {}".format(footprints.crs['init']))

    # TODO: Check coordinate systems and reproject if necessary
    # project = partial(
    #     pyproj.transform,
    #     pyproj.Proj(init='epsg:32618'),
    #     pyproj.Proj(init='epsg:32619'))
    # # If we wanted to reproject
    # transform(project, samples(test[1]['geometry']))

    # Set output schema and CRS
    schema = {'geometry': 'Polygon', 'properties': {'ID': 'int:10', 'PTRW': 'int:10'}}
    if samples.crs['init'] == footprints.crs['init']:
        crs = samples.crs
    else:
        sys.exit("Samples and footprints are not in the same coordinate system. Aborted")

    # Sorted list of scenes, in order to assign them sequentially
    scene_list = 358, 359, 457, 458, 459, 461, 462, 557, 558, 559, 560, 561, 658, 659, 660, 661, 758, 759, 760, 761, \
                 858, 859, 860, 959, 960

    # Iterate over points to check if they intersect with east or west zone, and write shapes
    # to do it in order, create a list with the order of the scenes and then iterate over that list, calling
    # the polygon by the properties.

    with collection(output_shp, "w", "ESRI Shapefile", schema, crs=crs) as output_west:
        for s in scene_list:
            scene = [scn for scn in footprints if scn['properties']['PTRW'] == s][0]
            #pr2=transform(project, pr)
            for point in samples:
                if shape(scene['geometry']).intersects(shape(point['geometry'])):
                    output_west.write({
                        'properties': {
                            'ID': point['properties']['ID'],
                            'PTRW': scene['properties']['PTRW']
                        },
                        'geometry': point['geometry']
                    })


if __name__ == '__main__':
    extract_pathrow()
