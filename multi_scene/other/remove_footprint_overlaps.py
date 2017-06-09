#!/usr/bin/env python
import fiona
from fiona import collection
from shapely.geometry import mapping, shape
from shapely.ops import cascaded_union
import os

# Script to properly assing the overlapping areas of the scene footprints to the correct scene, in order to
# make it easier to assign the path-row value to the biannual reference samples

os.chdir("/media/paulo/785044BD504483BA/OneDrive/Lab/Base layers")

# Define schema and output crs
schema = {'geometry': 'Polygon', 'properties': {'ID': 'int:10', 'PTRW': 'int:10'}}
crs = {'init': 'epsg:32618'}

# Sorted list of scenes, in order to process them sequentially
scene_list = 358, 359, 457, 458, 459, 461, 462, 557, 558, 559, 560, 561, 658, 659, 660, 661, 758, 759, 760, 761, \
            858, 859, 860, 959, 960
rev_scene_list = scene_list[::-1]

# Create shapefile of scenes without the overlaps, in the right order
with collection("WRS2_amazon_selection_UTM18_no-overlap.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for i in range(len(scene_list)):
        with fiona.open("WRS2_amazon_selection_UTM18.shp", "r") as footprints:
            # Get scene in order (matching the scene list)
            scene = [scn for scn in footprints if scn['properties']['PTRW'] == scene_list[i]][0]
            # Get the rest of the scenes left and merge them
            rest = [shape(scn['geometry']) for scn in footprints if scn['properties']['PTRW'] in rev_scene_list[:-i-1]]
            rest_merged = cascaded_union(rest)
            # Find the differences between them and write that to file
            a = shape(scene['geometry']).difference(rest_merged)
            output.write({
                'properties': {
                    'ID': scene['id'],
                    'PTRW': scene['properties']['PTRW']
                 },
                'geometry': mapping(a)
             })
