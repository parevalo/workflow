import os
import numpy as np
import fiona
import rasterio
import shapely
from shapely.geometry import mapping, shape
from shapely.ops import cascaded_union
from fiona import collection
from collections import defaultdict


os.chdir("/media/paulo/785044BD504483BA/test/")

with rasterio.open('2008_test_window.tif') as src:
    ds1 = src.read(1)
    out_meta = src.meta.copy()
    transform = out_meta['transform']

# Polygonize and store in list as tuples
polygonized = [(shapely.geometry.shape(g), value) for g, value in
               rasterio.features.shapes(ds1, transform=transform)]


# Or polygonize and store in dictionary. PAY ATTENTION TO USING THE TRANSFORM TO GET THE RIGHT COORDINATES!
polys_by_label = defaultdict(list)
for geom, value in rasterio.features.shapes(ds1, transform=transform):
    polys_by_label[value].append(shapely.geometry.shape(geom))


# Dissolve based on keys and store in regular dict, bc defauldict doesn't work when saving.
patches_by_label = {}
for keys in polys_by_label:
    patches_by_label[keys] = cascaded_union(polys_by_label[keys])


# Define schema manually, i.e. type of vector and field characteristics
schema = {'geometry': 'Polygon', 'properties': {'id': 'int:10'}}
crs = {'init': 'epsg:32618'}
# Iterate over features and write new file. If file exists it will be overwritten.
with collection("test_write_shape.shp", "w", driver='ESRI Shapefile', schema=schema, crs=crs) as output:
    # This loop doesn't work if it's a defaultdict!!
    for key in patches_by_label:
        output.write({
            'properties': {
                 'id': key
            },
            'geometry': mapping(shape(patches_by_label[key]))
        })
