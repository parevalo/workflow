# Requires specifying PYTHONPATH to be that of qgis, in this case
# export PYTHONPATH=/share/pkg/qgis/2.6.1/install/share/qgis/python
# and loading qgis and python/2.7.5 modules

import sys
from qgis.core import *
from PyQt4.QtCore import *

import os
import urllib
import zipfile
import tempfile

# Initialize QGIS Application
#QgsApplication.setPrefixPath("/share/pkg/qgis/2.6.1/install/bin/qgis", True)
#app = QgsApplication([], True)
#QgsApplication.initQgis()

# Add the path to Processing framework
#sys.path.append('/share/pkg/qgis/2.6.1/install/share/qgis/python/plugins')

# Import and initialize Processing framework
from processing.core.Processing import Processing
Processing.initialize()
import processing

# Load vector
layerpath="C:\OneDrive\Lab\Base layers\scene_overlap.shp"

# To extract basename from the path if needed
path = os.path.splitext(layerpath)[0]
filename=os.path.basename(path)

layer = QgsVectorLayer(layerpath, filename, "ogr")
if not layer.isValid():
  print "Layer failed to load!"

# Get layer data provider (starting edition session may be required after this)
dp = layer.dataProvider()

# Check existing fields and create a new one to store area if possible, FIX
for field in layer.fields():
    if field.name() == "area_ha":
        print "There's an existing field with the name area_ha, overwritting"

# Create field and populate it
#dp.addAttributes([QgsField("area_ha", QVariant.Double)])

#TODO Add geometry check to avoid failures due to bad geoms
for feat in layer.getFeatures():
    geom = feat.geometry()
    if geom:
        err = geom.validateGeometry()
        if not err:
            print "No geometry errors found"
            continue
        else:
            print '%d geometry errors detected (feature %d)' % (len(err), feature.id()) 
            break

# Calculate areas, currently calculating them but not saving them
area = 0
with edit(layer):
    for gFeat in layer.getFeatures():
        calculator = QgsDistanceArea()
        calculator.setEllipsoid('WGS84')
        calculator.setEllipsoidalMode(True)
        # here you need to put the code
        calculator.computeAreaInit()
        geom = gFeat.geometry()
        #landArea = gFeat['area_ha']
        #Check if polygon is multipart 
        if geom.isMultipart() is False: 
            polyg = geom.asPolygon() #transform to list of points
            if len(polyg) > 0:
                area = calculator.measurePolygon(polyg[0])
                # Convert sq km to ha
                gFeat['area_ha'] = area / 100000000
                #landArea = area /1000000
                layer.updateFeature(gFeat)
        else: 
            multi = geom.asMultiPolygon()
            for polyg in multi:
                area = area + calculator.measurePolygon(polyg[0])
                l#andArea = area * 100
                print "Calculated on a multipart polygon"

print "Done!"



