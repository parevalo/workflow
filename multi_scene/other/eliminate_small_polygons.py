
# Requires specifying PYTHONPATH to be that of qgis, in this case
# export PYTHONPATH=/share/pkg/qgis/2.6.1/install/share/qgis/python
# and loading qgis and python/2.7.5 modules
# Input layers must have a projected system in m

#TODO
# Allow passing arguments from terminal
# Properly check if field exists
# Check if output file exists

import sys
from qgis.core import *
from PyQt4.QtCore import *

import os
import urllib
import zipfile
import tempfile


# Initialize QGIS Application
QgsApplication.setPrefixPath("/share/pkg/qgis/2.6.1/install/bin/qgis", True)
qgs =QgsApplication([], False)
qgs.initQgis()

# Add the path to Processing framework
sys.path.append('/share/pkg/qgis/2.6.1/install/share/qgis/python/plugins')

# Import and initialize Processing framework
from processing.core.Processing import Processing
Processing.initialize()
import processing

# Load vector
layerpath="C:\OneDrive\Lab\Base layers\scene_overlap.shp"
layerpath=" /projectnb/landsat/projects/Colombia/vector/scene_overlap.shp"

# To extract basename from the path if needed (THIS IS NOR WORKING CORRECTLY)
path = os.path.splitext(layerpath)[0]
filename=os.path.basename(path)

layer = QgsVectorLayer(layerpath, 'scene_overlap', "ogr")
if not layer.isValid():
  print "Layer failed to load!"

# Get layer data provider (starting edition session may be required after this)
dp = layer.dataProvider()

# Check existing fields and create a new one to store area if possible, FIX
for field in layer.fields():
    if field.name() == "area_ha":
        print "There's an existing field with the name area_ha, overwritting"

# Create field and populate it
dp.addAttributes([QgsField("area_ha", QVariant.Double)])

 Check geometry and break if there are errors
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

# Calculate areas and save them in the area_ha column
area = 0
with edit(layer):
    for gFeat in layer.getFeatures():
        calculator = QgsDistanceArea()
        calculator.setEllipsoid('WGS84')
        calculator.setEllipsoidalMode(True)
        calculator.computeAreaInit()
        geom = gFeat.geometry()
        #landArea = gFeat['area_ha']
        
        #Check if polygon is multipart 
        if geom.isMultipart() is False: 
            polyg = geom.asPolygon() #transform to list of points
            if len(polyg) > 0:
                area = calculator.measurePolygon(polyg[0])
                # Convert sq km to ha
                gFeat['area_ha'] = area *0.0001
                #landArea = area *0.0001
                layer.updateFeature(gFeat)
        else: 
            multi = geom.asMultiPolygon()
            for polyg in multi:
                area = area + calculator.measurePolygon(polyg[0])
                gFeat['area_ha'] = area *0.0001
                #landArea = area * 0.0001
                print "Calculated on a multipart polygon"

# This is how to do selection but it's not required because the tool includes it
# Set expression and select based on it
#expr = QgsExpression( 'area_ha <= 14450754076.222845077514648' )
#it = layer.getFeatures( QgsFeatureRequest( expr ) )
# Build list of features
#ids = [i.id() for i in it]
# Select features with the ids
#layer.setSelectedFeatures( ids )
# Remove selection
#layer.setSelectedFeatures([])

# Check if output layer exists and warn 


# Run the eliminate sliver polygons
minarea = '15000' #15K ha in this example
processing.runalg("qgis:eliminatesliverpolygons", layer, 'False', 'area_ha', 5,minarea, 0, 'overlap_merged.shp')

print "Done!"
qgs.exitQgis()

