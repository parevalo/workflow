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
QgsApplication.setPrefixPath("/share/pkg/qgis/2.6.1/install/bin/qgis", True)
app = QgsApplication([], True)
QgsApplication.initQgis()

# Add the path to Processing framework
sys.path.append('/share/pkg/qgis/2.6.1/install/share/qgis/python/plugins')

# Import and initialize Processing framework
from processing.core.Processing import Processing
Processing.initialize()
import processing

# Load vector
layer = QgsVectorLayer(poly_out, layer_name, "ogr")
if not layer.isValid():
  print "Layer failed to load!"

# To extract basename from the path if needed
path = os.path.splitext("C:\OneDrive\Lab\Base layers\scene_overlap.shp")[0]
os.path.basename(path)

# Get layer data provider (starting edition session may be required after this)
dp = layer.dataProvider()

# Check existing fields and create a new one to store area if possible
for field in layer.fields():
    if field.name() == "area_ha":
        print "There's an existing field with the name area_ha, overwritting"
    else:
        # Create field and populate it
        prov.addAttributes([QgsField("area_ha", QVariant.Double)])
        area = 0
        with edit(layer):
            for gFeat in layer.getFeatures():
                calculator = QgsDistanceArea()
                calculator.setEllipsoid('WGS84')
                calculator.setEllipsoidalMode(True)

                # here you need to put the code
                calculator.computeAreaInit()

                geom = gFeat.geometry()
                landArea = gFeat['area_ha']
                    # Check if polygon is multipart 
                    if geom.isMultipart() is False: # if only simple polygon, calculate only for this
                        polyg = geom.asPolygon() #transform to list of points
                        if len(polyg) > 0:
                            area = calculator.measurePolygon(polyg[0])
                            # Convert sq km to ha
                            landArea = area * 100
                    else: 
                        multi = geom.asMultiPolygon()
                        for polyg in multi:
                            area = area + calculator.measurePolygon(polyg[0])
                            landArea = area * 100
                            print "Calculated on a multipart polygon)

