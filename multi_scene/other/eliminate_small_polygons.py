import sys
from qgis.core import *

import os
import urllib
import zipfile
import tempfile

# Initialize QGIS Application
QgsApplication.setPrefixPath("C:\\OSGeo4W64\\apps\\qgis", True)
app = QgsApplication([], True)
QgsApplication.initQgis()

# Add the path to Processing framework
sys.path.append('c:\\Users\\Ujaval\\.qgis2\\python\\plugins')

# Import and initialize Processing framework
from processing.core.Processing import Processing
Processing.initialize()
import processing
