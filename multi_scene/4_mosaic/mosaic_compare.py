#!/usr/bin/env python
from __future__ import division
import os
import csv

import numpy as np
from osgeo import gdal

gdal.UseExceptions()
gdal.AllRegister()

#  Read and open datasets
os.chdir('/projectnb/landsat/projects/Colombia/Mosaics')

#  Define lookup values for transitions of interest
allcl = [1, 2, 3, 4, 5, 6, 7]
lut = {}

#  Pixels classified as zero correspond to all the non-specified transitions
#  and they are the reason why the printed percentages don't add up to 100%
#  Forest
lut[1] = ([1], [1])
#  Grassland
lut[2] = ([2], [2])
#  Urban
lut[3] = ([3], [3])
#  Pastures/crops(PC)
lut[4] = ([4], [4])
#  Regrowth
lut[5] = ([5], [5])
#  Water
lut[6] = ([6], [6])
#  Other
lut[7] = ([7], [7])
#  Forest to PC
lut[8] = ([1], [4])
#  Forest to regrowth
lut[9] = ([1], [5])
#  Forest to others
lut[10] = ([1], [2, 3, 6, 7])
#  PC to regrowth
lut[11] = ([4], [5])
#  PC to forest
lut[12] = ([4], [1])
#  Regrowth to PC
lut[13] = ([5], [4])
#  Regrowth to forest
lut[14] = ([5], [1])
#  All classes to unclassified
lut[15] = (allcl, [0])

out_class = ['Other transition', 'Forest', 'Grassland', 'Urban', 'Pastures/Crops',
             'Regrowth', 'Water', 'Other', 'F -> PC', 'F -> Regr', 'F -> Other',
             'PC -> Regr', 'PC -> F', 'Regr -> PC', 'Regr -> F', 'All to UC']  

for y in range(0, 12):
	yr1 = str(y).zfill(2) 
	yr2 = str(y+1).zfill(2) 
	f1 = 'mosaic_20' + yr1 + '-masked.tif'
	f2 = 'mosaic_20' + yr2 + '-masked.tif'

	#  Open datasets
	ds1 = gdal.Open(f1, gdal.GA_ReadOnly)
	ds2 = gdal.Open(f2, gdal.GA_ReadOnly)
	if ds1 is None or ds2 is None:
    	    print 'Could not open'
    	    sys.exit(1)

	#  Load into arrays
	a1 = ds1.GetRasterBand(1).ReadAsArray()
	a2 = ds2.GetRasterBand(1).ReadAsArray()

	#  Init output
	map = np.zeros_like(a1)

	#  Housekeeping variables
	npix = (a1 != 0).sum()
	total = 0

	print 'Total number of pixels: {n}'.format(n=npix)

	for item in lut.iteritems():
		#  Unpack
		label = item[0]
		c1 = item[1][0]
		c2 = item[1][1]
		#  Matches
		m1 = np.in1d(a1, c1).reshape(a1.shape)
		m2 = np.in1d(a2, c2).reshape(a2.shape)
		m_1_2 = np.logical_and(m1, m2)

		#Housekeeping
		n = m_1_2.sum()
		
		p=round(n / npix*100,3)
		out_txt = 'stats_' + yr1 + '-' + yr2 + '_mask12.csv'
		
		with open(out_txt, 'a') as csvfile:
			stwriter = csv.writer(csvfile, delimiter=',')
			stwriter.writerow([out_class[label], n, p])
		
		total = total + n / npix * 100
		#  Add to map
		map = map + label * m_1_2

	print 'Applied LUT for {t}% of pixels'.format(t=total)

	#  Output
	#out_name = 'strata_map_' + yr1 + '-' + yr2 + '_mask12'
	#out_driver = gdal.GetDriverByName('ENVI')
	#out_ds = out_driver.Create(out_name,
	#		ds1.RasterXSize, ds1.RasterYSize, 1, gdal.GDT_Byte)
	#b = out_ds.GetRasterBand(1)
	#b.WriteArray(map)
	#b.SetCategoryNames(out_class)
	#out_ds.SetProjection(ds1.GetProjection())
	#out_ds.SetGeoTransform(ds1.GetGeoTransform())

	#out_ds = None

	#print 'Wrote map out to {f}'.format(f=out_name)

