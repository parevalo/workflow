Files in this folder include:

- original_lut.csv: Lookup table (LUT) containing the original stratification
used for creating the 2001-2016 map.
- defor_lut.csv: LUT with the original stratification but the change from forest
to secondary forest is also considered deforestation, and labeled as 8 (i.e
forest to pasture).
- for_nofor_lutA.csv: New LUT to map stable forest (1) , stable no forest (2), 
and forest loss (3) only. In this version, the secondary forest class is 
considered non forest. Also, no forest gain is not considered as a valid 
transition, so any pixel that shows this behavior is classified as no-forest.
- For_nofor_lutB.csv: Same as the previous but secondary forest is considered forest
- poststrat_buffer_lut.csv: LUT to preserve all values from original stratification
(final_strata_01_16_UTM18N.tif) except those of stable forest that intersect a 
buffer (e.g. final_strata_01_UTM18_buffer20.tif) that goes into the forest, for which 
a new value of 16 is assigned. Only possible combinations are included in the LUT.
- lut_strata_ref_match: Lookup table used to assign "match" or "mismatch" labels 
to the sample units to create a visualization and perform an exploratory analysis. 
