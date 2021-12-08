# PWF
Quick-and-dirty matlab script for replicating Fig. 12 in the paper by Piana Agostinetti et al. (2021) https://doi.org/10.5194/se-2021-37

The code serves for analyzing the DAS recording of the Mw 4.3 Hawtorne earthquake discussed in Piana Agostinetti et al. 
(2021); the final output is Figure 12 in that paper. 

Earthquake information :
M 4.3 - 23km ESE of Hawthorne, Nevada
Time: 2016-03-21 07:37:10 (UTC)
Location: 38.479 N 118.366 W
Depth: 9.9 km

Data are also provided: they are sac files cointaining stacked
DAS channels from the PoroTomo esperiment (University of Wisconsin, 2016). 
Brady's Geothermal Field DAS Earthquake Data [data set].  Retrieved from https://dx.doi.org/10.15121/1334285.
Station coordinates are in Stla, Stlo header fields.

Additional functions (posted separately):
Readsac -> to read sac data files
wgs2utm -> to convert (lat,lon) to UTM coordinates

G. Saccorotti, INGV. 8 December 2021
