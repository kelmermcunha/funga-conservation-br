## Spatial analysis presented in the manuscript entitled "Challenges and Perspectives for Fungal Conservation in Brazil"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Corresponding author: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


import warnings
import math
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib.colors

from shapely.geometry import Point, Polygon
from shapely.geometry import Polygon, box
from matplotlib_scalebar.scalebar import ScaleBar


# Loading harmonised occurrence dataset generated in data_precleaning.py
# Choose between scenarios presented in the manuscript. Remember to change code bellow accordingly.
occ_strict = './data/output/occurrences_harmonised_strict.csv'
# occ_relaxed = './data/output/occurrences_harmonised_relaxed.csv'

if occ_strict:
    occurrences_harmonised = pd.read_csv(occ_strict)
else:
    occurrences_harmonised = pd.read_csv(occ_relaxed)

# Loading Brazilian biome and municipalities shapefiles
biomes = gpd.read_file('data/dashboard_biomes-static-layer/dashboard_biomes-static-layer.shp')
biomes = biomes.to_crs("EPSG:4326")

muni = gpd.read_file('data/BR_Municipios_2022/BR_Municipios_2022.shp')
muni = muni.to_crs("EPSG:4326")

# Georreferenced data
## General idea: extract only occurrences where 'decimalLongitude' is not empty (meaning that there is
## coordinates associated), transform it into a dataframe with geometry, and check in which province (sensu Morrone et al., 2022) each of
## the occurrence points falls within.
CoordsMatrix = occurrences_harmonised[~occurrences_harmonised['decimalLongitude'].isna()]
CoordsMatrix = CoordsMatrix[CoordsMatrix['decimalLongitude']!=0]
CoordsMatrix = gpd.GeoDataFrame(
    CoordsMatrix, geometry=gpd.points_from_xy(CoordsMatrix.decimalLongitude, CoordsMatrix.decimalLatitude)
)
CoordsMatrix = CoordsMatrix.set_crs('epsg:4326')
joinedMatrix = gpd.sjoin(CoordsMatrix, biomes, how='inner', op='within')
CoordsMatrix = joinedMatrix[['current_name', 'name', 'geometry']].fillna('')

## Removing duplicates
print('Removing duplicates from georreferenced harmonised occurrences')
derep_occs = gpd.GeoDataFrame(columns=CoordsMatrix.columns)
duplicate_occs = gpd.GeoDataFrame(columns=CoordsMatrix.columns)

for sp in CoordsMatrix['current_name'].unique():
    subset = CoordsMatrix[CoordsMatrix['current_name'] == sp]

    duplicates_mask = subset.duplicated(subset=['geometry'], keep=False)
    duplicates = subset[duplicates_mask]
    duplicate_occs = pd.concat([duplicate_occs, duplicates])

    unique_subset = subset.drop_duplicates(subset=['geometry'])
    derep_occs = pd.concat([derep_occs, unique_subset])

print(f"Total duplicate records: {len(duplicate_occs)}")
print(f"Total unique records after deduplication: {len(derep_occs)}")

# Non-georreferenced data
## General idea: extract occurrences where the 'decimalLongitude' is empty (meaning that there is no associated
## coordinates), and use the county and municipality columns to match it with a shapefile of Brazil municipalities.
## Where there is a match, a centroid point is extracted from the municipality polygon, giving the occurrence an
## approximated coordinate. With the generated coordinate, the same method used in georreferenced data was used,
## to check in which province the occurrence point falls within.
noCoordsMatrix = occurrences_harmonised[occurrences_harmonised['decimalLongitude'].isna()]

# County column
NCcountyMatrix = noCoordsMatrix[~noCoordsMatrix['county'].isna()]

warnings.filterwarnings("ignore")

NCcountyMatrix['county'] = NCcountyMatrix['county'].replace(r'Mun\..*', '', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace(r'\d+', '', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace(r',.*', '', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace(r'-.*', '', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace(r'\(.*', '', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace('S\.', 'São', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace('Sta\.', 'Santa', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace('Ten\.', 'Tenente', regex=True)
NCcountyMatrix['county'] = NCcountyMatrix['county'].replace('/.*', '', regex=True)
NCcountyMatrix = NCcountyMatrix.rename(columns={'county':'NM_MUN'})
NCcountyMatrix = NCcountyMatrix.merge(muni, on='NM_MUN', how='inner')
NCcountyMatrix = NCcountyMatrix[NCcountyMatrix['geometry'].notna()]
NCcountyMatrix = gpd.GeoDataFrame(NCcountyMatrix, geometry=NCcountyMatrix['geometry'])
NCcountyMatrix['geometry'] = NCcountyMatrix['geometry'].centroid
NCcountyMatrix = gpd.GeoDataFrame(NCcountyMatrix, geometry=NCcountyMatrix['geometry'], crs="EPSG:4326")
joined = gpd.sjoin(NCcountyMatrix, biomes, how='inner', op='within')
NCcountyMatrix = joined[['current_name', 'name', 'geometry']].fillna('')

# Municipality column
NCmuniMatrix = noCoordsMatrix[~noCoordsMatrix['municipality'].isna()]
import warnings
warnings.filterwarnings("ignore")

NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace(r'Mun\..*', '', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace(r'\d+', '', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace(r',.*', '', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace(r'-.*', '', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace(r'\(.*', '', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace('S\.', 'São', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace('Sta\.', 'Santa', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace('Ten\.', 'Tenente', regex=True)
NCmuniMatrix['municipality'] = NCmuniMatrix['municipality'].replace('/.*', '', regex=True)
NCmuniMatrix = NCmuniMatrix.rename(columns={'municipality':'NM_MUN'})
NCmuniMatrix = NCmuniMatrix.merge(muni, on='NM_MUN', how='inner')
NCmuniMatrix = NCmuniMatrix[NCmuniMatrix['geometry'].notna()]
NCmuniMatrix = gpd.GeoDataFrame(NCmuniMatrix, geometry=NCmuniMatrix['geometry'])
NCmuniMatrix['geometry'] = NCmuniMatrix['geometry'].centroid
NCmuniMatrix = gpd.GeoDataFrame(NCmuniMatrix, geometry=NCmuniMatrix['geometry'], crs="EPSG:4326")
joined = gpd.sjoin(NCmuniMatrix, biomes, how='inner', op='within')
NCmuniMatrix = joined[['current_name', 'name', 'geometry']].fillna('')

# Combining both georreferenced and non-georreferenced dataframes
combinedDataMatrix = pd.concat([CoordsMatrix, NCcountyMatrix, NCmuniMatrix], ignore_index=True, sort=False)
combinedDataMatrix['current_name'].nunique() ## Checking total species number after data treatment.

# Plotting province species number and occurrence points spatial pattern
## Loading Morrone et al., (2022) Neotropical provinces map
neo = gpd.read_file('data/neotropicalBioregionsSHP/NeotropicMap_Geo.shp')
neo = neo.to_crs('EPSG:4326')

## Joining provinces with species occurrence points
provinces = gpd.sjoin(combinedDataMatrix, neo, how='inner', op='within')
provinces = provinces[['current_name', 'name', 'Provincias', 'geometry']].fillna('')

sppNeoMatrix = provinces.pivot_table(index='Provincias', columns='current_name', aggfunc='size', fill_value=0)
sppNeoMatrix[sppNeoMatrix > 0] = 1

## Mapping species numbers to each province and plotting it (draft map that was edited for publication)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#eab89f","#d28d6a","#b76237","#993600"])

points = gpd.GeoSeries(
    [Point(-60,-20), Point(-61,-20)], crs=4326
)
points = points.to_crs(32619)
distance_meters = points[0].distance(points[1])

spp_series = sppNeoMatrix.sum(axis=1, numeric_only=True)
neo['sppNumber'] = neo['Provincias'].map(spp_series).fillna(0)
variable = 'sppNumber'
vmin, vmax = 0, math.ceil(max(neo['sppNumber']/100))*100
fig, ax = plt.subplots(1, figsize=(30,10))
ax.set_xlim([-80, -30])
ax.set_ylim([-35, 8])
ax.add_artist(ScaleBar(distance_meters, units='km', location='lower right', length_fraction=0.1))
ax.set_xticks([])
ax.set_yticks([])



sm = plt.cm.ScalarMappable(cmap='Blues', 
                           norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
fig.colorbar(sm, shrink=0.2, pad=0.0005, ticks=mticker.MultipleLocator(1500), ax=ax)
neo.plot(column=variable, cmap='Blues', linewidth=0.8, ax=ax)
combinedDataMatrix.plot(ax=ax, markersize=0.1, color='#ebcc34', alpha=0.8)
plt.savefig('./data/output/provinceMap_plasma.png', dpi=300, bbox_inches='tight')

if occ_strict:
    combinedDataMatrix.to_csv('./data/output/combinedDataMatrix_strict.csv')
else:
    combinedDataMatrix.to_csv('./data/output/combinedDataMatrix_relaxed.csv')