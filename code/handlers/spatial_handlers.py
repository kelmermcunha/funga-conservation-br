import math
import logging
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from shapely.geometry import Point
from matplotlib_scalebar.scalebar import ScaleBar

from config import *
from modules.spatial_analysis import *


def read_harmonised_occurrences(use_strict_spatial:bool, occ_strict_harmonised:str, occ_relaxed_harmonised:str):
    """
    Load harmonised occurrence records based on strict or relaxed dataset setting.
    """

    logging.info("Reading harmonised occurrences")

    if use_strict_spatial:
        occurrences_harmonised = pd.read_csv(occ_strict_harmonised)
    else:
        occurrences_harmonised = pd.read_csv(occ_relaxed_harmonised)

    return occurrences_harmonised

def treat_georeferenced(occurrences_harmonised):
    """
    Convert georeferenced occurrence records into a GeoDataFrame.
    Filters out records with missing or zero coordinates.
    """

    logging.info("Processing georeferenced occurrences")

    CoordsMatrix = occurrences_harmonised[~occurrences_harmonised['decimalLongitude'].isna()]
    CoordsMatrix = CoordsMatrix[CoordsMatrix['decimalLongitude']!=0]
    CoordsMatrix = gpd.GeoDataFrame(
        CoordsMatrix, geometry=gpd.points_from_xy(CoordsMatrix.decimalLongitude, CoordsMatrix.decimalLatitude)
    )

    return CoordsMatrix

def perform_georeferenced_analysis(CoordsMatrix, biome_path):
    """
    Conduct biome-level analysis and deduplication on georeferenced occurrences.
    """

    logging.info("Starting georeferenced analysis")

    CoordsMatrix = georeferenced_analysis(CoordsMatrix, biome_path)

    derep_occs = gpd.GeoDataFrame(columns=CoordsMatrix.columns)
    duplicate_occs = gpd.GeoDataFrame(columns=CoordsMatrix.columns)

    for sp in CoordsMatrix['current_name'].unique():
        subset = CoordsMatrix[CoordsMatrix['current_name'] == sp]

        duplicates_mask = subset.duplicated(subset=['geometry'], keep=False)
        duplicates = subset[duplicates_mask]
        duplicate_occs = pd.concat([duplicate_occs, duplicates])

        unique_subset = subset.drop_duplicates(subset=['geometry'])
        derep_occs = pd.concat([derep_occs, unique_subset])

    logging.info(f"Unique georeferenced records: {len(derep_occs)}")

    return CoordsMatrix

def treat_nongeoreferenced_county(occurrences_harmonised, muni_path):
    """
    Process occurrences with only 'county' spatial information.
    Attempts centroid assignment using matched shapefile.
    """

    logging.info("Processing non-georeferenced data based on county")

    muni = gpd.read_file(muni_path)
    muni = muni.to_crs('epsg:4326')

    noCoordsMatrix = occurrences_harmonised[occurrences_harmonised['decimalLongitude'].isna()]

    # County column
    NCcountyMatrix = noCoordsMatrix[~noCoordsMatrix['county'].isna()]


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

    logging.info(f"County-level records processed: {len(NCcountyMatrix)}")

    return NCcountyMatrix

def treat_nongeoreferenced_muni(occurrences_harmonised, muni_path):
    """
    Process occurrences with only 'municipality' spatial information.
    Attempts centroid assignment using matched shapefile.
    """

    logging.info("Processing non-georeferenced data based on municipality")

    muni = gpd.read_file(muni_path)
    muni = muni.to_crs('epsg:4326')

    noCoordsMatrix = occurrences_harmonised[occurrences_harmonised['decimalLongitude'].isna()]

    NCmuniMatrix = noCoordsMatrix[~noCoordsMatrix['municipality'].isna()]
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

    logging.info(f"Municipality-level records processed: {len(NCmuniMatrix)}")

    return NCmuniMatrix

def perform_nongeoreferenced_analysis(NCcountyMatrix, NCmuniMatrix):
    """
    Perform biome-level analysis on non-georeferenced records derived from county or municipality data.
    """

    logging.info("Performing non-georeferenced spatial analysis")

    NCcountyMatrix, NCmuniMatrix = nongeoreferenced_analysis(NCcountyMatrix, NCmuniMatrix, biome_path=biome_path)

    return NCcountyMatrix, NCmuniMatrix

def join_gdfs(CoordsMatrix, NCcountyMatrix, NCmuniMatrix):
    """
    Combine georeferenced and non-georeferenced GeoDataFrames into a single dataset.
    """

    logging.info("Joining all geospatially-resolved occurrence datasets")

    combinedDataMatrix = pd.concat([CoordsMatrix, NCcountyMatrix, NCmuniMatrix], ignore_index=True, sort=False)

    logging.info(f"Estimate of the total number of known and accessible species with georeferenced data: {combinedDataMatrix['current_name'].nunique()} species")

    return combinedDataMatrix

def plot_results(combinedDataMatrix, neo_path):
    """
    Generate draft map showing species richness across biogeographical provinces.
    Save plot and write combined dataset to file.
    """
    
    logging.info("Plotting results and saving dataset")

    neo = gpd.read_file(neo_path)
    neo = neo.to_crs('epsg:4326')

    if 'level_0' in combinedDataMatrix.columns:
        combinedDataMatrix = combinedDataMatrix.drop(columns=['level_0'])

    provinces = gpd.sjoin(combinedDataMatrix, neo, how='inner', predicate='within')
    provinces = provinces[['current_name', 'name', 'Provincias', 'geometry']].fillna('')

    sppNeoMatrix = provinces.pivot_table(index='Provincias', columns='current_name', aggfunc='size', fill_value=0)
    sppNeoMatrix[sppNeoMatrix > 0] = 1

    ## Mapping species numbers to each province and plotting it (draft map that was edited for publication)
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
    plt.savefig(output+'provinceMap_plasma.png', dpi=300, bbox_inches='tight')

    if occ_strict:
        combinedDataMatrix.to_csv(output+'combinedDataMatrix_strict.csv')
    else:
        combinedDataMatrix.to_csv(output+'combinedDataMatrix_relaxed.csv')

    logging.info("Map and combined dataset saved")