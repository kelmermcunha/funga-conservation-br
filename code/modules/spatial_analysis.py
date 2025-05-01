## Functions to perform spatial analysis presented in the results of the manuscript entitled 
## "Brazil as a global player in Fungal Conservation: A rapid shift from neglect to Action"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Contact: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


import geopandas as gpd
from shapely.geometry import Point, Polygon
from matplotlib_scalebar.scalebar import ScaleBar


def georeferenced_analysis(CoordsMatrix, biome_path):
    """
    Perform spatial join between georeferenced occurrence records and biomes.
    
    Parameters:
        CoordsMatrix (GeoDataFrame): Georeferenced species occurrence records.
        biome_path (str): Path to the biome shapefile.
    
    Returns:
        GeoDataFrame: Filtered dataset with species names and biome names.
    """

    biomes = gpd.read_file(biome_path)
    biomes = biomes.to_crs('epsg:4326')

    CoordsMatrix = CoordsMatrix.set_crs('epsg:4326')

    if 'level_0' in CoordsMatrix.columns:
        CoordsMatrix = CoordsMatrix.drop(columns=['level_0'])

    joinedMatrix = gpd.sjoin(CoordsMatrix, biomes, how='inner', predicate='within')
    CoordsMatrix = joinedMatrix[['current_name', 'name', 'geometry']].fillna('')

    return CoordsMatrix

def nongeoreferenced_analysis(NCcountyMatrix, NCmuniMatrix, biome_path):
    """
    Perform spatial join for non-georeferenced records (from county and municipality),
    linking them with corresponding biomes.
    
    Parameters:
        NCcountyMatrix (GeoDataFrame): Occurrence records with only county-level spatial info.
        NCmuniMatrix (GeoDataFrame): Occurrence records with only municipality-level spatial info.
        biome_path (str): Path to the biome shapefile.
    
    Returns:
        Tuple[GeoDataFrame, GeoDataFrame]: Filtered datasets (county and municipality) with species and biome info.
    """

    biomes = gpd.read_file(biome_path)
    biomes = biomes.to_crs('epsg:4326')

    if 'level_0' in NCcountyMatrix.columns:
        NCcountyMatrix = NCcountyMatrix.drop(columns=['level_0'])

    joined_county = gpd.sjoin(NCcountyMatrix, biomes, how='inner', predicate='within')
    NCcountyMatrix = joined_county[['current_name', 'name', 'geometry']].fillna('')

    if 'level_0' in NCmuniMatrix.columns:
        NCmuniMatrix = NCmuniMatrix.drop(columns=['level_0'])

    joined_muni = gpd.sjoin(NCmuniMatrix, biomes, how='inner', predicate='within')
    NCmuniMatrix = joined_muni[['current_name', 'name', 'geometry']].fillna('')

    return NCcountyMatrix, NCmuniMatrix