import logging
import warnings
from handlers.spatial_handlers import read_harmonised_occurrences
from handlers.spatial_handlers import treat_georeferenced
from handlers.spatial_handlers import perform_georeferenced_analysis
from handlers.spatial_handlers import treat_nongeoreferenced_county
from handlers.spatial_handlers import treat_nongeoreferenced_muni
from handlers.spatial_handlers import perform_nongeoreferenced_analysis
from handlers.spatial_handlers import join_gdfs
from handlers.spatial_handlers import plot_results
from config import *

## Start log file
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("run2_spatial_analysis.log", mode='a'),
        logging.StreamHandler()
    ]
)

logging.captureWarnings(True)
warnings.simplefilter('default')

## Perform spatial analysis (including plots)
occurrences_harmonised = read_harmonised_occurrences(use_strict_spatial=True, occ_strict_harmonised=occ_strict_harmonised, occ_relaxed_harmonised=occ_relaxed_harmonised)

CoordsMatrix = treat_georeferenced(occurrences_harmonised)

CoordsMatrix = perform_georeferenced_analysis(CoordsMatrix=CoordsMatrix, biome_path=biome_path)

NCcountyMatrix = treat_nongeoreferenced_county(occurrences_harmonised, muni_path=muni_path)
NCmuniMatrix = treat_nongeoreferenced_muni(occurrences_harmonised, muni_path=muni_path)

NCcountyMatrix, NCmuniMatrix = perform_nongeoreferenced_analysis(NCcountyMatrix, NCmuniMatrix)

combinedDataMatrix = join_gdfs(CoordsMatrix=CoordsMatrix, NCcountyMatrix=NCcountyMatrix, NCmuniMatrix=NCmuniMatrix)

plot_results(combinedDataMatrix, neo_path=neo_path)