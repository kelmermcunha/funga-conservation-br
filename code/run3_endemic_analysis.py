import logging
import warnings
import pandas as pd
from handlers.endemic_handlers import read_combined_matrix
from handlers.endemic_handlers import gather_gbif_keys
from handlers.endemic_handlers import gather_gbif_keys_original
from handlers.endemic_handlers import create_queries_gbif
from handlers.endemic_handlers import send_requests
from handlers.endemic_handlers import perform_endemic_analysis
from config import *


## Start log file
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("run3_endemic_analysis.log", mode='a'),
        logging.StreamHandler()
    ]
)

logging.captureWarnings(True)
warnings.simplefilter('default')

combinedDataMatrix = read_combined_matrix(use_strict_endemic=True, comb_matrix_strict=comb_matrix_strict, comb_matrix_relaxed=comb_matrix_relaxed)

keys_gbif, spp_no_keys = gather_gbif_keys(combinedDataMatrix)

keys_gbif2 = gather_gbif_keys_original(spp_no_keys, use_strict_endemic=True)

key_gbif_conct = create_queries_gbif(keys_gbif, keys_gbif2, batch_size=batch_size)

send_requests(output=output)

br_endemics = perform_endemic_analysis(key_gbif_conct)