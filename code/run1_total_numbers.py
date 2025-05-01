import logging
import warnings
from handlers.taxonomic_handlers import format_occurrences
from handlers.taxonomic_handlers import shs_treatment
from handlers.taxonomic_handlers import join_df_shs
from handlers.taxonomic_handlers import perform_harmonisation
from config import *

## Start log file
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("run1_total_numbers.log", mode='a'),
        logging.StreamHandler()
    ]
)

logging.captureWarnings(True)
warnings.simplefilter('default')

## Perform taxonomic harmonisation and run analysis for total species estimates
df_species = format_occurrences(use_strict=True, occ_strict=occ_strict, occ_relaxed=occ_relaxed)

shs = shs_treatment(df_species)

df_species = join_df_shs(df_species, shs)

perform_harmonisation(df_species, mb_path)