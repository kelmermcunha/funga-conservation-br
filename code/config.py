## Configuration script for running associated code in the manuscript entitled 
## "Brazil as a global player in Fungal Conservation: A rapid shift from neglect to Action"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Contact: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)

data = '/Users/kelmermcunha/Documents/IUCN/Fungal Conservation in Brazil/submission/2025/code/data/'
output = '/Users/kelmermcunha/Documents/IUCN/Fungal Conservation in Brazil/submission/2025/code/output_test/'

## GBIF occurrences — choose between datasets to run each scenario presented in the manuscript. 
## For running with relaxed dataset, change to 'False' the 'use_strict' object bellow.
use_strict = True
occ_strict = data+'/gbif/occurrence_strict.csv'
occ_relaxed = data+'/gbif/occurrence.txt'

## Path to MycoBank
mb_path = data+'MBList_2025_2.xlsx'

## Jaro Winkler distance threshold for passing occurrences full scientific names to second fuzzy match iteration
jaro_threshold = 0.15

## Paths to manually verified names in a .xlsx file
verified_manually_strict = output+'verified_names_manually_strict.xlsx'
verified_manually_relaxed = output+'verified_names_manually_relaxed.xlsx'

## Boolean to define if spatial analysis will be performed with strict or relaxed dataset.
## Set according to dataset used for taxonomic harmonisation and total species estimates.
## For running with relaxed dataset, change to 'False' the 'use_strict_spatial' object bellow.
use_strict_spatial = True

## Paths for both strict and relaxed occurrences dataset after harmonisation (for usage in spatial analysis)
occ_strict_harmonised = output+'occurrences_harmonised_strict.csv'
occ_relaxed_harmonised = output+'occurrences_harmonised_relaxed.csv'

## Path for shapefiles used in the spatial analysis:
biome_path = data+'dashboard_biomes-static-layer/dashboard_biomes-static-layer.shp'
muni_path = data+'BR_Municipios_2022/BR_Municipios_2022.shp'
neo_path = data+'neotropicalBioregionsSHP/NeotropicMap_Geo.shp'

## Boolean to define if endemic analysis will be performed with strict or relaxed dataset. Also affects the second
## iteration of usageKey gathering behavior. Set according to dataset used for taxonomic harmonisation, 
## total species estimates, and spatial analysis. For running with relaxed dataset, change to 'False' the 'use_strict_endemic' object bellow.
use_strict_endemic = True

## Paths for both strict and combined matrix occurrences dataset after spatial analysis (for usage in endemic analysis)
comb_matrix_strict = output+'combinedDataMatrix_strict.csv'
comb_matrix_relaxed = output+'combinedDataMatrix_relaxed.csv'

## GBIF API url to gather species usageKeys
api_url = "https://api.gbif.org/v1/species/match"

## Paths for both strict and combined dereplicated occurrences (for usage in endemic analysis usageKey gathering second iteration)
derep_occs_strict = output+'derep_occs_strict.csv'
derep_occs_relaxed = output+'derep_occs_relaxed.csv'

## Batch size to divide species names into queries (recommended to leave at 5000)
batch_size = 5000

## GBIF credentials for sending requests 
## Change username:password to your personal GBIF account and password. These files will be available for download in the GBIF website (https://www.gbif.org)
username_password = 'username:password'