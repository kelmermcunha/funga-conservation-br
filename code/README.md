This file describes the scripts found in the `code` directory. They contain all code used in the pipeline from dataset preparation, taxonomic harmonisation, data cleaning, and spatial and endemicity analyses.

### 1) Data precleaning
Pre-clean data by removing duplicates and unwanted occurrences based on the type of record. Prepares the [MycoBank](https://mycobank.org) database for taxonomic harmonisation using the auxiliary function in `format_mb.py`. If Species Hypotheses are present, pull data on taxonomy via the [PlutoF API](https://plutof.docs.apiary.io/#). Perform taxonomic harmonisation of all occurrences gathered from GBIF based on exact matches and a two-step fuzzy matching.

### 2) Spatial analysis
Estimates the total number of valid and digitally accessible fungal species names for each Neotropical province sensu [Morrone et al. (2023)]( https://doi.org/10.1590/0001-3765202220211167) that intersects with Brazilian territory based on the harmonised and pre-cleaned GBIF occurrence dataset on: a) coordinates for georeferenced occurrences and b) by generating centroids for non-georeferenced occurrences with detailed metadata for the collection locality.

### 3) Endemicity analysis
Retrieve usage Keys for creating occurrence download requests via the GBIF API. Download available occurrences for all known, digitally accessible, and accepted fungal species names retrieved in `data_precleaning.py`. Based on available occurrences for each species name, infer its endemicity status based on the associated country for each occurrence under species names.
