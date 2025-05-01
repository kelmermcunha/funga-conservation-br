import os
import sys
import time
import math
import json
import logging
import pandas as pd
import requests
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from config import *
from modules.endemic_analysis import find_endemics

def read_combined_matrix(use_strict_endemic:bool, comb_matrix_strict:str, comb_matrix_relaxed:str):
    """
    Reads the combined data matrix depending on dataset option.
    """

    if use_strict_endemic:
        combinedDataMatrix = pd.read_csv(comb_matrix_strict)
        logging.info(f"Reading combined matrix from: {comb_matrix_strict}")
    else:
        combinedDataMatrix = pd.read_csv(comb_matrix_relaxed)
        logging.info(f"Reading combined matrix from: {comb_matrix_relaxed}")

    return combinedDataMatrix

def process_species(index, sp, api_url):
    """
    Query GBIF API for taxon keys and names for a given species.
    Returns tuple of (index, GBIF key, matched scientific name).
    """

    params = {'name': sp, 'verbose': 'true'}
    temp_key = 'NA'
    temp_name = 'NA'
    
    try:
        response = requests.get(api_url, params=params)
        
        if response.status_code == 200:
            data = response.json()
            
            if ('rank' in data
                and data.get('rank') not in {'GENUS', 'KINGDOM', 'FAMILY', 'PHYLUM'}
                and data.get('confidence', 0) >= 80
                and 'usageKey' in data):
                temp_key = data['usageKey']
                temp_name = data['scientificName']

            if temp_key == 'NA' and 'alternatives' in data:
                for alt in data['alternatives']:
                    if (alt.get('rank') not in {'GENUS', 'KINGDOM', 'FAMILY', 'PHYLUM'}
                        and alt.get('confidence', 0) >= 75
                        and 'usageKey' in alt):
                        temp_key = alt['usageKey']
                        temp_name = data['scientificName']
                        break
        else:
            print(f"\nError for {sp}: Status code {response.status_code}")
            return (index, 'NA', 'NA')
    
    except Exception as e:
        logging.warning(f"Exception for {sp}: {str(e)}")
        return (index, 'NA', 'NA')
    
    return (index, temp_key, temp_name)

def print_progress(current, total, start_time):
    """
    Displays a progress bar in the console.
    """

    elapsed = time.time() - start_time
    progress = current / total
    bar_length = 40
    filled = int(bar_length * progress)
    bar = '█' * filled + '-' * (bar_length - filled)
    
    # Calculate estimated time remaining
    if current > 0:
        remaining = (elapsed / current) * (total - current)
        remaining_str = f"{remaining:.1f}s"
    else:
        remaining_str = "calculating..."
    
    # Calculate processing rate
    rate = current / elapsed if elapsed > 0 else 0
    
    sys.stdout.write(
        f'\rProcessing: |{bar}| {current}/{total} ({progress:.1%}) '
        f'| Rate: {rate:.1f} sp/s | Elapsed: {elapsed:.1f}s | Remaining: {remaining_str}'
    )
    sys.stdout.flush()

def gather_gbif_keys(combinedDataMatrix):
    """
    Retrieve GBIF taxon keys for all unique species in the dataset.
    Returns DataFrame of matched keys and list of unmatched species.
    """

    logging.info("Starting GBIF key resolution")

    sppNames = combinedDataMatrix.drop_duplicates(subset=['current_name'])
    sppNames = list(sppNames['current_name'].unique())
    sppNames = sppNames[0:100]

    ## Initialize result lists
    keys = ['NA'] * len(sppNames)
    keys_names = ['NA'] * len(sppNames)

    batch_size = 100
    delay_between_batches = 0.7
    max_workers = 10
    total_species = len(sppNames)

    logging.info(f"Total species to process: {total_species}")

    ## Starting time and running progress and process species function
    start_time = time.time()
    processed_count = 0

    for i in range(0, total_species, batch_size):
        batch = sppNames[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_species, i+j, sp, api_url) 
                    for j, sp in enumerate(batch)]
            
            for future in as_completed(futures):
                index, temp_key, temp_name = future.result()
                keys[index] = temp_key
                keys_names[index] = temp_name
                processed_count += 1
                print_progress(processed_count, total_species, start_time)
        
        if i + batch_size < total_species:
            time.sleep(delay_between_batches)

    ## Clear the progress bar when done
    sys.stdout.write('\r' + ' ' * 120 + '\r')  # Clear the line
    sys.stdout.flush()

    spp_no_keys = [sppNames[i] for i in range(len(sppNames)) if keys[i] == 'NA']

    logging.info(f"Matched keys: {total_species - len(spp_no_keys)} | Unmatched: {len(spp_no_keys)}")

    keys_gbif = pd.DataFrame({'scientificName':sppNames, 'matchName':keys_names, 'gbifKey':keys})
    keys_gbif = keys_gbif[keys_gbif['gbifKey']!='NA']

    return keys_gbif, spp_no_keys

def gather_gbif_keys_original(spp_no_keys, use_strict_endemic:bool):
    """
    Attempt GBIF key resolution for species unmatched in initial pass.
    Uses original scientific names from dereplicated occurrences.
    """

    logging.info("Processing unmatched species with original names")

    if use_strict_endemic:
        derep_occs = pd.read_csv(derep_occs_strict)
    else:
        derep_occs = pd.read_csv(derep_occs_relaxed)

    keys_na_names = derep_occs[derep_occs['current_name'].isin(spp_no_keys)]

    sppNames2 = list(keys_na_names['scientificName'].unique()) ## .unique() not strictly necessary

    ## Initialize result lists
    keys2 = ['NA'] * len(sppNames2)
    keys_names2 = ['NA'] * len(sppNames2)
    spp_no_keys2 = []

    batch_size = 100
    delay_between_batches = 0.7
    max_workers = 10
    total_species = len(sppNames2)

    logging.info(f"Unmatched species to reprocess: {total_species}")

    ## Starting time and running progress and process species function
    start_time = time.time()
    processed_count = 0

    for i in range(0, total_species, batch_size):
        batch = sppNames2[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_species, i+j, sp, api_url) 
                    for j, sp in enumerate(batch)]
            
            for future in as_completed(futures):
                index, temp_key, temp_name = future.result()
                keys2[index] = temp_key
                keys_names2[index] = temp_name
                if temp_key == 'NA':
                    spp_no_keys2.append(sppNames2[index])
                
                processed_count += 1
                print_progress(processed_count, total_species, start_time)
        
        if i + batch_size < total_species:
            time.sleep(delay_between_batches)

    ## Clear the progress bar when done
    sys.stdout.write('\r' + ' ' * 120 + '\r')  # Clear the line
    sys.stdout.flush()

    spp_no_keys2 = [sppNames2[i] for i in range(len(sppNames2)) if keys2[i] == 'NA']

    logging.info(f"Reprocessing complete | New matches: {total_species - len(spp_no_keys2)} | Remaining unmatched: {len(spp_no_keys2)}")

    print(f"Completed! Processed {total_species} species in {time.time()-start_time:.1f} seconds")
    print(f"Species without keys found: {len(spp_no_keys2)}")

    keys_gbif2 = pd.DataFrame({'scientificName':sppNames2, 'matchName':keys_names2, 'gbifKey':keys2})

    # Create a mapping dictionary from derep to retrieve harmonised names based on MycoBank
    name_mapping = dict(zip(derep_occs['scientificName'], derep_occs['current_name']))
    # Update
    keys_gbif2['current_name'] = keys_gbif2['scientificName'].map(name_mapping)

    del keys_gbif2['scientificName']
    keys_gbif2.rename(columns={'current_name': 'scientificName'}, inplace=True)
    keys_gbif2 = keys_gbif2[['scientificName', 'matchName', 'gbifKey']]

    return keys_gbif2

def create_queries_gbif(keys_gbif, keys_gbif2, batch_size=batch_size):
    """
    Generate GBIF download query JSON files for matched taxon keys.
    Returns concatenated DataFrame of all resolved keys.
    """

    logging.info("Creating GBIF download queries")

    # Merge
    keys_gbif_conct = pd.concat([keys_gbif, keys_gbif2], ignore_index=True)

    # Creating queries to pass to GBIF API in order to download occurrence data for all species.
    taxon_keys = keys_gbif_conct['gbifKey'].to_list()
    num_keys = len(taxon_keys)
    batch_size = batch_size  # Number of taxon keys to include in each query

    for i in range(0, num_keys, batch_size):
        batch_keys = taxon_keys[i:i+batch_size]
        query_dict = {
            "creator": "YOUR_USERNAME", ## change username
            "notification_address": ["YOUR_EMAIL@EMAIL.COM"], ## change email
            "sendNotification": True,
            "format": "SIMPLE_CSV",
            "predicate": {
                "type": "and",
                "predicates": [
                    {
                        "type": "in",
                        "key": "TAXON_KEY",
                        "values": batch_keys
                    }
                ]
            }
        }
        query_file_name = f"{output}query_{i+1}_{i+len(batch_keys)}.json"
        with open(query_file_name, "w") as f:
            json.dump(query_dict, f)

    logging.info(f"Generated {math.ceil(num_keys / batch_size)} query files")
    return keys_gbif_conct

def send_requests(output):
    """
    Submit GBIF download requests for all query JSON files.
    """

    logging.info("Submitting GBIF download requests")

    json_files = [file for file in os.listdir(output) 
            if '.json' in file and os.path.isfile(os.path.join(output, file))]
    
    logging.info(f"Found {len(json_files)} query files to submit")
    
    for i in json_files:
        os.system(f'curl --include --user {username_password} --header "Content-Type: application/json" --data @{i} https://api.gbif.org/v1/occurrence/download/request')
        time.sleep(10)

    logging.info("All download requests submitted to GBIF")

def perform_endemic_analysis(keys_gbif_conct):
    """
    Process downloaded GBIF occurrence data to identify endemic species.
    """

    logging.info("Starting endemic species analysis")

    n_batches = math.ceil(len(keys_gbif_conct)/batch_size)

    zip_files = [file for file in os.listdir(data) 
        if '.zip' in file and os.path.isfile(os.path.join(data, file))]

    while len(zip_files) != n_batches:
        user_input = input("Please save download requests into data repository. Have you saved? Answer 'no' to abort. (y/n): ".strip().lower())
        if user_input == 'y':
            zip_files = [file for file in os.listdir(data) 
                if '.zip' in file and os.path.isfile(os.path.join(data, file))]
            if len(zip_files) == n_batches:
                break
        elif user_input == 'n':
            raise FileNotFoundError(f'Length of zipped files in data ({len(zip_files)}) is different from number of batches ({n_batches}), aborting.')
        else:
            print('Invalid input. Enter "y" or "n"')

    # Reading .zip files and inferring endemic species total number
    keys_gbif_conct['binomial'] = keys_gbif_conct['matchName'].str.split(' ').str[:2].str.join(' ')
    binomials = keys_gbif_conct['binomial']

    zip_files = [file for file in os.listdir(data) 
                if '.zip' in file and os.path.isfile(os.path.join(data, file))]
    
    logging.info(f"Found {len(zip_files)} GBIF download files")

    country_sets = defaultdict(set)

    for zipf in zip_files:
        logging.debug(f"Processing {zipf}")
        for chunk in pd.read_csv(
            data+zipf,
            sep="\t",
            chunksize=500_000,
            dtype={'species': 'category', 'countryCode': 'category'},
            usecols=['species', 'countryCode'],
            on_bad_lines='skip',
            low_memory=False,
            engine='c'
        ):
            # Filling NAs with 'UNK'
            chunk['countryCode'] = chunk['countryCode'].cat.add_categories(['UNK']).fillna('UNK')
            
            # Update country sets per species
            for sp, cc_set in chunk.groupby('species', observed=False)['countryCode']:
                country_sets[sp].update(cc_set.unique())

    br_endemics = find_endemics(country_sets, binomials)
    br_endemics.to_csv(output+'br_endemics.csv')
    
    logging.info(f"Identified {len(br_endemics)} endemic species")

    return br_endemics