## Analysis for endemic species inferences presented in the manuscript entitled "Challenges and Perspectives for Fungal Conservation in Brazil"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Corresponding author: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


## Packages

import geopandas as gpd
import pandas as pd
import os
import re
import sys
import time
import requests
import json
import numpy as np
import requests
import time
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from rapidfuzz.distance import JaroWinkler
from rapidfuzz.process import extract


def process_species(index, sp):
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
        print(f"\nException for {sp}: {str(e)}")
        return (index, 'NA', 'NA')
    
    return (index, temp_key, temp_name)
def print_progress(current, total, start_time):
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

def main():
    # Loading combined occurrences matrix (both georreferenced and non-georreferenced records) generated in spatial_analysis.py
    # Choose between scenarios presented in the manuscript. Remember to change code bellow accordingly.
    combmatrix_strict = './data/output/combinedDataMatrix_strict.csv'
    # combmatrix_relaxed = './data/output/combinedDataMatrix_relaxed.csv'
    
    if combmatrix_strict:
        combinedDataMatrix = pd.read_csv(combmatrix_strict)
    else:
        combinedDataMatrix = pd.read_csv(combmatrix_relaxed)
    
    sppNames = combinedDataMatrix.drop_duplicates(subset=['current_name'])
    sppNames = list(sppNames['current_name'].unique())
    
    # Connecting with GBIF api in order to gather occurrences for all names (= len(sppNames)) and check which species have occurrences only for Brazil.
    # The 'process_species' function connects with the API and takes the GBIF key number for each species, which will be used to download occurrences.
    api_url = "https://api.gbif.org/v1/species/match"
    
    ## Initialize result lists
    keys = ['NA'] * len(sppNames)
    keys_names = ['NA'] * len(sppNames)
    spp_no_keys = []
    
    batch_size = 100
    delay_between_batches = 0.7
    max_workers = 10
    total_species = len(sppNames)
    
    ## Starting time and running progress and process species function
    start_time = time.time()
    processed_count = 0
    
    for i in range(0, total_species, batch_size):
        batch = sppNames[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_species, i+j, sp) 
                      for j, sp in enumerate(batch)]
            
            for future in as_completed(futures):
                index, temp_key, temp_name = future.result()
                keys[index] = temp_key
                keys_names[index] = temp_name
                if temp_key == 'NA':
                    spp_no_keys.append(sppNames[index])
                
                processed_count += 1
                print_progress(processed_count, total_species, start_time)
        
        if i + batch_size < total_species:
            time.sleep(delay_between_batches)
    
    ## Clear the progress bar when done
    sys.stdout.write('\r' + ' ' * 120 + '\r')  # Clear the line
    sys.stdout.flush()
    
    ## Post-processing
    spp_no_keys = [sppNames[i] for i in range(len(sppNames)) if keys[i] == 'NA']
    
    print(f"Completed! Processed {total_species} species in {time.time()-start_time:.1f} seconds")
    print(f"Species without keys found: {len(spp_no_keys)}") ## These will be left out for a second iteration, 
                                                             ## as GBIF did not find any match with the valid name according to MycoBank
    
    keys_gbif = pd.DataFrame({'scientificName':sppNames, 'matchName':keys_names, 'gbifKey':keys})
    keys_gbif = keys_gbif[keys_gbif['gbifKey']!='NA']
    
    # For species with no keys, run a second fetching iteration by using the original taxon name associated with occurrences (original column 'scientificName')
    # This information will be gathered from 'derep_occs.csv'
    if combmatrix_strict:
        derep_occs = pd.read_csv('./data/output/derep_occs_strict.csv')
    else:
        derep_occs = pd.read_csv('./data/output/derep_occs_relaxed.csv')
    
    keys_na_names = derep_occs[derep_occs['current_name'].isin(spp_no_keys)]
    
    sppNames2 = list(keys_na_names['scientificName'].unique()) ## .unique() not strictly necessary
    ## Initialize result lists
    keys2 = ['NA'] * len(sppNames2)
    keys_names2 = ['NA'] * len(sppNames2)
    spp_no_keys2 = []
    
    ## Starting time and running progress and process species function
    start_time = time.time()
    processed_count = 0
    
    for i in range(0, total_species, batch_size):
        batch = sppNames2[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_species, i+j, sp) 
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
    
    spp_no_keys2 = [sppNames2[i] for i in range(len(sppNames2)) if keys[i] == 'NA']
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
    
    # Merge
    keys_gbif_conct = pd.concat([keys_gbif, keys_gbif2], ignore_index=True)
    
    # Creating queries to pass to GBIF API in order to download occurrence data for all species.
    taxon_keys = keys_gbif_conct['gbifKey'].to_list()
    num_keys = len(taxon_keys)
    batch_size = 5000  # Number of taxon keys to include in each query
    
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
        query_file_name = f"./code/data/output/query_{i+1}_{i+len(batch_keys)}.json"
        with open(query_file_name, "w") as f:
            json.dump(query_dict, f)
    
    json_files = [file for file in os.listdir('./code/data/output') 
                  if '.json' in file and os.path.isfile(os.path.join('./code/data/output', file))]
    
    
    ## Sending requests based on .json files
    ## Code used to request the download
    for i in json_files:
        os.system(f'curl --include --user username:password --header "Content-Type: application/json" --data @{i} https://api.gbif.org/v1/occurrence/download/request')## Change username:password to your account username and password
        time.sleep(10)
    ## Change username:password to your personal GBIF account and password. These files will be available for download in the GBIF website (https://www.gbif.org)
    ## Save it in /code/data/file.zip to continue to run the analysis
    ## IMPORTANT: if more than three downloads, it will fail. Need to tweak.
    
    print('All download requests were passed into GBIF. Download and save then in /code/data as .zip files.')
    while True:
        user_input = input('Have you already saved the files? (y/n): '.strip().lower())
        if user_input == 'y':
            print('Continuing...')
            break
        elif user_input == 'n':
            print('Save files to continue.')
            time.sleep(2)
        else:
            print('Invalid input. Enter "y" or "n"')
    print('Resuming...')
    
    # Reading .zip files and inferring endemic species total number
    keys_gbif_conct['binomial'] = keys_gbif_conct['matchName'].str.split(' ').str[:2].str.join(' ')
    binomials = keys_gbif_conct['binomial']
    
    zip_files = [file for file in os.listdir('./code/data') 
                  if '.zip' in file and os.path.isfile(os.path.join('./code/data', file))]
    br_endemics = []
    country_sets = defaultdict(set)
    
    for zipf in zip_files:
        for chunk in pd.read_csv(
            './code/data/'+zipf,
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
    
    ## Find all single-country species
    single_country_species = {
        sp: cc for sp, cc in country_sets.items() 
        if len(cc) == 1
    }
    
    ## Filter for match with 'binomials' and country set equal to 'BR'
    br_endemics = [
        sp for sp in binomials 
        if sp in single_country_species 
        and next(iter(single_country_species[sp])) == 'BR'
    ]

if __name__ == "__main__":
    main()
