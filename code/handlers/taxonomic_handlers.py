import os
import re
import numpy as np
import aiohttp
from tqdm import tqdm
import pandas as pd
import asyncio
import nest_asyncio
import logging

from config import *
from modules.taxonomic_harmonisation import *
from modules.format_mb import *


def format_occurrences(use_strict:bool, occ_strict:str, occ_relaxed:str):
    """
    Orchestrate relevant column selection and formatting in GBIF downloaded data.
    """

    logging.info("Starting format_occurrences")

    if use_strict:
        rawData = pd.read_csv(occ_strict, on_bad_lines='warn', low_memory=False)
    else:
        rawData = pd.read_csv(occ_relaxed, on_bad_lines='warn', low_memory=False)
                            ## Data was downloaded on 25/03/2023 from the GBIF database 
                            ## website. The following filters were used: 'Country or area' = Brazil,
                            ## 'Scientific name' = Fungi. Reference:
                            ## GBIF.org (21 March 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.pnxaj2

    ## Remove occurrences not identified at the species level, select relevant columns, and drop duplicates based on the selected columns
    logging.info(f"'Raw' dataset size: {len(rawData)}")

    df_species = rawData.dropna(subset=['species']).copy()

    logging.info(f"Records after dropping rows with missing species: {len(df_species)}")

    if occ_strict:
        df_species = df_species[['gbifID', 'institutionCode', 'collectionCode', 'catalogNumber', 'year', 'month', 'day', 
                                'continent', 'stateProvince', 'county', 'municipality', 'locality', 'decimalLatitude',
                                'decimalLongitude','scientificName','species', 'acceptedScientificName', 'recordedBy']]
        
        df_species.drop_duplicates(
        subset=['species','recordedBy','institutionCode','catalogNumber','decimalLatitude','decimalLongitude',
                'year', 'month', 'day'], inplace=True)

    else:
        df_species = df_species[['gbifID', 'publisher', 'type', 'institutionCode', 'collectionCode', 'basisOfRecord', 'occurrenceID', 'catalogNumber', 
                            'eventDate', 'year', 'month', 'day', 'higherGeography', 'continent', 'countryCode', 'stateProvince', 'county', 'municipality', 'locality', 
                            'verbatimLocality', 'verbatimElevation', 'decimalLatitude', 'decimalLongitude', 'coordinateUncertaintyInMeters', 
                            'coordinatePrecision', 'pointRadiusSpatialFit', 'verbatimCoordinateSystem', 'georeferencedDate', 'scientificName', 'issue', 
                            'hasCoordinate', 'hasGeospatialIssues', 'species', 'acceptedScientificName', 'recordedBy']]

        df_species.drop_duplicates(
        subset=['species','recordedBy','institutionCode','catalogNumber','decimalLatitude','decimalLongitude',
                'verbatimLocality','year', 'month', 'day'], inplace=True)

    df_species['scientificName'] = [re.sub(r'\d+', '', i) for i in df_species['scientificName']]

    logging.info("Finished formatting occurrences")
    return df_species

async def fetch_sh_data(session, sh_code, retries=3):
    """
    Asynchronous function for connecting with PlutoF API service to retrieve taxonomic information for Species Hypothesis.
    These are present only in the relaxed dataset. However, the function is called for both to maintain code integrity.
    """

    api_sh = "https://api.plutof.ut.ee/v1/public/dshclusters/search/?name="
    api_taxon = "https://api.plutof.ut.ee/v1/public/taxa/"

    for attempt in range(retries):
        try:
            timeout = aiohttp.ClientTimeout(total=600)
            async with session.get(api_sh + sh_code, timeout=timeout) as sh_response:
                if sh_response.status == 200:
                    sh_data = await sh_response.json()
                    if 'attributes' in sh_data['data']:
                        binomial_author = sh_data['data']['attributes']['name'] + ' ' + sh_data['data']['attributes']['epithet_author']
                    else:
                        taxon = sh_data['data'][0]['relationships']['taxon_node']['data']['id']
                        async with session.get(api_taxon + taxon + '/') as taxon_response:
                            if taxon_response.status == 200:
                                sh_taxon = await taxon_response.json()
                                binomial_author = sh_taxon['data']['attributes']['name'] + ' ' + sh_taxon['data']['attributes']['epithet_author']
                            else:
                                binomial_author = 'NA'
                    return sh_code, binomial_author
                else:
                    return sh_code, 'NA'
        except asyncio.TimeoutError:
            return sh_code, 'NA'
        except Exception as e:
            return sh_code, 'NA'
        await asyncio.sleep(1)

async def main(shs_list, shs_code, shs_binomial_author, wait_time=1.5):
    """
    Asynchronous function for connecting with PlutoF API service to retrieve taxonomic information for Species Hypothesis (SH).
    These are present only in the relaxed dataset. However, the function is called for both to maintain code integrity.
    """
    
    if isinstance(shs_list, np.ndarray):
        shs_list = shs_list.tolist()

    timeout = aiohttp.ClientTimeout(total=600)  
    async with aiohttp.ClientSession(timeout=timeout) as session:
        tasks = [fetch_sh_data(session, sh_code) for sh_code in shs_list]
        results = []
        for task in tqdm(asyncio.as_completed(tasks), total=len(tasks)):
            result = await task
            results.append(result) 
            await asyncio.sleep(wait_time)

        sorted_results = sorted(results, key=lambda x: shs_list.index(x[0]))

        for sh_code, binomial_author in sorted_results:
            shs_code.append(sh_code)
            shs_binomial_author.append(binomial_author)

def shs_treatment(df_species):
    """
    Orchestrate SHs asynchronous functions output handling along with non-SHs information present in GBIF data.
    """

    logging.info("Starting SHs treatment")

    shs = df_species[df_species['scientificName'].str.contains('SH')]

    if len(shs) > 0:
        shs_list = shs['acceptedScientificName'].unique()
        
        nest_asyncio.apply()

        shs_code = []
        shs_binomial_author = []
    
        asyncio.run(main(shs_list, shs_code, shs_binomial_author, wait_time=1.5))

        sh_dict = {'code':shs_code, 'scientificName':shs_binomial_author}
        sh_api = pd.DataFrame(sh_dict)

        shs.drop('scientificName', axis=1, inplace=True)
        shs = pd.merge(shs, sh_api, how='left', left_on='acceptedScientificName', right_on='code')
        shs.drop('code', axis=1, inplace=True)
        shs = shs[df_species.columns]

        logging.info(f"SHs processed: {len(shs)}")

        return shs
    
    else:
        logging.info('No SHs detected.')

        return shs

def join_df_shs(df_species, shs):
    """
    Simply merge non-SH and SH data together.
    """

    logging.info("Joining SH and non-SH occurrences")

    df_species['scientificName'] = [re.sub(r'SH.FU', '', i) for i in df_species['scientificName']]
    df_species = df_species[df_species.scientificName != '']
    df_species = pd.concat([df_species, shs], axis=0, ignore_index=True)
    df_species.reset_index(inplace=True)

    logging.info(f"Total records after join (plus second filter for empty scientificNames): {len(df_species)}")

    return df_species

def perform_harmonisation(df_species, mb_path):
    """
    Handle all taxonomic harmonisation steps, orchestrate data along the pipeline.
    """

    logging.info("Starting taxonomic harmonisation")

    mycobank = format_mb(mb_path)

    ## Exact match between names
    df_species['current_name'] = exact_matches(df_species, mycobank)

    ## Extract mismatches to pass to fuzzy match function
    mismatches = df_species[df_species['current_name']=='NA'].copy()
    mismatches.reset_index(inplace=True)

    ## Remove identified mismatches from object containing direct matches. Will be add later after fuzzy matching.
    harmonised = df_species[df_species['current_name']!='NA'].copy()
    harmonised.reset_index(inplace=True)

    logging.info(f"Exact matches found: {len(harmonised)} ({(len(harmonised)/len(df_species))*100:.1f}% of total occurrences)")

    ## Fuzzy matching between names
    mismatches['fuzzname'], mismatches['fuzzscore'], mismatches['current_name'] = fuzzy_match(mismatches, mycobank, score_threshold=jaro_threshold)
    fuzzymatched = mismatches[mismatches['current_name']!='NA']

    logging.info(f"Fuzzy matches (first iteration): {len(fuzzymatched)} ({(len(fuzzymatched)/len(df_species))*100:.1f}% of total occurrences)")

    fuzzymismatched = mismatches[mismatches['current_name']=='NA']

    ## Fuzzy matching with those names that did not match above. Now, selecting potential matches based on genus
    fuzzymismatched['fuzzname'],fuzzymismatched['fuzzscore'],fuzzymismatched['current_name'],fuzzymismatched['epithet_score'],fuzzymismatched['author_score']=fuzzy_match_genera(fuzzymismatched, mycobank)
    fuzymatchedfinal = fuzzymismatched[fuzzymismatched['fuzzscore']!='NA']
    fuzymatchedfinal['fuzzscore']=fuzymatchedfinal['fuzzscore'].astype(float)

    logging.info(f"Fuzzy matches (second iteration): {len(fuzymatchedfinal)} ({(len(fuzymatchedfinal)/len(df_species))*100:.1f}% of total occurrences)")

    ## Creating .csv with cases where epithet fuzzscore >= 0.07 to check manually
    manual_check = fuzymatchedfinal[fuzymatchedfinal['epithet_score'] >= 0.07]

    logging.info(f"Manual check matches: {len(manual_check)} ({(len(manual_check)/len(df_species))*100:.1f}% of total occurrences)")

    if occ_strict:
        manual_check.drop_duplicates(subset='scientificName').to_csv(output+'manual_check_taxonomy_strict.csv')
    else:
        manual_check.drop_duplicates(subset='scientificName').to_csv(output+'manual_check_taxonomy_relaxed.csv')

    ## Reading back verified names
    if occ_strict:
        filename = verified_manually_strict
    else:
        filename = verified_manually_relaxed

    while not os.path.exists(filename):
        print(f'{filename} not found.')
        action = input(f'Please create {filename} or check if it is accessible. [A]bort or [T]ry again? (A/T): ').strip().upper()

        if action == 'A':
            raise FileNotFoundError(f'File {filename} not found and user aborted')
        elif action == 'T':
            if os.path.exists(filename):
                break
    
    verified_names = pd.read_excel(filename)
    verified_names_dict = verified_names.set_index('scientificName')['current_name'].to_dict()

    manual_check['current_name'] = manual_check['scientificName'].map(verified_names_dict).combine_first(manual_check['current_name'])
    manual_check = manual_check[manual_check['current_name']!='NOTFOUND']

    ## Merging all dataframes
    fuzzymatchedfinal = fuzymatchedfinal[fuzymatchedfinal['epithet_score']<0.07]
    fuzzymatched = fuzzymatched.drop(['fuzzname', 'fuzzscore'], axis=1)
    fuzzymatchedfinal = fuzzymatchedfinal.drop(['fuzzname', 'fuzzscore', 'epithet_score', 'author_score'], axis=1)

    manual_check = manual_check.drop(['fuzzname', 'fuzzscore', 'epithet_score', 'author_score'], axis=1)

    occurrences_harmonised = pd.concat([harmonised, fuzzymatched, fuzzymatchedfinal, manual_check])
    occurrences_harmonised = pd.concat([harmonised, fuzzymatched, fuzzymatchedfinal, manual_check])
    occurrences_harmonised = occurrences_harmonised.reset_index(drop=True)

    derep_occs = occurrences_harmonised.drop_duplicates(subset='current_name')

    spp_n = occurrences_harmonised['current_name'].nunique()

    gen = occurrences_harmonised['current_name'].str.split().str[0]
    gen_uni = gen.unique()

    logging.info(f"Total species harmonised (known and accepted total species estimate): {spp_n}")
    logging.info(f"Total genera harmonised (known and accepted total genera estimate): {len(gen_uni)}")

    if occ_strict:
        occurrences_harmonised.to_csv(output+'occurrences_harmonised_strict.csv')
        derep_occs.to_csv(output+'derep_occs_strict.csv')
    else:
        occurrences_harmonised.to_csv(output+'occurrences_harmonised_relaxed.csv')
        derep_occs.to_csv(output+'derep_occs_relaxed.csv')
