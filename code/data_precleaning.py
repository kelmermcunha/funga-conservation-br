## Precleaning data for all analysis presented in the manuscript entitled 
## "Brazil as a global player in Fungal Conservation: A rapid shift from neglect to Action"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Contact: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


# Loading modules
import re
import asyncio
import nest_asyncio
import unicodedata
import aiohttp
import pandas as pd
from tqdm import tqdm
from format_mb import *
from rapidfuzz.distance import JaroWinkler
from rapidfuzz.process import extract



def extract_current_name(text):
    
    current_name_string = re.search(r"Current name: .*?(?= synonym(s)?|\[|Basionym)", text)
    if current_name_string:
        current_name = current_name_string.group()
    else:
        current_name_string = re.search(r"^[A-Z][a-z]+.*?(?= synonym(s)?|\[|Basionym)", text)
        if current_name_string:
            current_name = current_name_string.group()
        else:
            return text
        
    excluded_words = ['champignons', 'Floræ Scandinaviæ', 'Fungi Europaei']
    if "&" in current_name and all(word not in current_name for word in excluded_words):
        pattern = r"(?<=Current name: ).*?&.*?(?=,)"
        match = re.search(pattern, current_name)
        if match:
            return match.group()
        else:
            return current_name
    else:
        pattern = r"(?<=Current name: ).*?(?=,)"
        match = re.search(pattern, current_name)
        if match:
            return match.group()
        else:
            pattern = r"(?<=Current name: ).*?(?=\[|:)|^(?!Current\s)[A-Z][a-z]+.*?(?=\[|:)"
            match = re.search(pattern, text)
            if match:
                return match.group()
            else:
                return current_name
def exact_matches(occurrences, mycobank):
    
    occ_species = occurrences['scientificName'].tolist()
    occ_species = [name.upper().replace(' ', '') for name in occ_species]
    
    mycobank_names = mycobank['binomial_authors_syn'].tolist()
    mycobank_names = [name.upper().replace(' ', '') for name in mycobank_names]
    
    mycobank_binomial_authors = mycobank['binomial_authors'].tolist()
    mycobank_name_status = mycobank['Name status'].tolist()
    
    mb_dict = {}
    for i, taxon in enumerate(mycobank_names):
        if taxon not in mb_dict:
            mb_dict[taxon] = []
        mb_dict[taxon].append((mycobank_binomial_authors[i], mycobank_name_status[i]))
    
    current_names = []
    for index, taxon in enumerate(occ_species):
        
        if taxon in mb_dict:
            matches = mb_dict[taxon]
            if len(matches) > 1:
                legit=False
                for i in matches:
                    if 'legitimate' in i[1].lower():
                        current_names.append(i[0])
                        legit=True
                        break
                if not legit:
                    current_names.append(matches[0][0])
            else:
                current_names.append(mb_dict[taxon][0][0])
        else:
            current_names.append('NA')
    
    return current_names

def fuzzy_match(mismatches, mycobank):
    mismatches_names = mismatches['scientificName'].tolist()
    
    mycobank_taxon_names = mycobank['Taxon name'].tolist()
    mycobank_binomial_authors_syn = mycobank['binomial_authors_syn'].tolist()
    mycobank_binomial_authors = mycobank['binomial_authors'].tolist()

    mb_dict = {}
    for i, taxon in enumerate(mycobank_taxon_names):
        if taxon not in mb_dict:
            mb_dict[taxon] = []
        mb_dict[taxon].append((mycobank_binomial_authors_syn[i], mycobank_binomial_authors[i]))
    
    fuzznames = []
    fuzzscores = []
    current_names = []
    
    for index, species in enumerate(mismatches_names):
        if ' var. ' in species:
            fuzznames.append('NA')
            fuzzscores.append('NA')
            current_names.append('NA')
        else:
            sp = ' '.join(species.split()[:2])
            if sp in mb_dict:
                dict_entries = mb_dict[sp]
                binomial_synonyms = [entry[0] for entry in dict_entries]
                try:
                    match = extract(species, binomial_synonyms, scorer=JaroWinkler.distance)[0]
                    fuzzname = match[0]
                    score = match[1]
                    if score >= 0.15:
                        fuzznames.append('NA')
                        fuzzscores.append('NA')
                        current_names.append('NA')
                    else:
                        fuzznames.append(fuzzname)
                        fuzzscores.append(score)

                        for entry in dict_entries:
                            if entry[0]==fuzzname:
                                current_names.append(entry[1])
                                break
                except Exception:
                    fuzznames.append('NA')
                    fuzzscores.append('NA')
                    current_names.append('NA')
            else:
                fuzznames.append('NA')
                fuzzscores.append('NA')
                current_names.append('NA')

    return fuzznames, fuzzscores, current_names
def fuzzy_match_genera(mismatches, mycobank):
    mismatches_names = mismatches['scientificName'].tolist()
    
    mycobank_taxon_names = mycobank['Taxon name'].tolist()
    mycobank_binomial_authors_syn = mycobank['binomial_authors_syn'].tolist()
    mycobank_binomial_authors = mycobank['binomial_authors'].tolist()

    mb_dict = {}
    for i, taxon in enumerate(mycobank_taxon_names):
        if taxon not in mb_dict:
            mb_dict[taxon] = []
        mb_dict[taxon].append((mycobank_binomial_authors_syn[i], mycobank_binomial_authors[i]))

    mb_dict = {
        unicodedata.normalize('NFKD', k).encode('ASCII', 'ignore').decode(): v
        for k, v in mb_dict.items()
    }
    
    fuzznames = []
    fuzzscores = []
    current_names = []
    epithet_scores = []
    author_scores = []
    
    # Add tqdm to monitor progress
    for index, species in tqdm(enumerate(mismatches_names), desc="Matching species", total=len(mismatches_names)):
        gen = unicodedata.normalize("NFKD", ' '.join(species.split()[:1])).encode("ASCII", "ignore").decode()
        
        gen_matches = []
        for key in mb_dict:
            if key.startswith(gen):
                gen_matches.extend(mb_dict[key])

        binomial_synonyms = [entry[0] for entry in gen_matches]

        try:
            if any(substring in species for substring in [' var. ', ' .f ', ' .subsp. ', ' subgen. ', ' sect. ']):
                sep = species.split()
                epivar = ' '.join([sep[1], sep[2], sep[3]])
                aut = ' '.join(sep[4:])

                match = extract(species, binomial_synonyms, scorer=JaroWinkler.distance)[0]
                fuzznames.append(match[0])
                fuzzscores.append(match[1])

                bi_syn_sep = match[0].split()
                bi_syn_epivar = ' '.join([bi_syn_sep[1], bi_syn_sep[2], bi_syn_sep[3]])
                bi_syn_aut = ' '.join(bi_syn_sep[4:])

                match1 = JaroWinkler.normalized_distance(epivar, bi_syn_epivar)
                epithet_scores.append(match1)
                match2 = JaroWinkler.normalized_distance(aut, bi_syn_aut)
                author_scores.append(match2)

                for entry in gen_matches:
                    if entry[0]==match[0]:
                        current_names.append(entry[1])
                        break
            else:
                sep = species.split()
                epi = sep[1]
                aut = ' '.join(sep[2:])

                match = extract(species, binomial_synonyms, scorer=JaroWinkler.distance)[0]
                fuzznames.append(match[0])
                fuzzscores.append(match[1])

                bi_syn_sep = match[0].split()
                bi_syn_epi = bi_syn_sep[1]
                bi_syn_aut = ' '.join(bi_syn_sep[2:])

                match1 = JaroWinkler.normalized_distance(epi, bi_syn_epi)
                epithet_scores.append(match1)
                match2 = JaroWinkler.normalized_distance(aut, bi_syn_aut)
                author_scores.append(match2)

                for entry in gen_matches:
                    if entry[0]==match[0]:
                        current_names.append(entry[1])
                        break

        except Exception as e:
            fuzznames.append('NA')
            fuzzscores.append('NA')
            current_names.append('NA')
            epithet_scores.append('NA')
            author_scores.append('NA')
            

    return fuzznames, fuzzscores, current_names, epithet_scores, author_scores

# Loading data and databases

## GBIF occurrences — choose between datasets to run each scenario presented in the manuscript. Remember to change bellow accordingly.
occ_strict = './data/gbif/occurrence_strict.csv'
#occ_relaxed = './data/gbif/occurence.txt'

print('Loading GBIF data')
if occ_strict:
    rawData = pd.read_csv(occ_strict, on_bad_lines='warn', low_memory=False)
else:
    rawData = pd.read_csv(occ_relaxed, on_bad_lines='warn', low_memory=False)
                        ## Data was downloaded on 25/03/2023 from the GBIF database 
                        ## website. The following filters were used: 'Country or area' = Brazil,
                        ## 'Scientific name' = Fungi. Reference:
                        ## GBIF.org (21 March 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.pnxaj2

## Remove occurrences not identified at the species level, select relevant columns, and drop duplicates based on the selected columns
df_species = rawData.dropna(subset=['species']).copy()
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
            'verbatimLocality','year', 'month', 'day','verbatimLocality'], inplace=True)


## MycoBank
print('Loading MycoBank')
mycobank = format_mb('./data/MBList_2025_2.xlsx')

df_species['scientificName'] = [re.sub(r'\d+', '', i) for i in df_species['scientificName']]

## Treating Species Hypothesis (only for relaxed total species estimative)
shs = df_species[df_species['scientificName'].str.contains('SH')]
if len(shs) > 0:
    shs_list = shs['acceptedScientificName'].unique()
    
    nest_asyncio.apply()

    api_sh = "https://api.plutof.ut.ee/v1/public/dshclusters/search/?name="
    api_taxon = "https://api.plutof.ut.ee/v1/public/taxa/"

    shs_code = []
    shs_binomial_author = []

    async def fetch_sh_data(session, sh_code, retries=3):
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

    async def main(shs_list, wait_time=1.5):
        
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

                
    asyncio.run(main(shs_list, wait_time=1.5))

    sh_dict = {'code':shs_code, 'scientificName':shs_binomial_author}
    sh_api = pd.DataFrame(sh_dict)

    shs.drop('scientificName', axis=1, inplace=True)
    shs = pd.merge(shs, sh_api, how='left', left_on='acceptedScientificName', right_on='code')
    shs.drop('code', axis=1, inplace=True)
    shs = shs[df_species.columns]

    df_species['scientificName'] = [re.sub(r'SH.FU', '', i) for i in df_species['scientificName']]
    df_species = df_species[df_species.scientificName != '']
    df_species = pd.concat([df_species, shs], axis=0, ignore_index=True)
    df_species.reset_index(inplace=True)

# Taxonomic harmonisation

## Exact match between names
print('Exact match between GBIF and MycoBank taxonomy')
df_species['current_name'] = exact_matches(df_species, mycobank)

## Extract mismatches to pass to fuzzy match function
mismatches = df_species[df_species['current_name']=='NA'].copy()
mismatches.reset_index(inplace=True)

## Remove identified mismatches from object containing direct matches. Will be add later after fuzzy matching.
harmonised = df_species[df_species['current_name']!='NA'].copy()
harmonised.reset_index(inplace=True)

## Fuzzy matching between names

print('Fuzzy match between GBIF and MycoBank taxonomy')
mismatches['fuzzname'], mismatches['fuzzscore'], mismatches['current_name'] = fuzzy_match(mismatches, mycobank)
fuzzymatched = mismatches[mismatches['current_name']!='NA']
fuzzymismatched = mismatches[mismatches['current_name']=='NA']

## Fuzzy matching with those names that did not match above. Now, selecting potential matches based on genus
print('Fuzzy match second round using genera')
fuzzymismatched['fuzzname'],fuzzymismatched['fuzzscore'],fuzzymismatched['current_name'],fuzzymismatched['epithet_score'],fuzzymismatched['author_score']=fuzzy_match_genera(fuzzymismatched, mycobank)
fuzymatchedfinal = fuzzymismatched[fuzzymismatched['fuzzscore']!='NA']
fuzymatchedfinal['fuzzscore']=fuzymatchedfinal['fuzzscore'].astype(float)

## Creating .csv with cases where epithet fuzzscore >= 0.07 to check manually
manual_check = fuzymatchedfinal[fuzymatchedfinal['epithet_score'] >= 0.07]
if occ_strict:
    manual_check.drop_duplicates(subset='scientificName').to_csv('./data/outputs/manual_check_taxonomy_strict.csv')
else:
    manual_check.drop_duplicates(subset='scientificName').to_csv('./data/outputs/manual_check_taxonomy_relaxed.csv')

## Reading back verified names
if occ_strict:
    verified_names = pd.read_excel('./data/output/verified_names_manually_strict.xlsx')
else:
    verified_names = pd.read_excel('./data/output/verified_names_manually_relaxed.xlsx')

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

print(f'Total species names found: {spp_n}')

gen = occurrences_harmonised['current_name'].str.split().str[0]
gen_uni = gen.unique()

print(f'Total genera found: {len(gen_uni)}')

if occ_strict:
    occurrences_harmonised.to_csv('./data/output/occurrences_harmonised_strict.csv')
    derep_occs.to_csv('./data/output/derep_occs_strict.csv')
else:
    occurrences_harmonised.to_csv('./data/output/occurrences_harmonised_relaxed.csv')
    derep_occs.to_csv('./data/output/derep_occs_relaxed.csv')