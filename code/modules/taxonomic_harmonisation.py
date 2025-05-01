## Functions to perform taxonomic harmonisation of occurrences presented in the results of the manuscript entitled 
## "Brazil as a global player in Fungal Conservation: A rapid shift from neglect to Action"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Contact: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


import unicodedata
from tqdm import tqdm
from rapidfuzz.process import extract
from rapidfuzz.distance import JaroWinkler


def exact_matches(occurrences, mycobank):
    """
    Use full scientific name (binomial + authorship) to gather exact matches between MycoBank entries and species names associated with 
    occurrences.
    """

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

def fuzzy_match(mismatches, mycobank, score_threshold=0.15):
    """
    For those full scientific names that were not matched directly in the exact match step, this functions performs a fuzzy matching based on the
    Jaro Winkler distance, based on all possible matches.
    
    The score_threshold parameter is setted to 0.15 based on tests and it refers to a threshold for the Jaro Winkler distance to
    define which occurrences full scientific names will be flagged as NA and passed to the second fuzzy matching iteration.
    """

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
                    if score >= score_threshold:
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
    """
    Perform a second iteration of fuzzy matching for those full scientific names that have a Jaro Winkler distance >= 0.15 (defined above) and
    those that did not match. Use a subset of MycoBank entries based on genus. Perform fuzzy matching also on epithet and authorship.
    """

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