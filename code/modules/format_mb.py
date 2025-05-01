## Auxiliary function for Precleaning data for all analysis presented in the manuscript entitled 
## "Challenges and Perspectives for Fungal Conservation in Brazil"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Corresponding author: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)

import re
import pandas as pd

def format_mb(path):
    mycobank = pd.read_excel(path)

    keys = ['gen.', 'fam.', 'ordo', 'subgen.', 'subfam.', 'sect.', 'tr.', 'subsect.', 'subcl.',
    'subdiv.', 'cl.', 'div.', 'ser.', 'subdivF.', 'regn.', 'subregn.', 'subordo', 'subtr.',
    'stirps', 'subser.']

    joint_keys = '|'.join(map(re.escape, keys))
    mycobank = mycobank[~mycobank['Rank.Rank name'].str.contains(joint_keys, case=False, na=False)]

    suffixes = ['Fungi','carpi','ceae','cei','ceteae','cetes','da','dae','dea','deae','dei','dia','dyeae','eae','eta',
    'etes','fungi','gae','gastri','ida','ini','lea','leae','les','nae','nea','neae','nei','ora','osa','ota','phoria',
    'pila','pori','rea','ria','riae','rix','rya','sea','sta','tae','tea','teae','tei','tes','theae','tia','tina',
    'tres','tria','xa','xea','xia','zeae','zoa']    

    pattern = r'^(?:\S+)$'

    mycobank['Taxon name'] = mycobank['Taxon name'].str.strip()
    mycobank = mycobank[~(
        mycobank['Taxon name'].str.match(pattern, na=False) &  # Ensures it's a single word
        mycobank['Taxon name'].str.endswith(tuple(suffixes), na=False)  # Checks if it ends with a suffix
    )]

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


    current_names = []
    for i in mycobank['Synonymy']:
        current_names.append(extract_current_name(i))

    for i in range(len(current_names)):
        if 'Current name:' in current_names[i]:
            pat = r"(?<=Current name: ).*?(?=\(\?\)|\{\?\}|:|\(\d+\))"
            treated = re.search(pat, current_names[i])
            if treated:
                current_names[i] = treated.group()      
            else:
                current_names[i] = mycobank['Synonymy'][i]
    
    current_names = [i.replace('(?)', '') for i in current_names]
    current_names = [i.replace('anon. ined.', '') for i in current_names]
    current_names = [re.sub(r"\(\d+\)", "", i) for i in current_names]
    current_names = [re.sub(r",\s*[^,]*\d+", "", i) for i in current_names]
    current_names = [re.sub(r":.*", "", i) for i in current_names]
    current_names = [re.sub(r"\d+", "", i) for i in current_names]

    mycobank['binomial_authors'] = current_names
    mycobank['Authors'] = mycobank['Authors'].fillna('')
    mycobank['binomial_authors_syn'] = mycobank['Taxon name']+' '+mycobank['Authors']

    return mycobank



