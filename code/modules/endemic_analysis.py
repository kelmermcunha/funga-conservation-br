## Functions to perform endemic species inferences presented in the manuscript entitled "Challenges and Perspectives for Fungal Conservation in Brazil"
## Authors: Domingos Cardoso & Kelmer Martins-Cunha
## Corresponding author: Kelmer Martins-Cunha (kelmermartinscunha@gmail.com)


def find_endemics(country_sets, binomials):
    """
    Identify endemic species restricted to Brazil (country code 'BR').
    """

    br_endemics = []

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

    return br_endemics