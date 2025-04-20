This repository contains all the code used in the results presented in the manuscript "Brazil as a global player in Fungal Conservation: A rapid shift from neglect to Action". It is organized into two main parts:

1. `code` directory: all scripts used to pre-clean data, prepare databases and perform spatial and endemicity analyses to estimate the total accepted known and accessible fungal species occurring in Brazil and infer the percentage of putatively endemic species within retrieved names. The README file in this subdirectory contains further details about the scripts.
2. The `data` directory contains data needed to run scripts and outputs that the pipeline generates. It is also used for downstream analyses and figures.

Scripts are written in Python and interact with the PlutoF API (https://plutof.docs.apiary.io/#) to retrieve fungal Species Hypothesis data, and with th Global Biodiversity Information Facility (GBIF —https://techdocs.gbif.org/en/openapi/) API to retrieve usage Keys and send occurrence download requests.
