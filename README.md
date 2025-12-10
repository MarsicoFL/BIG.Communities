# BIG Communities

Data and code for the manuscript **“Identity‑by‑descent captures shared environmental factors at biobank scale.”**  
Pre‑print: <https://doi.org/10.1101/2025.05.03.652048>

This repository is part of the **Biorepository, Informatics & Genomics (BIG) initiative at the University of Tennessee Health Science Center**.  

Project home: <https://uthsc.edu/cbmi/big/>

Data resource publication: <https://www.nature.com/articles/s41467-025-59375-0>

**Live dashboard:** <https://francomarsico.shinyapps.io/BIG_Communities/>

## Highlights
- Identity‑by‑descent (IBD) clustering at biobank scale, summarising four continent‑level communities and seventeen sub‑communities.
- Aggregated environmental indicators (for example PM2.5) aligned to the same ZIP code geography used for the communities.
- Interactive Shiny application to compare communities, exposures, and health outcomes side by side.

## Contents
- `app.r` – Shiny dashboard that ties together the data layers and visualisations.
- `data/` – Aggregated data products used by the dashboard (no individual‑level records).
- `rsconnect/` – Deployment metadata for shinyapps.io.

## Getting started locally
1. Install R (≥ 4.3) and RStudio or your preferred IDE.
2. Install the required packages:

```r
install.packages(c(
  "shiny", "shinyWidgets", "shinydashboard", "shinydashboardPlus",
  "leaflet", "leaflet.extras", "sf", "geojsonio", "dplyr", "tidyr",
  "readr", "DT", "ggplot2", "viridisLite", "plotly", "tibble",
  "forcats", "stringr", "RColorBrewer"
))
```

3. Launch the dashboard from the repository root:

```r
shiny::runApp(".")
```

The app reads the pre‑computed objects in `data/` by default. No credentials are required.

## Model
The joint analysis of genetic (G) and environmental factors (E) and their impact on health outcomes (H) enables probabilistic modeling of causal pathways. The diagram illustrates a gene-environment mediation framework where genetic factors G influence health outcome H through direct effects and indirect effects mediated by environmental factors E, which comprise both measured (Em) and unmeasured (Eu) components. The indirect pathway reflects the phenomenon of fine-scale population structure, whereby individuals sharing greater identity-by-descent (IBD) segments also tend to share environmental exposures due to geographic and social proximity. This framework facilitates the decomposition of total genetic effects into direct and environmentally-mediated components, enabling identification of potentially modifiable pathways in the genotype-phenotype relationship.

<img width="255" height="217" alt="model" src="https://github.com/user-attachments/assets/06572967-5661-4971-bc45-94c45470073f" />


## Data notes
- All data are aggregated to protect participant privacy; counts are filtered so that ZIP codes with fewer than 100 individuals are removed from the map.
- Environmental layers (for example PM2.5) come from public sources referenced in the manuscript; see the pre‑print for complete methodology.
- If you need to regenerate any tables, follow the scripts bundled in `data/` or reach out via the contact below.

## Screenshots

### Community overview
<img width="1909" height="957" alt="Community overview" src="https://github.com/user-attachments/assets/315b5737-8dfd-40a8-a591-d10583791a2e" />

### Mapping BIG
<img width="1909" height="957" alt="Geospatial comparison" src="https://github.com/user-attachments/assets/65221633-06b4-4abf-8fc5-b6222b07db6d" />

### Health conditions across communities
<img width="1909" height="957" alt="Health comparison" src="https://github.com/user-attachments/assets/6642f5c3-aadd-4cf4-a535-0f6c235f1cf5" />

### Communities, environment, and health
<img width="1909" height="957" alt="Communities, environment, and health" src="https://github.com/user-attachments/assets/7a1109d6-5383-4d60-8497-a7af2b1e7ac2" />

## Citation
If you use this code or dashboard, please cite the pre‑print above and acknowledge the BIG Initiative.

## Contact
Questions or feedback: fmarsico@uthsc.edu
