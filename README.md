# Providing Adaptive knowledge for Ratcheting up the EU Biodiversity strategy for Sustainable landscapes and protected areas
## PAREUS


The aim of PAREUS is to perform a spatially explicit cross-country assessment of existing land-use planning and conservation practices to identify improvements needed for reaching the EU Biodiversity Strategy 2030 targets. Specifically, PAREUS will apply an innovative landscape approach to integrate policy and practice for multiple land uses to ensure equitable and sustainable use of land within a spatially explicit framework and in close collaboration with stakeholders from multiple sectors to ensure the co-creation of knowledge. 
Within the framework of the project, PAREUS will create an engagement and co-creation arena involving the whole spectrum of stakeholder perspectives with often conflicting social and economic interests regarding PAs, OECM and sustainable landscapes. Jointly, the project partners will explore different ways to create a coherent network of PAs for biodiversity where also the wider countryside is included. For this, PAREUS will perform a spatially explicit social-ecological accounting of regional landscapes to enhance the potential for biodiversity conservation through OECM approaches. We will also evaluate the policy landscape regarding the consequences of current land-use planning practice. A geo-prospective planning tool will be developed to create participatory visualizations of the future orientation of landscapes. Our communication strategy is to create increased awareness about the possibility to achieve biodiversity targets through a combination of land sparing (PAs) and land sharing (OECM) for a broad target audience. Through transnational and transdisciplinary research and national case studies, PAREUS will identify factors that prevent targets from being implemented at the national scale, and how they may be crystallized into national and regional strategies and policies. By tackling the central challenge of biodiversity and land degradation from the perspective of transformative governance, PAREUS will be able to complement PA conservation with OECM approaches and help EU member states and associate members to implement their own national biodiversity strategy.

This repository provides the documentation and navigation structure for PAREUS. The project is organised into six work packages (WP1–WP6). Each work package contains reproducible code, documentation, metadata, and links to archived datasets.

Research data are **not stored in this repository**. Instead, all datasets are archived in a trusted data repository and referenced here using persistent identifiers (DOIs).

---

## Objectives

- Enhance stakeholder engagement, participation, and co-creation for improved understanding and adaptive implementation of transformative protected and conserved areas for sustainable land-use planning.
- Map and account the social-ecological landscape to quantify the current state of the network of protected and conserved areas in the wider countryside.
- Assess the policy landscape to identify opportunities and hindrances of existing legislative and policy documents for protected and conserved areas in the landscape.
- Develop an integrated and inclusive geoprospective tool to adaptively visualize and design coherent networks of protected and conserved areas in real landscapes.
- Synthesize the insights acquired from the regional social- ecological and policy mapping for exploring pathways for upscaling towards protected and conserved areas in the wider countryside.
- Operationalize the integrated landscape approach consistently in the different national case study regions in Europe.

---

## Project structure

```

PROJECT/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── project_metadata.yml
├── data_catalogue.csv
│
├── code/
│   ├── WP1/
│   ├── WP2/
│   ├── ...
│
├── data/
│   ├── WP1/
│   ├── WP2/
│   ├── ...
│
├── docs/
│   ├── WP1/
│   ├── WP2/
│   ├── ...
│
└── outputs/
│   ├── WP1/
│   ├── WP2/
│   ├── ...
└── Shared resources
```
Code: Contains source code for processing and analysing data, organised by work package.
Data: Contains metadata and links to archived datasets.
Docs: Documentation for each work package including method, workflow and processes
Outputs: Project outputs organised by work package such as maps, figures and cookbooks
---

# Work packages

| WP | Description | Lead | Repository | Data | Documentation |
|----|-------------|------|------------|------|---------------|
| WP1 | Co-creating transformative PCA landscapes | Sigrid Engen (NINA) | GitHub | DOI | Link |
| WP2 | Mapping and accounting of social-ecological PCA landscapes | Philip Roche (INRAE) | GitHub | DOI | Link |
| WP3 | Evaluating policies for PCA landscapes | Robert Kanka (ILESAS) | GitHub | DOI | Link |
| WP4 | Developing the geoprospective PCA landscape tool| Reto Spielhofer (NINA) | GitHub | DOI | Link |
| WP5 | Synthesizing PCA landscapes for sustainability | Philip Roche (INRAE) | GitHub | DOI | Link |
| WP6 | Coordination of National Case Studies | Robert Kanka (ILESAS) | GitHub | DOI | Link |

---


# Shared resources

Datasets

- National DEM
- Land cover
- Administrative boundaries
- Climate data

Software

- Shared R packages
- Shared Python utilities

Standards

- Coordinate reference system
- Naming conventions
- Metadata standard
- Coding guidelines

---

## Data access

Research datasets are archived in a trusted repository and assigned persistent DOIs.

| WP | Repository | DOI |
|----|------------|-----|
| WP1 | Zenodo | DOI |
| WP2 | Dataverse | DOI |
| WP3 | Zenodo | DOI |
| WP4 | Institutional repository | DOI |

Detailed metadata and download instructions are provided in the corresponding `data/WP*/README.md`.

---

## Reproducibility

Each work package includes:

- documented workflows
- version-controlled code
- software dependencies (`renv`)
- metadata
- archived datasets
- reproducible outputs

---

## Licence

- Code: MIT
- Documentation: CC BY 4.0
- Data: see individual dataset licences.

---

## Citation

Please cite:

1. the project,
2. the relevant work package,
3. the archived dataset,
4. and the associated publication where applicable.

---

# Citation

If using this project, please cite

- project software
- relevant work package
- archived dataset
- associated publication

Citation files (`CITATION.cff`) are available in each repository.

---

# Licence

Code: MIT

Documentation: CC-BY 4.0

Datasets: see individual dataset licences.

---

# Contact

Project lead

Roel May

Norwegian Institute for Nature Research (NINA)

roel.may@nina.no

Technical contact

Reto Spielhofer

reto.spielhofer@nina.no

---