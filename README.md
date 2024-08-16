# Exeter Oncology Model: RCC edition produced as part of NICE's pathways pilot

[![Valid CITATION.cff](https://github.com/pythonhealthdatascience/stars-eom-rcc/actions/workflows/cff_validation.yaml/badge.svg)](https://github.com/pythonhealthdatascience/stars-eom-rcc/actions/workflows/cff_validation.yaml)

## About the model

The Exeter Oncology Model: RCC edition is a platform cost-effectiveness model encompassing each decision node in the disease area for advanced renal cell carcinoma.

This model has been created as part of a NICE pilot aimed both at reducing the long-term resource requirements for appraisal within crowded treatment pathways and at assessing the feasibility of incorporating treatment sequences into decision making.

The Exeter Oncology Model: RCC edition has been constructed as part of pathways pilot appraisal ID6186 and the appraisal of cabozantinib plus nivolumab ID6184. No data is contained in the code. All input data is contained in the data folder, dummy data is provided where the data used in the appraisal was marked as confidential.

A user interface was originally planned be added to this model at a later stage in the project. A tests section had also been planned for a later phase of the project.

## How to install and run the model and vary model parameters

Instructions on how to use the model can be found in the document: `ID6184 Using the R decision model (EAG instructions) noACIC.docx`.

## Repository overview

```bash
├── 1_Data
│   └──  ...
├── 2_Scripts
│   └──  ...
├── 3_Functions
│   └──  ...
├── 4_Output
│   └──  ...
├── renv
│   └──  ...
├── tests
│   └──  ...
├── .gitignore
├── .Rprofile
├── CHANGELOG.md
├── CITATION.cff
├── eom-rcc.Rproj
├── ID6184 Using the R decision model (EAG instructions) noACIC.docx
├── LICENSE
├── README.md
└── renv.lock
```

The folder structure for the model loosely follows the Zorginstituut Nederland (ZIN) "*Guideline for building cost-effectiveness models in R*":

* **Data:** Externally derived parameters should be saved here
* **Scripts:** One working script is provided (`Model_Structure.R`) which contains the following sections:
    * **Installation**: Containing all information to make the model operational. This section also states the version of R and the packages used at the time of submission.
    * **Loading functions**
    * **Model inputs** including: loading of input parameters, conduct of survival analysis, application of relative effectiveness and formatting of inputs for costs and quality of life
    * **Population settings**
    * **Patient flow**: This is where the patient flow is produced dependent on the selected model structure and results are produced and compiled
* **Functions**
* **Output:** Intermediate and final output data can be saved to here

The other files are folders are:

* `renv/`: Instructions for creation of R environment
* `tests/`: Tests
* `.gitignore`: Untracked files
* `.Rprofile`: Activates R environment
* `CHANGELOG.md`: Description of changes between GitHub releases
* `CITATION.cff`: Instructions for citing this repository
* `eom-rcc.Rproj`: R project settings
* `ID6184 Using the R decision model (EAG instructions) noACIC.docx`: Instructions for running the model
* `LICENSE`: MIT license
* `README.md`: This file!
* `renv.lock`: Lists R version and all packages in the R environment

## Citation

If you re-use this model please appropriately credit PenTAG for the work and refer to it as the Exeter Oncology Model: RCC edition:

> Lee D., Muthukumar M., Lovell A., Farmer C., Burns D., Matthews J., Coelho H., O'Toole B., Trigg L., Snowsill T., Barnish M., Nikoglou T., Brand A., Ahmad Z., Abdelsabour A., Robinson S., Wilson E., Melendez-Torres G. Exeter Oncology Model: RCC edition URL: https://github.com/nice-digital/NICE-model-repo

The author ORCID IDs (where available) are:

[![ORCID: Lee](https://img.shields.io/badge/Dawn_Lee-0000--0003--4027--8456-brightgreen)](https://orcid.org/0000-0003-4027-8456)
[![ORCID: Burns](https://img.shields.io/badge/Darren_Burns-0000--0002--5209--8041-brightgreen)](https://orcid.org/0000-0002-5209-8041)
[![ORCID: Coelho](https://img.shields.io/badge/Helen_Coelho-0000--0002--4799--4300-brightgreen)](https://orcid.org/0000-0002-4799-4300)
[![ORCID: Trigg](https://img.shields.io/badge/Laura_Trigg-0000--0002--8447--2616-brightgreen)](https://orcid.org/0000-0002-8447-2616)
[![ORCID: Snowsill](https://img.shields.io/badge/Tristan_Snowsill-0000--0001--7406--2819-brightgreen)](https://orcid.org/0000-0001-7406-2819)
[![ORCID: Barnish](https://img.shields.io/badge/Maxwell_Barnish-0000--0003--0139--6548-brightgreen)](https://orcid.org/0000-0003-0139-6548)
[![ORCID: Abdelsabour](https://img.shields.io/badge/Ahmed_Abdelsabour-0009--0007--2532--4676-brightgreen)](https://orcid.org/0009-0007-2532-4676)
[![ORCID: Robinson](https://img.shields.io/badge/Sophie_Robinson-0000--0003--0463--875X-brightgreen)](https://orcid.org/0000-0003-0463-875X)
[![ORCID: Wilson](https://img.shields.io/badge/Edward_Wilson-0000--0002--8369--1577-brightgreen)](https://orcid.org/0000-0002-8369-1577)
[![ORCID: Melendez-Torres](https://img.shields.io/badge/GJ_Melendez--Torres-0000--0002--9823--4790-brightgreen)](https://orcid.org/0000-0002-9823-4790)

## Acknowledgements

This modified repository was developed by [**Amy Heather**](https://github.com/amyheather) as part of work package 3 on the project "STARS: Sharing Tools and Artefacts for Reproducible Simulations". Changes from the original repository are described in the `CHANGELOG.md`.

[![ORCID: Heather](https://img.shields.io/badge/Amy_Heather-0000--0002--6596--3479-brightgreen)](https://orcid.org/0000-0002-6596-3479)

## License

This repository is licensed under an MIT license.

## Funding

The development of the EOM-RCC model for NICE, as part of the pathways pilot, was funded by the National Institute for Health and Care Research (NIHR) Evidence Synthesis Programme as project number [NIHR136008](https://www.dev.fundingawards.nihr.ac.uk/award/NIHR136008).

STARS is supported by the Medical Research Council [grant number [MR/Z503915/1](https://gtr.ukri.org/projects?ref=MR%2FZ503915%2F1)].