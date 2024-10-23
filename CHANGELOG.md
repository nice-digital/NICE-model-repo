# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates formatted as YYYY-MM-DD as per [ISO standard](https://www.iso.org/iso-8601-date-and-time-format.html).

## v1.4.0 - 2024-10-23

Add interactive sequence explorer to the shiny application, and made amendments elsewhere to the documentation, as per feedback from PenTAG team.

### Added

* Interactive sequence explorer tab within the shiny app
* GitHub action to convert `CITATION.cff` into a `.bib` citation that can easily download and use with reference managers

### Changed

* Amendments to the text and diagrams across several different pages, as per feedback from the PenTAG team
* Add Amy Heather to citation
* Add PenTAG email to licence
* Add reader mode so sidebar and table of contents can be toggled to hide
* Improved image quality of the draw.io diagrams
* Update list of packages in renv for shiny app

### Removed

* Deleted old draw.io diagrams

### Fixed

* Corrected display of survival and hazard curves

## v1.3.0 - 2024-09-27

Produced pilot shiny web application to find valid treatment sequences.

### Added

* Shiny web application (in `shinyapp/`), which is also hosted on shinyapps.io

### Changed

* Add required packages for app to the `renv`
* Add description of app to the main `README` file

## v1.2.0 - 2024-09-23

Detailed model documentation via a quarto site hosted with GitHub pages.

### Added

* Quarto site with detailed documentation for the model including:
  * Acronyms
  * Context on the associated NICE appraisals, articles and reports
  * A detailed summary of the analysis
  * A plain english summary of the analysis
  * Installation instructions
  * A step-by-step walkthrough of the code in `Model_Structure.R`
  * Descriptions of the probabilistic analysis and scenario analysis
  * Details about the license, citation instructions and the changelog

### Changed

* Add `folder` input to `f_res_ProduceWordDoc()` to allow specification of the output folder to save the word document too (as previously hard coded `./4_Output`)
* Add information about quarto website and run time to `README.md`
* Moved images from `img/` to `docs/images/`

### Removed

* Removed unused template README and citation files from `1_Data/` `3_Functions/` and `4_Output/`. The paragraph that had been entered into the template README in the data folder was moved onto the probabilistic analysis page in the quarto documentation

## v1.1.1 - 2024-08-20

Embellished README

### Added

* Add files produced by model run to `.gitignore`

### Changed

* Extended `README.md` to include banner image with logos, extra badges, table of contents, all the information from `ID6184 Using the R decision model..` (e.g. installation guide with images, overview of input files, future versions)

## v1.1.0 - 2024-08-16

Implemented the essential components of the STARS framework (exc. open science archive).

### Added

* `CITATION.cff` and GitHub action to check validity (`cff_validation.yaml`)

### Changed

* Extended `README.md` to include some instructions for installing and running the model, more detailed repository overview, citation information, ORCID IDs, acknowledgements, license and funding information

### Fixed

* Formatting of copyright statement in `LICENSE`

## v1.0.0 - 2024-08-16

🌱 First release. EOM-RCC model as shared by the PenTAG team, with minor changes that enabled us to run the model.

### Added

* Code from the original repository - <https://github.com/nice-digital/NICE-model-repo>
* R environment for dependency management (`renv`)

### Changed

* Set model to run sequentially (in `2_Scripts/Model_Structure.R`)
