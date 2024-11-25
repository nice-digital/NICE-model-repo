# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates formatted as YYYY-MM-DD as per [ISO standard](https://www.iso.org/iso-8601-date-and-time-format.html).

## v1.1.0 - 2024-11-25

Merged modifications by Amy Heather made as part of the project [STARS (Sharing Tools and Artefacts for Reproducible Simulations in healthcare)](https://github.com/pythonhealthdatascience/stars-eom-rcc). These changes included improvements to the README, creating detailed documentation hosted via GitHub pages, and creation of a pilot shiny web application to find valid treatment sequences (including an interactive sequence explorer).

### Added

* Documentation (`docs/`) hosted on GitHub pages
* Web application (`shinyapp/`) hosted on shinyapps.io
* `CHANGELOG.md` to provide reader-friendly explanation of differences between versions
* `CITATION.cff` to provide instructions on how to cite the repository

### Changed

* Extended and improved details in README
* Add PenTAG email to licence
* Add `folder` input to `f_res_ProduceWordDoc()` to allow specification of the output folder to save the word document too (as previously hard coded `./4_Output`)

### Fixed

* Formatting of copyright statement in `LICENSE`

### Removed

* Removed unused template README and citation files from `1_Data/` `3_Functions/` and `4_Output/`. The paragraph that had been entered into the template README in the data folder was moved onto the probabilistic analysis page in the quarto documentation

## v1.0.0 - 2024-07-30

🌱 First release of EOM-RCC model.