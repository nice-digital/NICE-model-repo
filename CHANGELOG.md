# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates formatted as YYYY-MM-DD as per [ISO standard](https://www.iso.org/iso-8601-date-and-time-format.html).

## v1.1.0 - 2024-08-16

Implemented the essential components of the STARS framework (exc. open science archive).

### Added

* `CITATION.cff` and GitHub action to check validity (`cff_validation.yaml`)

### Changed

* Extended `README.md` to include some instructions for installing and running the model, more detailed repositoriy overview, citation information, ORCID IDs, acknowledgements, license and funding information

### Fixed

* Formatting of copyright statement in `LICENSE`

## v1.0.0 - 2024-08-16

🌱 First release. EOM-RCC model as shared by the PenTAG team, with minor changes that enabled us to run the model.

### Added

* Code from the original repository - <https://github.com/nice-digital/NICE-model-repo>
* R environment for dependency management (`renv`)

### Changed

* Set model to run sequentially (in `2_Scripts/Model_Structure.R`)
