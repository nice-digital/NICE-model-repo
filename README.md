# Exeter Oncology Model: RCC edition produced as part of NICE's pathways pilot

[![Valid CITATION.cff](https://github.com/pythonhealthdatascience/stars-eom-rcc/actions/workflows/cff_validation.yaml/badge.svg)](https://github.com/pythonhealthdatascience/stars-eom-rcc/actions/workflows/cff_validation.yaml)


The Exeter Oncology Model: RCC edition is a platform cost-effectiveness model encompassing each decision node in the disease area for advanced renal cell carcinoma.

This model has been created as part of a NICE pilot aimed both at reducing the long-term resource requirements for appraisal within crowded treatment pathways and at assessing the feasibility of incorporating treatment sequences into decision making.

The Exeter Oncology Model: RCC edition has been constructed as part of pathways pilot appraisal ID6186 and the appraisal of cabozantinib plus nivolumab ID6184. No data is contained in the code. All input data is contained in the data folder, dummy data is provided where the data used in the appraisal was marked as confidential. A user interface was originally planned be added to this model at a later stage in the project.

The folder structure loosely follows ZIN guidelines:

    1. Data: externally derived parameters should be saved here
    2. Scripts - one working script is provided (model_structure) which contains the following sections:
        i. Installation - containing all information to make the model operational. This section also states the version of R and the packages used at the time of submission.
        ii. Loading functions
        iii. Model inputs including: loading of input parameters, conduct of survival analysis, application of relative effectiveness and formatting of inputs for costs and quality of life
        iv. Population settings
        v. Patient flow: this is where the patient flow is produced dependent on the selected model structure and results are produced and compiled
    3. Functions
    4. Output (intermediate and final output data can be saved to here)

   A tests section had been planned for a later phase of the project.

  Instructions on how to use the model can be found in document: ID6184 Using the R decision model (EAG instructions) noACIC

## Citation

If you re-use this model please appropriately credit PenTAG for the work and refer to it as the Exeter Oncology Model: RCC edition:

> Lee D., Muthukumar M., Lovell A., Farmer C., Burns D., Matthews J., Coelho H., O'Toole B., Trigg L., Snowsill T., Barnish M., Nikoglou T., Brand A., Ahmad Z., Abdelsabour A., Robinson S., Wilson E., Melendez-Torres G. Exeter Oncology Model: RCC edition URL: https://github.com/nice-digital/NICE-model-repo

The author ORCID IDs (where available) are:

| Author | ORCID |
| - | - |
| Dawn Lee | [![ORCID: Lee](https://img.shields.io/badge/ORCID-0000--0003--4027--8456-brightgreen)](https://orcid.org/0000-0003-4027-8456) |
| Madhusubramanian Muthukumar | - |
| Alan Lovell | - |
| Caroline Farmer | - |
| Darren Burns | [![ORCID: Burns](https://img.shields.io/badge/ORCID-0000--0002--5209--8041-brightgreen)](https://orcid.org/0000-0002-5209-8041) |
| Justin Matthews | - |
| Helen Coelho | [![ORCID: Coelho](https://img.shields.io/badge/ORCID-0000--0002--4799--4300-brightgreen)](https://orcid.org/0000-0002-4799-4300) |
| Brian O’Toole | - |
| Laura Trigg | [![ORCID: Trigg](https://img.shields.io/badge/ORCID-0000--0002--8447--2616-brightgreen)](https://orcid.org/0000-0002-8447-2616) |
| Tristan Snowsill | [![ORCID: Snowsill](https://img.shields.io/badge/ORCID-0000--0001--7406--2819-brightgreen)](https://orcid.org/0000-0001-7406-2819) |
| Maxwell Barnish | [![ORCID: Barnish](https://img.shields.io/badge/ORCID-0000--0003--0139--6548-brightgreen)](https://orcid.org/0000-0003-0139-6548) |
| Thalia Nikoglou | - |
| Amanda Brand | - |
| Zain Ahmad | - |
| Ahmed Abdelsabour |  [![ORCID: Abdelsabour](https://img.shields.io/badge/ORCID-0009--0007--2532--4676-brightgreen)](https://orcid.org/0009-0007-2532-4676) |
| Sophie Robinson | [![ORCID: Robinson](https://img.shields.io/badge/ORCID-0000--0003--0463--875X-brightgreen)](https://orcid.org/0000-0003-0463-875X) |
| Edward Wilson | [![ORCID: Wilson](https://img.shields.io/badge/ORCID-0000--0002--8369--1577-brightgreen)](https://orcid.org/0000-0002-8369-1577) |
| G.J. Melendez-Torres | [![ORCID: Melendez-Torres](https://img.shields.io/badge/ORCID-0000--0002--9823--4790-brightgreen)](https://orcid.org/0000-0002-9823-4790) |

## Acknowledgements

This modified repository was developed by:

| Author | ORCID | GitHub |
| - | - | - |
| Amy Heather | [![ORCID: Heather](https://img.shields.io/badge/ORCID-0000--0002--6596--3479-brightgreen)](https://orcid.org/0000-0002-6596-3479) | <https://github.com/amyheather>