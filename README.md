# Exeter Oncology Model: RCC edition produced as part of NICE's pathways pilot

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

  If you re-use this model please appropriately credit PenTAG for the work and refer to it as the Exeter Oncology Model: RCC edition; Authors: Dawn Lee, Madhusubramanian Muthukumar, Alan Lovell, Caroline Farmer, Darren Burns, Justin Matthews, Helen Coelho, Brian O’Toole, Laura A Trigg, Tristan Snowsill, Maxwell S Barnish, Thalia Nikoglou, Amanda Brand, Zain Ahmad, Ahmed Abdelsabour, Sophie Robinson, Edward CF Wilson, G.J. Melendez-Torres.
