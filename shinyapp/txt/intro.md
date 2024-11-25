The ambition for EOM:RCC is to create a shiny interface that allows users to interact and run the full economic model from a web application.

This is a **pilot example** of just one part of the analysis, which is to create a **table of possible treatment sequences**. This will vary depending on the treatments being considered, the included populations, maximum number of treatment lines. The included populations are:

* **Population 1:** Favourable risk and more than 1 year since immuno-oncology (IO) treatment
* **Population 2:** Intermediate/poor risk and more than 1 year since IO
* **Population 3:** Favourable risk and less than 1 year since IO
* **Population 4:** Intermediate/poor risk and less than 1 year since IO

The "favourable" and "intermediate/poor" risk statuses are as defined by the International Metastatic Renal Cell Carcinoma Database Consortium (IMDC).

The app has two tabs:

* **Run analysis:** To create a table of possible treatment sequences, as described above.
* **Sequence explorer:** Select first, second, third and fourth line treatments, and see how the treatments available at each line vary depending on your choices.

This app has been developed by [Amy Heather](https://github.com/amyheather) [![ORCID 0000-0002-6596-347](../www/ORCIDiD_icon16x16.png)](https://orcid.org/0000-0002-6596-3479) as part of the project STARS (Sharing Tools and Artifacts for Reusable Simulations in healthcare). The app code can be viewed at <https://github.com/nice-digital/NICE-model-repo/tree/main/shinyapp>.

Other resources related to this model include:

* [Model documentation](https://nice-digital.github.io/NICE-model-repo): A website provides detailed documentation for the model. Among other things, this includes a detailed and plain english summary of the analysis, installation instructions, a step-by-step code walkthrough, and descriptions of the probabilistic analysis and scenario analysis.
* [GitHub repository](https://github.com/nice-digital/NICE-model-repo): All code for the model, documentation and app are contained in this repository. The `README.md` file (which is displayed when you open the repository on GitHub) provides key information about the model and repository.