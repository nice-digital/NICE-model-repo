# Model documentation

This folder contains documentation for the Exeter Oncology Model: Renal Cell Carcinoma edition (EOM:RCC). It is structured as a **quarto website** and hosted on **GitHub pages**, available at: <https://pythonhealthdatascience.github.io/stars-eom-rcc/>

## View documentation locally

To view this quarto site locally you will need to:

1. **Ensure Quarto and R are installed**
2. **Restore the renv** - if you have not used renv previously, you will first need to install it using `install.packages("renv")`. You can then run `renv::restore()`, which will create an R environment based on the instructions provided in the `renv.lock` file. You will need to run this from the parent folder, which contains the renv.
3. **Run the model and save the results** - this is required by the quarto site, as the documentation uses pre-run results, but the result file is too large to be synced with GitHub. For further instructions, see `pages/walkthrough/README.md`.
4. **Build the quarto site** - you can then use `quarto::render()` to build the quarto site

## Repository overview

```bash
├── images
│   └──  ...
├── pages
│   └──  ...
├── pdf
│   └──  ...
├── scripts
│   └──  ...
├── .gitignore
├── README.md
├── _quarto.yml
├── elsvan.csl
├── index.qmd
├── references.bib
└── styles.css
```

* `images/` - image files used in README or quarto site
* `pages/` - pages of the quarto site, except the home page
* `pdf/` - PDF files used in the quarto site
* `scripts/` - contains a single script, used to generate some illustrative plots for the quarto site
* `.gitignore` - untracked files
* `README.md` - this file!
* `_quarto.yml` - set-up instructions for the quarto site
* `elsvan.csl` - instructions for formatting of citations
* `index.qmd` - quarto site home page
* `references.bib` - list of BibTeX-formatted references
* `styles.css` - css styling used to alter banner height of quarto site to accomodate logo