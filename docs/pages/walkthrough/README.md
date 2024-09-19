# README

The walkthrough pages in this folder are all created based on `code0a_walkthrough.qmd`.

To reduce the run time of this walkthrough, it loads **pre-run results for the economic model**. However, these results are a **large file size**, meaning they **cannot be easily synced on GitHub**, and so this walkthrough needs to be **run offline* rather than through GitHub actions.

## How to generate the pre-run results

As the results are not upload to GitHub, if you are trying to re-run these `.qmd` files, you will need to produce your own version of this file.

**TODO: Add instructions on how**

## How to update the walkthrough pages

As these are run offline, if you make any changes to `code0a_walkthrough.qmd`, you will need to go through the following steps to ensure they are updated on the code1-code5 pages in the quarto book.

1. Render `code0a_walkthrough.qmd`, producing `code0a_walkthrough.html.md`
2. Run `code0b_splitmd.qmd`, which will split the results from that markdown file into five seperate files (`code1.md` to `code5.md`)
3. Sync the `code1.md` to `code5.md` files to remote GitHub repository