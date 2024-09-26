# Shiny application

This folder generates a very basic shiny application which allows the user to select a few inputs and run some basic code related to the model.

It is hoped that this work can be further developed into a comprehensive user interface for people to interact with the model, instead of via the excel input sheet, but funding for this has not been confirmed currently.

## View the app online

This app is hosted with `shinyapps.io` and can be viewed at <https://amyheather.shinyapps.io/shinyapp/>. With [ShinyApps.io set up via rsconnect](https://shiny.posit.co/r/articles/share/shinyapps/), the app was deployed (and is updated) via `deployApp()`:

```
library(rsconnect)
deployApp()
```

## View the app locally

If you want to view this app on your local machine, you'll need to:

1. Set up R environment, using `renv.lock` from parent folder.

```
install.packages("renv")
renv::restore()
```

2. Run the app in this folder, either via console:

```
runApp()
```

Or, if in RStudio, open `app.R` and then select the "Run App" button.