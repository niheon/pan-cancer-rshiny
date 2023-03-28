# pan_cancer_shiny

## R Shiny dashboard app for visualizing Pan Cancer data from Lodi, et al. 2022

This Shiny app accepts as input a list of genes interest and calculates a signature score for each cell for the selected cancer type. Additionally, the app
creates a series of dynamically generated customizable diagnostic plots using Plotly including: bar plot of cell frequencies for select cancers, and a scatter plot correlating cell frequencies with signature scores. The resulting figures can be downloaded as .png or .pdf files. The data is also displayed in a dynamic HTML table that allows specific entries in the dataset to be queried.

The genes of interest are input but users by uploading a CSV file in the following format:

- First column must contain a unique list of genes
- A column header is not strictly necessary, but we suggest using the header "genes"
- The CSV file must contain no row names
- An example CSV input file can be found in the /data directory.

This app is released under GPL-3.0 License.
