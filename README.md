# MetaNetwork
A Weighted Gene Co-expression Network Analysis Rshiny Application

## Installing MetaNetwork for the First Time
MetaNetwork requires the following packages to be loaded in an R session: 
WGCNA, heatmap3, tidyverse, readxl, openxlsx, qtl, corrplot, BiocManager, allez. Further, we recommend using RStudio to run MetaNetwork due to its user-friendly folder navigation interface.

After downloading the full repositroy from GitHub, open RStudio. Go to the file menu and choose "Existing Directory." Click browse and navigate to the folder containing the MetaNetwork repository. Run the following code in the Console to install MetaNetwork's required packages:  

install.packages("WGCNA")
install.packages("heatmap3")
install.packages("tidyverse")
install.packages("readxl")
install.packages("openxlsx")
install.packages("qtl")
install.packages("corrplot")
install.packages("BiocManager")
install.packages("allez")
install.packages("devtools")

Run the following code in the console to launch MetaNetwork: 

shinyAppDir(getwd())

## Step 1: WGCNA Analysis
### Uploading data file
Upload your data to MetaNetwork using the included template excel sheet. Currently, MetaNetwork requires protein identifiers to be Uniprot Accessions. In future updates, this may be expanded. The second column, for gene symbols, is not required. 
### Uploading groups file
Using the included GroupsTemplate.xlsx file, enter the names of your samples in the experiment column. Make sure they match the names of your samples in the DataTemplate.xlsx file. In the SampleID column, enter what experimental condition corresponds to each one of your samples. 

### Uploading Uniprot Database File
The MetaNetwork repository will contain a reviewed, human uniprot database and reviewed, mouse uniprot database. After clicking the upload button, navigate the UniprotDB folder to select the database for your sample. If you require a different species, please download the database from Uniprot. 

### User Input Values
The modifiable values in Step 1: WGCNA Analysis enable custom user control of the clustering, module merging, and minimum module size. They are automatically set to default values commonly used in WGCNA analysis. "Enter how many columns to ignore" is set automatically at two. If you only have a column of acessions with no corresponding gene, set this equal to 1. If your data includes a column of p-values derived from a significance testing indicating significant changes, ensure the column is last, and check the box. 

### Submit Job
After uploading your data, MetaNetwork will display a preview. Make sure everything looks correct, and then click the "Submit Job" button. 

### Output
After the completion of the WGCNA workflow, there will be a new folder entitled "Results" in your working directory. It will contain quality control and diagnostic figures, as well as the module memberships and module eigenproteins. The module memberships can be found in the ResultsWGCNA.xlsx file, while the module eigenproteins can be found in Eigenproteins_by_sample.csv. 

In the future, detailed explanations of each plot will be included in the MetaNetwork repository. Until then, please refer to the section "For Further Reading."

## Step 2: GO Enrichment Analysis
Gene Ontology Enrichment Analysis is performed by the allez package (insert citation). This step only requires uploading the ResultsWGCNA.xlsx file from Results folder and selecting the correct organism from the dropdown selection box. This version of MetaNetwork supports only mouse and human databases.

## Step 3: Phenotype Correlation
This step is not yet fully functional. 

# For further reading

# Citations
