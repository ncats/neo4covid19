Source code and data repository associated with the study disclosed in manuscript:


Zahoranszky-Kohalmi et al., A Workflow of Integrated Resources to Catalyze Network Pharmacology Driven COVID-19 Research. 

Preprint DOI:
DOI:


#LICENSE REMARKS


This repository contains source code, and data and results files which are organized into various subdirectories.

Source code subdirectories:

code/
ci/

Data and results subdirectory:

data/


*Source Code License of `neo4covid19` Repository*

The applicable license to source code can be found in the root of the repository under filename: `LICENSE.txt` . This license is applicable to all files recursively in the Source code subdirectories as defined above. The file `NOTES.txt` in the root of the repository lists source code modules that were utilized and their respective licenses. These modules have their own licenses which might be different from the Source Code License of this repository, and they need to be respected accordingly.

*Data License of `neo4covid19` Repository*

The applicable license to data and results can be found under filename: `data/LICENSE.txt` . This license is applicable to all files recursively in the Data and results subdirectory  as defined above. The file `data/NOTES.txt` lists input files and resources utilized to compile the Neo4COVID19 and results files that can be considered as derivative work of those resources. These input files and resources have they own licenses which might be different from the Data License of this repository, and they need to be respected accordingly. In the same file we also list which results files can be considered as derivative works, and we also list the the respective ascendent data source(s).




#Reproducing the Neo4COVID19 workflow and database




Process:

*Step 1*

Assembling VIHPs.

`python prepare.py`


If you get API related error, then you can use the 'test' switch, i.e. `python prepare.py test`, which will retrieve mapping from a saved results of a successful API cal in the past when using this workflow.

*Step 2*

Assembling a SmartGraph network of HATs.

Go to SmartGraph (https://smartgraph.ncats.io).

Clear the fields "Start Nodes" and "End Nodes" then click on "clear graph".


Copy the content of 'uniprot_id' column in file `data/input/HATs.tsv`. This set of UniProt IDs will be your "Start Nodes" in SmartGraph (https://smartgraph.ncats.io).

Take the output of Step 1 located at `data/output/unique_host_proteins_prestring.txt`. Copy the UniProt IDs from column 'uniprot_id' and use them as "End Nodes" in SmartGraph.

Set the "Max Distance" parameter to 3.

Leave the "PPI Confidence Level" at its default value, i.e. 0.00.

Click on "find shortest path".

Once the network is assmebled in SmartGraph, click on "Download graph", select "Cytoscape JSON", then rename the downloaded file to `SG_HATs_dist_3_conf_0.00.json` and place the file into `data/input/`.




Repeat the SmartGraph analysis by swapping the starting and end nodes to create a "reverse" network.

Save the results into `data/input/` as `SG_HATs_reverse_dist_3_conf_0.00.json`.


*Step 3*

Compiling the Neo4j database.


Make sure you have the Neo4j database server runing and that you have the connection configuration file set up properly (plese refer to code/README.md). Beware, that this step will wipe the Neo4j database clean. In case you have different data loaded to your Neo4j database (that have different node/relation types), those may not be wiped by the code. In order to assure a clean Neo4COVID19 deployment, you need to make sure that you start from a completely empty database. Please refer to Neo4j documentaion for details (https://neo4j.com/). Once you took care of the database you can proceed with the workflowas shown below.

`python compile.py`

or, in case of API call related error:

`python compile.py test`

(this will use the results of a previous API call saved as a temporary dataset).





The process to build the Neo4COVID19 database concludes here.

We wish you success in your research to find a cure for COVID-19!
