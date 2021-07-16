Source code and data repository associated with the study disclosed in manuscript:


Zahoranszky-Kohalmi et al., A Workflow of Integrated Resources to Catalyze Network Pharmacology Driven COVID-19 Research. 

Preprint DOI:
DOI:


#LICENSE REMARKS


This repository contains source code, and data and results files which are organized into various subdirectories.

Source code subdirectories:

`code/`
<BR>
`ci/`

Data and results subdirectory:

`data/`


*Source Code License of `neo4covid19` Repository*

The applicable license to source code can be found in the root of the repository under filename: `LICENSE.txt` . This license is applicable to all files recursively in the Source code subdirectories as defined above. The file `NOTES.txt` in the root of the repository lists source code modules that were utilized and their respective licenses. These modules have their own licenses which might be different from the Source Code License of this repository, and they need to be respected accordingly.

*Data License of `neo4covid19` Repository*

The applicable license to data and results can be found under filename: `data/LICENSE.txt` . This license is applicable to all files recursively in the Data and results subdirectory  as defined above. The file `data/NOTES.txt` lists input files and resources utilized to compile the Neo4COVID19 and results files that can be considered as derivative work of those resources. These input files and resources have they own licenses which might be different from the Data License of this repository, and they need to be respected accordingly. In the same file we also list which results files can be considered as derivative works, and we also list the the respective ascendent data source(s).




#Reproducing the Neo4COVID19 workflow and database

*Environment*


> Workflow was tested successfully on Python version 3.6.9, py2neo version 4, and Neo4j version 3.5.6.


> Python (conda, see: https://docs.conda.io/projects/conda/en/latest/) environment needs to have the following libraries installed:

	- py2neo
	
		`pip install py2neo`

	- pandas

		`pip install pandas`

		or (if you have conda environent)
	
		`conda install pandas`

	- xlrd

		`conda install xlrd`

	- requests

		`conda install requests`

	- uuid

		`pip install uuid`


> Git Large File Storage support (https://git-lfs.github.com/)

	- Follow instructions on https://git-lfs.github.com/ to get it installed.


> Clone this repository
	
	- Create a directory under which you want to store the repository, navigate to that directory, then execute this command:

	`git clone https://github.com/ncats/neo4covid19.git`


> Configure Neo4j database

	- Make sure to have a username/password set up for your Neo4j server.


> Create necessary configuration file with Neo4j database access credentials

In the root of the checked out (local) repository create a directory 'cfg' that is at the same level as the directory 'code'.
Then, create a file 'db.cfg' inside the freshly created 'cfg' directory, and populate it with the Neo4j connection parameters, in this format (one parameter per line):

<BR>
<BR>
neo4jhost<BR>
7687<BR>
username<BR>
password<BR>
<BR>
<BR>


You need to replace neo4jhost with the hostname of the server running the Neo4j database server. This can be 'localhost' if working on locally on the server.
The 7687 is the default port for the Neo4j BOLT connection.
Replace 'username' and 'password' according to your Neo4j serve configuration.

Once the db.cfg file is created, you can compile the database. Beware, that the script first deletes ALL existing data from a Neo4j database, so that it can perform a clean build everytime.

Also, make sure you have the Neo4j database server running and that you have the connection configuration file set up properly. Beware, that the workflow will wipe the Neo4j database clean. In case you have different data loaded to your Neo4j database (that have different node/relation types), those may not be wiped by the code. In order to assure a clean Neo4COVID19 deployment, you need to make sure that you start from a completely empty database. Please refer to Neo4j documentaion for details (https://neo4j.com/). Once you took care of the database you can proceed with the workflow as shown below.


*Process of building the Neo4COVID19 database*

> Step 1

Assembling VIHPs.

`python harmonize.py`


Please be patient, this step might take a while, mainly due to the UniProt API.



> Step 2

Assembling a SmartGraph network of HATs.

Go to SmartGraph (https://smartgraph.ncats.io).

Clear the fields "Start Nodes" and "End Nodes" then click on "clear graph".


Copy the content of 'uniprot' column in file `data/output/sg_proteins_a.tsv`. This set of UniProt IDs will be your "Start Nodes" in SmartGraph (https://smartgraph.ncats.io).

Copy the content of 'uniprot' column in file `data/output/sg_proteins_b.tsv` and use them as "End Nodes" in SmartGraph.

Set the "Max Distance" parameter to 3.

Leave the "PPI Confidence Level" at its default value, i.e. 0.00.

Click on "find shortest path".

Once the network is assmebled in SmartGraph, click on "Download graph", select "Cytoscape JSON", then rename the downloaded file to `SG_HATs_dist_3_conf_0.00.json` and place the file into `data/input/`.


Repeat the SmartGraph analysis by swapping the starting and end nodes to create a "reverse" network.

Save the results into `data/input/` as `SG_HATs_reverse_dist_3_conf_0.00.json`.

Process data generated by SmartGraph:

`python process_sg.py`


> Step 3

Compiling the Neo4j database.


`python compile.py`


Please be patient, this step might take a while, mainly due to the UniProt API.


The process to build the Neo4COVID19 database concludes here.



*REMARKS*



Alternatively you can import the precompiled Neo4COVID19 database that is also distributed in this repository (please make sure you have Git LFS support enabled on your system to actually clone the database dump file).



*Importing the Neo4COVID19 Database*

Ref: https://neo4j.com/docs/operations-manual/current/tools/dump-load/

1. Stop the Neo4j database.

2. Import the Neo4COVID19 database from the database dump file.

Example command on a Linux OS:


`neo4j-admin load --from=/home/user/path_to/neo4covid19/data/db/Neo4COVID19_PUB_rolling_neo4j.dump --database=graph.db --force`

Please note that you need to replace the `/home/user/path_to/` part of the file path so that it reflects your environment. Also, here we assumed that the name of the Neo4j database is the default `graph.db`.

<BR>
If you're using a Mac, replace the `/home/user/` of the path by `/Users/user/`.

<BR>
<BR>

3. Start the Neo4j database.




In case you'd make changes in the Neo4COVID19 database, you can export the Neo4j database into a database dump file as shown below.


*Exporting the Neo4j database*


Ref: https://neo4j.com/docs/operations-manual/current/tools/dump-load/

1. Stop the Neo4j database.

2. Dumping the data into a file.

`neo4j-admin dump --database=graph.db --to=/home/user/path_to/neo4covid19/data/db/Neo4COVID19_PUB_rolling_neo4j.dump`

Please note that you need to replace the `/home/user/path_to/` part of the file path so that it reflects your environment. Also, here we assumed that the name of the Neo4j database is the default `graph.db`.

<BR>
If you're using a Mac, replace the `/home/user/` of the path by `/Users/user/`.

<BR>
<BR>



We wish you success in your research to find a cure for COVID-19!
