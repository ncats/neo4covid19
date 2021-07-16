# Neo4COVID19 Source Code and Data Repository


This repository is associated with the study disclosed in manuscript:
<BR>
<BR>
*Zahoranszky-Kohalmi et al.*, A Workflow of Integrated Resources to Catalyze Network Pharmacology Driven COVID-19 Research. 

Preprint https://www.biorxiv.org/content/10.1101/2020.11.04.369041v1. (DOI: 10.1101/2020.11.04.369041)
<BR>
<BR>
- A functional instance of the **Neo4COVID19** database is accessible here: https://aspire.covid19.ncats.io:7473/browser/
	
- For detailed information on how to access the database please refer to the instructions posted at https://neo4covid19.ncats.io/ under **ACCESS** tab. 


## Instructions to Reprodue the Workflow Used to Build the Neo4COVID19 Database
<BR>


### Prerequisites

The workflow was tested successfully on Python version 3.6.9, py2neo version 4, and Neo4j version 3.5.6.

1. Install **Git Large File Support (Git LFS)**


	- Follow instructions to install Git LFS *before* cloning the Git repository, see: https://git-lfs.github.com/




2. Clone repository

	- Create a working directory "**n4ctest**" which is in the user's home directory 
	(*/home/user/n4ctest* in linux, */Users/user/n4ctest on Mac*)

	- Change directory to working directory:

		`cd ~/n4ctest`


	- Clone repo:


		`git clone https://github.com/ncats/neo4covid19`

	- Change to source code directory:

		`cd neo4covid19/code`




3. Setting up Conda environment

	- Follow OS specific instructions to install Conda, see: https://docs.conda.io/en/latest/miniconda.html

	- Once Conda is installed, create a clean conda environment

		`conda create -y --name n4c python==3.6.9`

	- Activate the Conda environment

		`conda activate n4c`

	- Install Python dependencies

		`pip install py2neo`
		<BR>
	
		`conda install pandas`
		<BR>
		
		`pip install xlrd==1.2.0`
		<BR>
		
		`conda install requests`
		<BR>
		
		`pip install uuid`
		<BR>



4. Setting up Neo4j database

	- Follow instructions on setting up Neo4j for your OS, see: https://neo4j.com/docs/operations-manual/current/installation/

	- Example on Ubuntu (see: https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-neo4j-on-ubuntu-20-04)

		`sudo apt update`
		<BR>

		`sudo apt install apt-transport-https ca-certificates curl software-properties-common`
		<BR>

		`curl -fsSL https://debian.neo4j.com/neotechnology.gpg.key | sudo apt-key add -`
		<BR>

		`sudo add-apt-repository "deb https://debian.neo4j.com stable 3.5"`
		<BR>

		`sudo apt install neo4j`
		<BR>

		`sudo systemctl enable neo4j.service`
		<BR>

		`sudo systemctl status neo4j.service`


	- (Neo4j database web interface will be running on http://localhost:7474/ and Bolt at 127.0.0.1:7687)


	- Update default neo4j password to get the credentials required to build the Neo4COVID19 database with the workflow:

		The default credentials of any Neo4j installation:
		> **username: neo4j**<BR>**password: neo4j**

	- To change the default password, either visit http://localhost:7474/ where you will be prompted for changing the password, or connect to command-line client from terminal as described below:


		`cypher-shell`

		If not prompted for new password, then enter this into client:


		`CALL dbms.changePassword('mypassword');`

		then enter this:

		`:exit`

		


	- Create Neo4j database connection configuration file:

		`mkdir -p ~/n4ctest/neo4covid19/cfg`

		Assuming that the Neo4j is running locally and that you changed the default password to let's say **mypassword**, create a **db.cfg** file in *~/n4ctest/neo4covid19/cfg* directory with this content and format (one entry per line):


		>localhost<BR>7687<BR>neo4j<BR>mypassword


		If your Neo4j server is not running on *localhost* ten you need to replace *localhost* with the hostname of the server running the Neo4j database server. The *7687* is the default port for the Neo4j BOLT connection. **Replace** *username* and *password* according to your Neo4j serve configuration.

		Once the **db.cfg** file is created, you can compile the database. Beware, that the script first deletes ALL existing data from a Neo4j database, so that it can perform a clean build everytime.

		Also, make sure you have the Neo4j database server running and that you have the connection configuration file set up properly. The workflow will wipe the Neo4j database clean. In case you have different data loaded to your Neo4j database (that have different node/relation types), those may not be wiped by the code. In order to assure a clean Neo4COVID19 deployment, you need to make sure that you start from a completely empty database. Please refer to Neo4j documentaion for details (https://neo4j.com/). Once you took care of the database you can proceed with the workflow as shown below.

### Replicate workflow


- Change to source code directory

	`cd ~/n4ctest/neo4covid19/code/`

>Step 1. Data harmonization step

	`python harmonize.py`

	

>Step 2. Assembling a SmartGraph Network of HATs

**If you just want to rebuild the DB skip to Step 3**, otherwise continue as follows.
			
- Go to SmartGraph (https://smartgraph.ncats.io).

- Clear the fields "Start Nodes" and "End Nodes" then click on "clear graph".


- Copy the content of **uniprot** column in file `data/output/sg_proteins_a.tsv`. This set of UniProt IDs will be your *Start Nodes* in SmartGraph (https://smartgraph.ncats.io).

- Copy the content of *uniprot* column in file `data/output/sg_proteins_b.tsv` and use them as *End Nodes* in SmartGraph.

- Set the *Max Distance* parameter to 3.

- Leave the *PPI Confidence Level* at its default value, i.e. 0.00.

- Click on "find shortest path".

- Once the network is assmebled in SmartGraph, click on "Download graph", select "Cytoscape JSON", then rename the downloaded file to `SG_HATs_dist_3_conf_0.00.json` and place the file into `data/input/`.


- Repeat the SmartGraph analysis by swapping the starting and end nodes to create a "reverse" network.

- Save the results into `data/input/` as `SG_HATs_reverse_dist_3_conf_0.00.json`.

>Step 3. Process Data Generated by SmartGraph (in Step 2) 
			
	`python process_sg.py`
	
			
>Step 4. DB Compilation

	`python compile.py`


<BR>
<BR>
This concludes the replication procedure, the resultant Neo4COVID19 database should be availabel at http://localhost:7474 in your browser, if you used the localhost as your database URL during the replication process.



### REMARKS


Alternatively you can import the precompiled Neo4COVID19 database that is also distributed in this repository (please make sure you have Git LFS support enabled on your system to actually clone the database dump file).


- Importing the Neo4COVID19 Database


1. Stop the Neo4j database.

2. Import the Neo4COVID19 database from the database dump file.

	Example command on a Linux OS:


	`neo4j-admin load --from=/home/user/path_to/neo4covid19/data/db/Neo4COVID19_PUB_rolling_neo4j.dump --database=graph.db --force`

	Please note that you need to replace the */home/user/path_to/* part of the file path so that it reflects your environment. Also, here we assumed that the name of the Neo4j database is the default `graph.db`.

	<BR>
	If you're using a Mac, replace the */home/user/* of the path by */Users/user/*.

<BR>
<BR>

3. Start the Neo4j database.


In case you'd make changes in the Neo4COVID19 database, you can export the Neo4j database into a database dump file as shown below.


- Exporting the Neo4j database


1. Stop the Neo4j database.

2. Dumping the data into a file.

	`neo4j-admin dump --database=graph.db --to=/home/user/path_to/neo4covid19/data/db/Neo4COVID19_PUB_rolling_neo4j.dump`

	Please note that you need to replace the */home/user/path_to/* part of the file path so that it reflects your environment. Also, here we assumed that the name of the Neo4j database is the default `graph.db`.

	<BR>
	If you're using a Mac, replace the */home/user/* of the path by */Users/user/*.

<BR>
<BR>


### LICENSE REMARKS


This repository contains source code, and data and results files which are organized into various subdirectories.

Source code subdirectories:

`code/`
<BR>
`ci/`

Data and results subdirectory:

`data/`


- Source Code License of **neo4covid19** Repository

	The applicable license to source code can be found in the root of the repository under filename: `LICENSE.txt` . This license is applicable to all files recursively in the Source code subdirectories as defined above. The file `NOTES.txt` in the root of the repository lists source code modules that were utilized and their respective licenses. These modules have their own licenses which might be different from the Source Code License of this repository, and they need to be respected accordingly.

- Data License of **neo4covid19** Repository

	The applicable license to data and results can be found under filename: `data/LICENSE.txt` . This license is applicable to all files recursively in the Data and results subdirectory  as defined above. The file `data/NOTES.txt` lists input files and resources utilized to compile the Neo4COVID19 and results files that can be considered as derivative work of those resources. These input files and resources have they own licenses which might be different from the Data License of this repository, and they need to be respected accordingly. In the same file we also list which results files can be considered as derivative works, and we also list the the respective ascendent data source(s).






We wish you success in your research to find a cure for COVID-19!
	
	
	
### References

https://stackoverflow.com/questions/56713744/how-to-create-conda-environment-with-specific-python-version
https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-neo4j-on-ubuntu-20-04
https://debian.neo4j.com/
https://www.markdownguide.org/cheat-sheet/
https://neo4j.com/docs/operations-manual/current/tools/dump-load/



