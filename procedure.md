# Ref: https://stackoverflow.com/questions/56713744/how-to-create-conda-environment-with-specific-python-version
# Ref: https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
# Ref: https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-neo4j-on-ubuntu-20-04
# Ref: https://debian.neo4j.com/





Reproducing the workflow

Prerequisites




1. Git Large File Support (Git LFS)

Follow instructions to install Git LFS *before* cloning the Git repository, see: https://git-lfs.github.com/






2. Clone repository

Create a working directory "n4ctest" which is in the user's home directory 
(/home/user/n4ctest in linux, /Users/user/n4ctest on Mac)

Change directory to working directory:

cd ~/n4ctest


Clone repo:


git clone https://github.com/ncats/neo4covid19

cd neo4covid19






3. Setting up Conda environment

- Follow OS specific instructions to install Conda, see: https://docs.conda.io/en/latest/miniconda.html

- Once Conda is installed, create a clean conda environment

conda create -y --name n4c python==3.6.9

- Activate the Conda environment

conda activate n4c

- Install Python dependencies

pip install py2neo
conda install pandas
pip install xlrd==1.2.0
conda install requests
pip install uuid






4. Setting up Neo4j database

Follow instructions on setting up Neo4j for your OS, see:

Example on Ubuntu:

https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-neo4j-on-ubuntu-20-04

sudo apt update

sudo apt install apt-transport-https ca-certificates curl software-properties-common

curl -fsSL https://debian.neo4j.com/neotechnology.gpg.key | sudo apt-key add -

sudo add-apt-repository "deb https://debian.neo4j.com stable 3.5"

sudo apt install neo4j

sudo systemctl enable neo4j.service

sudo systemctl status neo4j.service


(Neo4j database web interface will be running on http://localhost:7474/ and Bolt at 127.0.0.1:7687)


Update default neo4j password to get the credentials required to build the Neo4COVID19 database with the workflow:

The default credentials of any Neo4j installation:

	username: neo4j
	password: neo4j

To change the default password, either visit http://localhost:7474/ where you will be prompted for changing the password, or connect to command-line client from terminal as described below:


cypher-shell

(If not prompted for new password, then enter this into client:


CALL dbms.changePassword('mypassword');



then enter this:

:exit

)


Create Neo4j database connection configuration file:

mkdir -p ~/n4ctest/neo4covid19/cfg

Assuming that the Neo4j is running locally and that you changed the default password to let's say "mypassword", create a "db.cfg" file in ~/n4ctest/neo4covid19/cfg directory with this content and format (one entry per line):


localhost
7687
neo4j
mypassword




5. Replicate workflow


`cd ~/n4ctest/neo4covid19/code/`


# Data harmonization step

`python harmonize.py`

# Perform SmartGraph extension as detailed in Step 2 in the the README.md file of this reposiotry ( https://github.com/ncats/neo4covid19/blob/master/README.md ), then perform next step:

`python process_sg.py`

# DB Compilation

`python compile.py`
