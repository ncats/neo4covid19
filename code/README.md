* Building the Neo4j database


Workflow tested successfully on Python version 3.6.9, py2neo version 4, and Neo4j version 3.5.6.

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
The process of compiling the database is as follows:


python compile.py

Sometimes UniProt server can be slow, or the connection may fail, so there is a method to use a local maping between genesymbols and UniProt IDs. Use this command to use local gene to UniPRot mapping:

python compile.py test



* Exporting the Neo4j database

(May require stopping the Neo4j db instance, try without stopping it first.)

# Ref: https://neo4j.com/docs/operations-manual/current/tools/dump-load/


neo4j-admin dump --database=graph.db --to=neo4covid-19_neo4j.dump


* Importing the Neo4j database

(May require stopping the Neo4j db instance, try without stopping it first.)


neo4j-admin load --from=neo4covid-19_neo4j.dump --database=graph.db --force




Libraries/components that do not fall under MIT license:

py2neo	Apache License, Version 2.0	[ https://pypi.org/project/py2neo/ ]	[ https://www.apache.org/licenses/LICENSE-2.0 ]

