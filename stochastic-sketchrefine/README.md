﻿# stochastic-sketchrefine

Fill up the database credentials in `Data/database.ini'

And the Gurobi WLS license information in `Utils/GurobiLicense.py'

Ensure that the files 'portfolio.csv' and 'lineitem.tbl' are in the 'Data' folder. The TPC-H dbgen script can be used to build an instance of the lineitem table. Then Setup the database by running:

./dataset_creator.sh

Load all the queries in the directory:

Workload/PortfolioWorkload

To run the algorithms, use:

python main.py
