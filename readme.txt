How to run and customise the yeast colony growth model:

O. setting the parameters:
1. open and edit the .csv configuration file
	- each column represents a strain in the "cells parameters" section (max 24 strains) 
2. save the file

I. running the code (basic use):
1. open run_yeast_model.m matlab script
2. set variable "inputname" to the configuration .csv file
3. set variable "outputname" to the name of the folder you would like to save the results (prefix identifying the experiment can be also added)
4. run the script: the function yeast_model.m is called with the parameters of the configuration file
5. click after the initial cell and nutrient layout appears in order to start the simulation (if visualization of the simulation process is enabled)

II. inhomogeneous initial nutrient distribution:
1. open yeast_model.m
2. enable section "inhomogeneous initial nutrient distribution"
3. set the nutrient content of the layers
4. save and run the simulation based on I.

III. nutrient supply during the simulation:
1. open yeast_model.m
2. enable section "nutrient supply at given steps"
3. set the time intervals of nutrient supply, the amount and position of nutrient refill
	- nutrient appears at the edges of the lower layer (p2) by default
4. save and run the simulation based on I.

IV. changing ager viscosity during the simulation (wet-dry):
1. open yeast_model.m
2. enable section "wet-dry change"
3. set the time of change and the new nutrient diffusion and division distance parameters
4. save and run the simulation based on I.

V. putting new agents on the plate during the simulation:
1. open yeast_model.m
2. enable section "put new cells no the plate"
3. set the time of cell insertion, cell type (i) and initial positions for the new cells (position of the drop)
4. save and run the simulation based on I.

VI. save current stae at given time points:
1. open yeast_model.m
2. enable section "save at given time points"
3. set the time points of saving
	- snapshot of the simulation plate is generated and saved
	- regular .mat outputs are saved
4. save and run the simulation based on I.