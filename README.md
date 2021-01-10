# OMSES

This repository contains the OMSES_DE.py and OPlots.py programs, and Microsoft Excel input file templates for the "design envelope" formulation of optimal mine site energy supply design envelope (OMSES_DE). It also contains model00.xlsx, which contains the mine energy demand and wind speed data used in the MASc thesis study of Jeff Eastick.

The python packages listed below are requried:
- pandas
- numpy
- matplotlib
- gurobipy (properietary)

The python program OMSES_DE.py formulates a design envelope OMSES problem using the Gurobi solver Python API. A license for the Gurobi solver can be purchased or requested from Gurobi at www.gurobi.com. 

The python program OPlots.py uses matplotlib to generate interactive or static plots of the results from an OMSES_DE object. 

A forewarning: this program was developed for personal academic use for research works. Software development best practices were not particularly followed.

Acknowledgements:
- Laurentian University and the Bharti School of Engineering
- Mining Innovation Rehabilitation and Applied Research Corporation (MIRARCO)
- Dr. Dean Millar
- Dr. Monica Carvalho
- Dr. Alberto Romero
