#%%
# -*- coding: utf-8 -*-
"""
    Optimal Mine Site Energy Supply, main.py (sample)
    Copyright (C) 2020, Jeff Eastick

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact:
    jeffeastick@gmail.com
    jx_eastick@laurentian.ca
"""
#%%



#%% Run Optimization Models:
# from OMSES_LT import *
from OMSES_DE import *
from OPlots import *  

# Create and run model00. Change 
# User must change the raw string path to the correct input file location.
m00 = OMSES_DE(r'C:\Users\jeffe\Documents\GitHub\OMSES\input\model00.xlsx')
m00.run_all() 

#use OPlots module to graph various results:
# plot_structure_at_i(m00,7)
# plt_util_flow(m00,'EE')

# delete the model and garbage collect after each model if you are going to run several models as a batch, otherwise you may run out of memory.
# del m00
# gc.collect()





