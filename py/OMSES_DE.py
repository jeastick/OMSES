#%%
# -*- coding: utf-8 -*-
"""
    Optimal Mine Site Energy Supply, Design Envelope (OMSES_DE)
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

from gurobipy import *
import pandas as pd
import numpy as np
from OPlots import *
#%%
class OMSES_DE:
  

  #%% This function returns the total number of representative hours in a period 'j'
  def dt_j(self,j):
    return 24


  #%% These are economic discrete compounding factor formulae

  def econ_PAIN(self,i,N):
    if i == 0:
      return N
    else:
      return ((1+i)**N-1)/(i*(1+i)**N)
  
  def econ_PFIN(self,i,N):
    return 1/(1+i)**N
  
  #%% This function takes a (correctly Excel-formatted) dataframe of timeseries data and parses it to a Gurobi tupledict of form dict[(key,i,j)] = value
  def timeseries_df_to_tupledict(self,df):
    temp_dict = df.to_dict(orient='list')
    ti = temp_dict.pop('i')
    tj = temp_dict.pop('j')
    tupdic = tupledict()
    for key in temp_dict:
      for step in range(self.ij):
        tupdic[tuple((key,ti[step],tj[step]))] = temp_dict[key][step]
    return tupdic

  #%% This function returns the completed (solved) superstructure results for period i,j from a solved model
  def get_supr_at(self,i,j):
    util_rows = len(self.u)
    tech_cols = len(self.v)
    wind_cols = len(self.w)
    stor_cols = len(self.x)
    
    imp_mat = np.zeros((util_rows,1))
    tec_mat = np.zeros((util_rows,tech_cols))
    wnd_mat = np.zeros((util_rows,wind_cols))
    sto_mat = np.zeros((util_rows,stor_cols))
    
    exp_mat = np.zeros((util_rows,1))
    was_mat = np.zeros((util_rows,1))
    dem_mat = np.zeros((util_rows,1))
    
    for util in range(util_rows):
      imp_mat[util][0] =  self.dv_util_flow_imp[(self.u[util],i,j)].getAttr('X')
      exp_mat[util][0] = -self.dv_util_flow_exp[(self.u[util],i,j)].getAttr('X')
      was_mat[util][0] = -self.dv_util_flow_was[(self.u[util],i,j)].getAttr('X')
      dem_mat[util][0] = -self.dv_util_flow_dem[(self.u[util],i,j)].getAttr('X')

    for util in range(util_rows):
      for tech in range(tech_cols):
        tec_mat[util][tech] = self.dv_util_tech[(self.u[util],self.v[tech],i,j)].getAttr('X')
    
    for util in range(util_rows):
      for wind in range(wind_cols):
        wnd_mat[util][wind] = self.dv_util_wind[(self.u[util],self.w[wind],i,j)].getAttr('X')
        
    for util in range(util_rows):
      for stor in range(stor_cols):
        sto_mat[util][stor] = self.dv_util_stor[(self.u[util],self.x[stor],i,j)].getAttr('X')
    
    ss_mat = np.concatenate((imp_mat, tec_mat, wnd_mat, sto_mat, exp_mat, was_mat, dem_mat), axis = 1)   
      
    return ss_mat
  
  #%% This function returns all (solved) superstructure matrices from a solved model, in the form of a dictionary with tuple keys (i,j)
  def gen_supr_at_all(self):
    self.d_suprs = {}
    for i in self.i:
      for j in self.j:
        self.d_suprs[(i,j)] = self.get_supr_at(i,j)
        
  #%% This function returns the superstructure sum for period i
  
  def get_supr_at_i(self,i):
    ssi = np.zeros(self.d_suprs[(i,1)].shape)
    
    for j in self.j:
      ssi = np.add(ssi,self.d_suprs[i,j])
      
    return ssi
      
  def gen_supr_at_i(self):
    
    self.d_suprs_i = {}
    
    for i in self.i:
      self.d_suprs_i[i] = np.zeros(self.d_suprs[(i,1)].shape)
      for j in self.j:
        self.d_suprs_i[i] = np.add(self.d_suprs_i[i],self.d_suprs[i,j])
        

  #%% This function returns the maximum value of all dv_util_flow_xxx variables
  def get_max_dv_util_flow(self):
    m = 0
    for i in self.i:
      for j in self.j:
        abs_max = np.abs(np.amax(self.d_suprs[i,j]))
        abs_min = np.abs(np.amin(self.d_suprs[i,j]))
        lm = max(abs_max,abs_min)
        if lm > m:
          m = lm
    return m
        
  def get_max_dv_util_flow_i(self):
    m = 0
    for i in self.i:
      abs_max = np.abs(np.amax(self.d_suprs_i[i]))
      abs_min = np.abs(np.amin(self.d_suprs_i[i]))
      lm = max(abs_max,abs_min)
      if lm > m:
        m = lm
    return m        
  
  def get_max_chrg_dsch_flow(self):
    m = 0
    
    for x in self.x:
      array = self.get_storage_ch_array(x)
      abs_max = np.abs(np.amax(array[1])) #charge value
      abs_min = np.abs(np.amin(array[2])) #discharge value
      lm = max(abs_max,abs_min)
      if lm > m:
        m = lm
    return m  

#%%

  def util_to_excel(self,u):
    data = self.get_util_array(u)
    df = pd.DataFrame(data).T
    labels = self.v.copy() + self.w.copy() + self.x.copy()
    labels.insert(0,"Imports")
    labels.append("Exports")
    labels.append("Waste")
    labels.append("Demand")
    df.columns = labels
    filename = "_utilflow_" + self.modelname + "_" + u + ".xlsx"
    df.to_excel(filename, index = False)
    
    
  def bigM_is_valid(self):
    return ((self.bigM > self.get_max_chrg_dsch_flow()))
        
  #%% This function returns the values of all linear expressions used in calculating the Total Cost
  
  def get_costs(self):
    costs = {}
    costs["c_fix_techCAPEX"] = self.c_fix_techCAPEX.getValue()
    costs["c_fix_windCAPEX"]=self.c_fix_windCAPEX.getValue()
    costs["c_fix_storCAPEX"]=self.c_fix_storCAPEX.getValue()
    costs["c_fix_utilCAPEX_ext"] = self.c_fix_utilCAPEX_ext.getValue()
    costs["c_fix_utilCAPEX_int"] = self.c_fix_utilCAPEX_int.getValue()
    costs["c_fix_techOPEX_fix"] = self.c_fix_techOPEX_fix.getValue()
    costs["c_fix_windOPEX_fix"]=self.c_fix_windOPEX_fix.getValue()
    costs["c_fix"] = self.c_fix.getValue()
    
    costs["c_var_imports"] = self.c_var_imports.getValue()
    costs["c_var_exports"] = self.c_var_exports.getValue()
    costs["c_var_techOPEX_var"] = self.c_var_techOPEX_var.getValue()
    costs["c_var_windOPEX_var"] = self.c_var_windOPEX_var.getValue()
    costs["c_var"] = self.c_var.getValue()
    
    costs["c_tot"] = self.c_tot.getValue()
    return costs
  
  #%% This function returns the number of technologies selected dv_num_techs
  
  def get_dv_num_techs(self):
    return self.dv_num_techs
  
  
  #%% This function returns the number of wind technologies selected dv_num_winds
  def get_dv_num_winds(self):
    return self.dv_num_winds
  
    #%% This function returns the max storage capacity of each storage technology 
  def get_dv_max_store(self):
    return self.dv_max_store
  
      #%% This function returns the max storage capacity of each storage technology 
  def get_dv_min_store(self):
    return self.dv_min_store
  
  #%% This function returns the utility capacity v_util_cpct
  
  def get_dv_util_cpct_int(self):
    return self.dv_util_cpct_int
  
  def get_dv_util_cpct_ext(self):
    return self.dv_util_cpct_ext  
   
    
  
  #%% Class init function
  def __init__(self,datafile):
  
    self.datafile = datafile
    
    
    print("-----------------------------------------------------------------------")
    print("----------------------------OMSES_DE-----------------------------------")
    print("-----------------------------------------------------------------------")
    print("File: " + self.datafile + "\n")
    print("-----------------------------------------------------------------------")
    print('')
    
    ### Load General Data and Set-up for future data imports
    self.file_in_general = pd.read_excel(self.datafile,
                                    sheet_name="General"
                                    ).set_index("VARIABLE").fillna(0)
    print("General data imported;"+"\n")
    
    #Time Scale Temporal Division
    self.d_gen  = self.file_in_general.T.to_dict()
    
    
    #Note: this file name string manipulation to get the filename will only work with Windows systems:
    self.modelname = self.datafile.split('\\')[len(self.datafile.split('\\'))-1].strip('.xlsx')
        
    print("Model name is " + self.modelname)
    
    #for i 
    self.ni = int(self.d_gen['i']['VALUE'])
    if self.ni == 1:
      self.i = self.ni
    else:
      self.i = range(1,self.ni+1)
    
    #for j
    self.nj = int(self.d_gen['j']['VALUE'])
    if self.nj == 1:
      self.j = self.nj
    else:
      self.j = range(1,self.nj+1)
    
    self.ij = self.ni*self.nj
    
    self.discount = self.d_gen['discount']['VALUE']
    
    self.bigM = self.d_gen['bigM']['VALUE']

    self.econ_PAIN_at = {}
    
    self.econ_PFIN_at = {}
    
    
    for i in self.i:
      self.econ_PAIN_at[i] = self.econ_PAIN(self.discount,i)
      self.econ_PFIN_at[i] = self.econ_PFIN(self.discount,i)
      
    # These self.months and self.seasons aren't used but may be used one day.

    self.months = {
        1 : 31, #JAN
        2 : 28, #FEB
        3 : 31, #MAR
        4 : 30, #APR
        5 : 31, #MAY
        6 : 30, #JUN
        7 : 31, #JUL
        8 : 31, #AUG
        9 : 30, #SEP
        10: 31, #OCT
        11: 30, #NOV
        12: 31, #DEC
        }
    
    self.seasons = {
        1 : 90, # WINTER, JAN/FEB/MAR
        2 : 91, # SPRING, APR/MAY/JUN
        3 : 92, # SUMMER, JUL/AUG/SEP
        4 : 92, # FALL,   OCT/NOV/DEC
        }
    
    print("General data setup complete.\n")
    
  #%%Read Excel file and Parse to Dictionaries

  def import_Data(self):
 
    #-----------------
    #----UTILITIES----
    #-----------------
        
    self.file_in_utilities = pd.read_excel(self.datafile,
                                      sheet_name="Utilities",
                                      header=1
                                      ).set_index("UTILITY").fillna(0)
    print("Utility data imported;")
    self.d_util = self.file_in_utilities.T.to_dict()
    print("Utility data setup complete.\n")
    
    self.file_in_import_prices = pd.read_excel(self.datafile,
                                          sheet_name="ImportPrices",
                                          header=1,
                                          nrows=self.ij,
                                          usecols=(list(range(len(self.d_util)+3)))
                                          ).set_index('step').fillna(0)
    print("ImportPrices data imported;")
    self.d_import_prices = self.timeseries_df_to_tupledict(self.file_in_import_prices)  
    print("ImportPrices data setup complete.\n")
    
    
    self.file_in_export_prices = pd.read_excel(self.datafile,
                                          sheet_name="ExportPrices",
                                          header=1,
                                          nrows=self.ij,
                                          usecols=(list(range(len(self.d_util)+3)))
                                          ).set_index('step').fillna(0)
    print("ExportPrices data imported;")
    self.d_export_prices = self.timeseries_df_to_tupledict(self.file_in_export_prices)  
    print("ExportPrices data setup complete.\n")    
    
    #------------
    #----TECH----
    #------------
    
    self.file_in_technologies = pd.read_excel(self.datafile,
                                         sheet_name="Technologies",
                                         header=1
                                         ).set_index("TECHNOLOGY").fillna(0)
    print("Technology data imported;")
    self.d_tech = self.file_in_technologies.T.to_dict()
    print("Technology data setup complete.\n")
    
    self.file_in_superstructure = pd.read_excel(self.datafile,
                                           sheet_name="Superstructure",
                                           nrows=len(self.d_util),
                                           usecols=list(range(len(self.d_tech)+1)),
                                           ).set_index("SUPERSTRUCTURE").fillna(0)
    print("Superstructure data imported;")
    self.mat_supr = self.file_in_superstructure.to_numpy()[:len(self.d_util),:len(self.d_tech)]
    self.mat_supr_pro = np.sign(self.mat_supr).clip(min=0)
    self.mat_supr_con = np.sign(self.mat_supr).clip(max=0)
    print("Superstructure data setup complete.\n")   
    
    #------------
    #----WIND----
    #------------
    
    self.file_in_wind_technologies = pd.read_excel(self.datafile,
                                         sheet_name="WindTechnologies",
                                         header=1
                                         ).set_index("WINDTECHNOLOGY").fillna(0)
    print("WindTechnology data imported;")
    self.d_windtech = self.file_in_wind_technologies.T.to_dict()
    print("WindTechnology data setup complete.\n")
    
    self.file_in_windsuperstructure = pd.read_excel(self.datafile,
                                           sheet_name="WindSuperstructure",
                                           nrows=len(self.d_util),
                                           usecols=list(range(len(self.d_windtech)+1)),
                                           ).set_index("SUPERSTRUCTURE").fillna(0)
    print("Wind Superstructure data imported;")
    self.mat_wind_supr = self.file_in_windsuperstructure.to_numpy()[:len(self.d_util),:len(self.d_windtech)]
    self.mat_wind_supr_pro = np.sign(self.mat_wind_supr).clip(min=0)
    self.mat_wind_supr_con = np.sign(self.mat_wind_supr).clip(max=0)
    print("Wind Superstructure data setup complete.\n")        
    

    self.file_in_wind_data = pd.read_excel(self.datafile,
                                          sheet_name="Wind",
                                          header=1, 
                                          nrows=self.ij,
                                          usecols=(list(range(1+3)))
                                          ).set_index('step').fillna(0)
    print("Wind speed data imported;")
    self.d_ws = self.timeseries_df_to_tupledict(self.file_in_wind_data)  
    print("Wind speed data setup complete.\n")
    
    #---------------
    #----STORAGE----
    #---------------
    
    self.file_in_storage = pd.read_excel(self.datafile,
                                         sheet_name="Storage",
                                         header=1
                                         ).set_index("STORAGETECHNOLOGY").fillna(0)
    print("Storage Technology data imported;")
    self.d_stor = self.file_in_storage.T.to_dict()
    print("Storage data setup complete.\n")
    
    self.file_in_storagesuperstructure = pd.read_excel(self.datafile,
                                           sheet_name="StorageSuperstructure",
                                           nrows=len(self.d_util),
                                           usecols=list(range(len(self.d_stor)+1)),
                                           ).set_index("SUPERSTRUCTURE").fillna(0)
    print("Storage Superstructure data imported;")
    self.mat_stor_supr = self.file_in_storagesuperstructure.to_numpy()[:len(self.d_util),:len(self.d_stor)]
    self.mat_stor_supr_pro = np.sign(self.mat_supr).clip(min=0)
    self.mat_stor_supr_con = np.sign(self.mat_supr).clip(max=0)
    print("Storage Superstructure data setup complete.\n")    
    
    
    #---------------
    #----DEMANDS----
    #---------------
    
    self.file_in_demands = pd.read_excel(self.datafile,
                                          sheet_name="Demands",
                                          header=1,
                                          nrows=self.ij,
                                          usecols=(list(range(len(self.d_util)+3)))
                                          ).set_index('step').fillna(0)
    print("Demand data imported;")
    self.d_demands = self.timeseries_df_to_tupledict(self.file_in_demands)  
    print("Demand data setup complete.")
    
    self.file_in_peaks = pd.read_excel(self.datafile,
                                          sheet_name="Peaks",
                                          header=1,
                                          nrows=self.ij,
                                          usecols=(list(range(len(self.d_util)+3)))
                                          ).set_index('step').fillna(0)
    print("Peak data imported;")
    self.d_peaks = self.timeseries_df_to_tupledict(self.file_in_peaks)  
    print("Peak data setup complete.")
    
    
        
    
    
    print("All data importing is complete.\n")
  
    #%% Define Sets and Input Data
    print("Defining sets u, v, w, x...")
    

    ###### To access a conversion factor,   using iteration, use   mat_supr[row][col]
    ######                                  using string names     mat_supr[u.index("utility_name")[v.index("tech_name")]]
    
    
    #Utilities (u)
    self.u = list(self.d_util.keys())
    print("u is:")
    print(self.u)
    
    #Technologies (v)
    self.v = list(self.d_tech.keys())
    print("v is:")
    print(self.v)
    print()
    #Superstructure Matrix
    print("Tech Superstructure Matrix:")
    print(self.mat_supr)
    print()
    
    #WindTechnologies (w)
    self.w = list(self.d_windtech.keys())
    print("w is:")
    print(self.w)
    print()
    print("Wind Superstructure Matrix:")
    print(self.mat_wind_supr)
    print()  
    
    # Wind speed x-coordinate for wind power coefficient interpolation function. Y-values will come from the imported d_windtech data
    self.wdr = 25 # wind speed data range for turbine performance coefficients, 0 m/s to 25 m/s
    self.xws = np.linspace(0,self.wdr,self.wdr+1)
    
    for wt in self.d_windtech:
      cp_array = np.zeros(self.wdr+1)
      for data in range(self.wdr+1):
        access_key = 'p' + str(data)
        cp_array[data] = self.d_windtech[wt][access_key]
      self.d_windtech[wt]['yp'] = cp_array
      
      
    #Wind power penetration limit
    self.wind_pen_limit = 1
    
    #StorageTechnologies (x)
    
    self.x = list(self.d_stor.keys())
    print("x is:")
    print(self.x)
    print()
    print("Storage Superstructure Matrix:")
    print(self.mat_stor_supr)
    print()      
    
    ###NOTES###
    ###### To access a technology/utility dictionary item through a numerical iteration you use d_tech[v[index]]["Key"]]
    ###### To access a technology/utility dictionary item with text use d_tech["technology_name"]["key"]
    

    print("Sets u, v, and w have been defined.\n")
    
    
  #%% Define the model and decicion variables
  def define_DVs(self):
    
    # Define Model
    print("Defining Model...")
    self.model = Model('OMSES')
    print("Model has been defined.\n")
    self.logfilename = "logfile_" + self.modelname + ".log"
    self.model.setParam("LogFile", self.logfilename)
    self.model.setParam("MIPGap", 0)    # default is 1e-4, set to 0 so we find the absolute optimal solution.
    self.model.setParam("MIPGapAbs", 0) # default is 1e-10, set to 0 so we find the absolute optimal solution.
    # Define Variables
    print("Defining variables...")
    
    # Variables prefixed with dv = variables/decision variables for the purposes of optimization.
    # Variables prefixed with pk = variables/decision variables for the purposes of checking peak demands against supply system capacity
    #------------
    #----TECH----
    #------------
    
    # Nominal:
    self.dv_num_techs = self.model.addVars(self.v,                                        vtype=GRB.INTEGER,      name="num_techs") # mu_v - Number of units of each technology
    self.dv_prod_tech = self.model.addVars(self.v,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="prod_tech") # PI_vij - Production by technology (MWh)    
    self.dv_kprod_tech= self.model.addVars(self.v,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="kprod_tech") # PI_vij - Production by technology (TWh)
    self.dv_util_tech = self.model.addVars(self.u,self.v,self.i,self.j, lb=-1e21,         vtype=GRB.CONTINUOUS,   name="util_tech") # X_uvij - Utility flows by technology
    
    # Peak:
    self.pk_prod_tech = self.model.addVars(self.v,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="pk_prod_tech") # PI^p vij - Peak Capacity by technology (MW)  
    self.pk_util_tech = self.model.addVars(self.u,self.v,self.i,self.j, lb=-1e21,         vtype=GRB.CONTINUOUS,   name="pk_util_tech") # X^p uvij - Utility flows by technology (MW)

    #------------
    #----WIND----
    #------------
     
    # Nominal
    self.dv_num_winds = self.model.addVars(self.w,                                        vtype=GRB.INTEGER,      name="num_winds") # mu_w - Number of units of each wind technology    
    self.dv_prod_wind = self.model.addVars(self.w,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="prod_wind") # PI_wij - Production by wind technology   
    self.dv_util_wind = self.model.addVars(self.u,self.w,self.i,self.j,                   vtype=GRB.CONTINUOUS,   name="util_wind") # X_uwij - Utility flows by wind technology

    # Peak
    self.pk_prod_wind = self.model.addVars(self.w,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="pk_prod_wind") # PI^p_wij - Production by wind technology (MW)   
    self.pk_util_wind = self.model.addVars(self.u,self.w,self.i,self.j,                   vtype=GRB.CONTINUOUS,   name="pk_util_wind") # X^p_uwij - Utility flows by wind technology (MW)  

    #for penetration factor limiting:
    if self.wind_pen_limit < 1:
      self.dv_util_tech_pro = self.model.addVars(self.u,self.v,self.i,self.j,               vtype=GRB.CONTINUOUS,   name="util_tech_pro") # X_uvij - Utility flows by technology
      self.dv_util_wind_pro = self.model.addVars(self.u,self.w,self.i,self.j,               vtype=GRB.CONTINUOUS,   name="util_wind_pro") # X_uwij - Utility flows by wind technology
      self.dv_util_prod_tec = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_prod_tec") # Utility production, just wind technologies
      self.dv_util_prod_wnd = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_prod_wnd") # Utility production, just wind technologies

    
    #---------------
    #----STORAGE----
    #---------------
    
    self.dv_max_store = self.model.addVars(self.x,                                        vtype=GRB.CONTINUOUS,   name="max_store")  # Psi^max_x storage capacity maximum 
    self.dv_min_store = self.model.addVars(self.x,                                        vtype=GRB.CONTINUOUS,   name="min_store")  # Psi^min_x storage capacity minimum    
    self.dv_lvl_store = self.model.addVars(self.x,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="lvl_store")  # Psi_xij storage level 
    
    self.dv_dir_store = self.model.addVars(self.x,self.i,self.j,                          vtype=GRB.BINARY,       name="dir_store")  # delta_xij charge flow direction (0 -> discharge, 1 -> charge)
  
    self.dv_util_chrg = self.model.addVars(self.x,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="util_chrg")  # Phi^ch_xij - Utility charge by storage tech
    self.dv_util_dsch = self.model.addVars(self.x,self.i,self.j,                          vtype=GRB.CONTINUOUS,   name="util_dsch")  # Phi^ds_xij - Utility discharge by storage tech
    
    self.dv_util_stor = self.model.addVars(self.u,self.x,self.i,self.j,    lb=-1e21,      vtype=GRB.CONTINUOUS,   name="util_stor")  # X_uxij - Utility flows by storage tech (= charge - discharge)


    # Peak - Set this to zero in constraints
    self.pk_util_stor = self.model.addVars(self.u,self.x,self.i,self.j,    lb=-1e21,      vtype=GRB.CONTINUOUS,   name="pk_util_stor")  # X^p_uwij - Peak utility flows by storage tech (= charge - discharge)


    #--------------------------
    #----NET ENERGY BALANCE----
    #--------------------------
    
    # Nominal
    self.dv_util_flow_imp = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_flow_imp") # Phi_uij - Utility flow - imported (MWh)
    self.dv_kutil_flow_imp = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="kutil_flow_imp") # Phi_uij - Utility flow - imported (GWh)
    self.dv_util_flow_exp = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_flow_exp") # Phi_uij - Utility flow - exported
    self.dv_util_flow_tec = self.model.addVars(self.u,self.i,self.j,    lb=-1e21,         vtype=GRB.CONTINUOUS,   name="util_flow_tec") # Phi_uij - Utility flow - technologies
    self.dv_util_flow_wnd = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_flow_wnd") # Phi_uij - Utility flow - wind technologies
    self.dv_util_flow_sto = self.model.addVars(self.u,self.i,self.j,    lb=-1e21,         vtype=GRB.CONTINUOUS,   name="util_flow_sto") # Phi_uij - Utility flow - storage
    self.dv_util_flow_was = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_flow_was") # Phi_uij - Utility flow - waste
    self.dv_util_flow_dem = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="util_flow_dem") # Phi_uij - Utility flow - demand


    # Peak
    self.pk_util_flow_imp = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_imp") # Phi_uij - Utility flow - imported (MW)
    self.pk_util_flow_exp = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_exp") # Phi_uij - Utility flow - exported
    self.pk_util_flow_tec = self.model.addVars(self.u,self.i,self.j,    lb=-1e21,         vtype=GRB.CONTINUOUS,   name="pk_util_flow_tec") # Phi_uij - Utility flow - technologies
    self.pk_util_flow_wnd = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_wnd") # Phi_uij - Utility flow - wind technologies
    self.pk_util_flow_sto = self.model.addVars(self.u,self.i,self.j,    lb=-1e21,         vtype=GRB.CONTINUOUS,   name="pk_util_flow_sto") # Phi_uij - Utility flow - storage
    self.pk_util_flow_was = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_was") # Phi_uij - Utility flow - waste
    self.pk_util_flow_dem = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_dem") # Phi_uij - Utility flow - demand

    # Peak, positive energy flows only. Used for grid capacity constraints. Sum of all positive flows to a utility in any given increment is the grid capacity required for that increment.
    
    self.pk_util_tech_p = self.model.addVars(self.u,self.v,self.i,self.j,                   vtype=GRB.CONTINUOUS,   name="pk_util_tech_p")
    self.pk_util_wind_p = self.model.addVars(self.u,self.w,self.i,self.j,                   vtype=GRB.CONTINUOUS,   name="pk_util_wind_p")
    self.pk_util_stor_p = self.model.addVars(self.u,self.x,self.i,self.j,                   vtype=GRB.CONTINUOUS,   name="pk_util_stor_p")
    self.pk_util_flow_tec_p = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_tec_p") # Phi_uij - Utility flow - technologies
    self.pk_util_flow_wnd_p = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_wnd_p") # Phi_uij - Utility flow - wind technologies
    self.pk_util_flow_sto_p = self.model.addVars(self.u,self.i,self.j,                      vtype=GRB.CONTINUOUS,   name="pk_util_flow_sto_p") # Phi_uij - Utility flow - storage    
    


    #-----------------
    #----UTILITIES----
    #-----------------
    
    self.dv_util_cpct_int = self.model.addVars(self.u,                                        vtype=GRB.CONTINUOUS,   name="util_cpct_int") # Gamma_u - Utility connection capacity (MW)
    self.dv_util_cpct_ext = self.model.addVars(self.u,                                        vtype=GRB.CONTINUOUS,   name="util_cpct_ext") # Gamma_u - Utility connection capacity (MW)

    self.dv_util_ext_conx = self.model.addVars(self.u,                                        vtype=GRB.BINARY,       name="util_ext_conx") # gamma_u
    self.dv_util_int_conx = self.model.addVars(self.u,                                        vtype=GRB.BINARY,       name="util_int_conx") # gamma_u

    print("Decision Variables have been defined.\n")
    
    
  
  
  
  #%% Define Constraints
  
  
  def define_Constrs(self):
    print("Defining constraints...")
    
    # Constraints are generated by in-line for loops. 
    # The iterables used in the for loops should be a selection of applicable sets u,v,w,i,j, for whichever decision variable is being constrained.
    

    #------------
    #----TECH----
    #------------
    
    # Installed Capacity and Upper Production Limit
    # Nominal
    self.con_prod_limit_upper = self.model.addConstrs((self.dv_prod_tech[v,i,j] <= self.dv_num_techs[v]*self.d_tech[v]["CAPACITY"]*self.dt_j(j) for v, i, j in self.dv_prod_tech),"con_prod_limit_upper")
    
    self.con_prod_rescale = self.model.addConstrs((self.dv_prod_tech[v,i,j] == self.dv_kprod_tech[v,i,j]*1000*1000 for v, i, j in self.dv_prod_tech),"con_prod_rescale")
    
    # Peak
    self.con_prod_limit_upper_pk = self.model.addConstrs((self.pk_prod_tech[v,i,j] <= self.dv_num_techs[v]*self.d_tech[v]["CAPACITY"] for v, i, j in self.pk_prod_tech),"con_prod_limit_upper_pk")    
    
    
    #Technology Module Groups    
    for v in self.dv_num_techs:
      if self.d_tech[v]['GROUPWITH'] == 0:
        pass
      else:
        self.model.addConstr(self.dv_num_techs[v] == self.dv_num_techs[self.d_tech[v]['GROUPWITH']], "con_tech_groups")

    #Pre-existing Technologies
    for v in self.v:
      if self.d_tech[v]['NUM_EX'] == 0:
        pass
      else:
        self.model.addConstr(self.dv_num_techs[v] == self.d_tech[v]['NUM_EX'], "con_tech_existing")
    
    #Constant Power Outputs
    for v in self.v:
      for i in self.i:
        for j in self.j:
          if self.d_tech[v]['CONST_PROD'] == 'Yes':
            self.model.addConstr(self.dv_prod_tech[v,i,j] == self.dv_num_techs[v]*self.d_tech[v]["CAPACITY"]*self.dt_j(j), "con_prod_constant")

    # Energy Balance (Vertical) Across Technologies
    # Nominal (MWh)
    self.con_tech_bal = self.model.addConstrs((self.dv_util_tech[u,v,i,j] == self.mat_supr[self.u.index(u),self.v.index(v)]*self.dv_prod_tech[v,i,j] for u, v, i, j in self.dv_util_tech), "con_tech_bal")
    # Peak (MW)
    self.con_tech_bal_pk = self.model.addConstrs((self.pk_util_tech[u,v,i,j] == self.mat_supr[self.u.index(u),self.v.index(v)]*self.pk_prod_tech[v,i,j] for u, v, i, j in self.pk_util_tech), "con_tech_bal_pk")    
    
    
    
    # Energy Balance (Horizontal) Across Technologies
    # Nominal (MWh)
    self.con_util_tec = self.model.addConstrs((self.dv_util_flow_tec[u,i,j] == self.dv_util_tech.sum(u,'*',i,j) for u, i ,j in self.dv_util_flow_tec),"con_util_tec") 

    # Peak (MW)
    self.con_util_tec_pk = self.model.addConstrs((self.pk_util_flow_tec[u,i,j] == self.pk_util_tech.sum(u,'*',i,j) for u, i ,j in self.pk_util_flow_tec),"con_util_tec_pk") 
   
    
    
    #------------    
    #----WIND----
    #------------    
    
    #Installed Capacity and Production Limits (Wind Technologies)
    
    self.wind_farm_losses = 0.15
    
    # Nominal
    self.con_wind_limit_upper = self.model.addConstrs((self.dv_prod_wind[w,i,j] <= (1-self.wind_farm_losses)*self.dv_num_winds[w]*np.interp(self.d_ws['v_rmc',i,j],self.xws,self.d_windtech[w]['yp'])/1000*self.dt_j(j) for w, i, j in self.dv_prod_wind), "con_wind_limit_upper")
    
    # Peak
    self.con_wind_limit_upper_pk = self.model.addConstrs((self.pk_prod_wind[w,i,j] == 0 for w, i, j in self.pk_prod_wind), "con_wind_limit_upper_pk")

    
    # Energy Balance (Vertical) Across Wind Technologies 
    # Nominal 
    self.con_wind_bal = self.model.addConstrs((self.dv_util_wind[u,w,i,j] == self.mat_wind_supr[self.u.index(u),self.w.index(w)]*self.dv_prod_wind[w,i,j] for u, w, i, j in self.dv_util_wind),"con_wind_bal")
    
    # Peak 
    self.con_wind_bal_pk = self.model.addConstrs((self.pk_util_wind[u,w,i,j] == self.mat_wind_supr[self.u.index(u),self.w.index(w)]*self.pk_prod_wind[w,i,j] for u, w, i, j in self.pk_util_wind),"con_wind_bal_pk")    
    
    
    
    # Energy Balance (Horizontal) Across Wind Technologies
    # Nominal
    self.con_util_wnd = self.model.addConstrs((self.dv_util_flow_wnd[u,i,j] == self.dv_util_wind.sum(u,'*',i,j) for u, i ,j in self.dv_util_flow_wnd),"con_util_wnd") 
   
    # Peak 
    self.con_util_wnd_pk = self.model.addConstrs((self.pk_util_flow_wnd[u,i,j] == self.pk_util_wind.sum(u,'*',i,j) for u, i ,j in self.pk_util_flow_wnd),"con_util_wnd_pk")     
    
    
    
    
    
    
    #------------------------------   
    #----WIND CAPACITY LIMITING----
    #------------------------------
    
    if self.wind_pen_limit < 1:
      
      self.con_tech_pro = self.model.addConstrs((self.dv_util_tech_pro[u,v,i,j] == self.mat_supr_pro[self.u.index(u),self.v.index(v)]*self.dv_prod_tech[v,i,j] for u, v, i, j in self.dv_util_tech), "con_tech_pro")
      
      self.con_util_tec_pro = self.model.addConstrs((self.dv_util_prod_tec[u,i,j] == self.dv_util_tech_pro.sum(u,'*',i,j) for u, i ,j in self.dv_util_flow_tec),"con_util_tec_pro")     
  
      self.con_wind_pro = self.model.addConstrs((self.dv_util_wind_pro[u,w,i,j] == self.mat_wind_supr_pro[self.u.index(u),self.w.index(w)]*self.dv_prod_wind[w,i,j] for u, w, i, j in self.dv_util_wind),"con_wind_pro")
      
      self.con_util_wnd_pro = self.model.addConstrs((self.dv_util_prod_wnd[u,i,j] == self.dv_util_wind_pro.sum(u,'*',i,j) for u, i ,j in self.dv_util_flow_wnd),"con_util_wnd_pro")
      
      
      
      self.con_wind_pen_limit = self.model.addConstrs((self.dv_util_prod_wnd[u,i,j] <= (self.dv_util_prod_tec[u,i,j]+self.dv_util_flow_imp[u,i,j])/(1/self.wind_pen_limit-1) for u, i, j in self.dv_util_prod_wnd), "con_wind_pen_limit")

                                                

    #---------------
    #----STORAGE----
    #---------------
    
    # Net charge/discharge flow for a given storage technology for a given util for a given time period
    # Nominal
    self.con_sto_net = self.model.addConstrs((self.dv_util_stor[u,x,i,j] == self.mat_stor_supr[self.u.index(u),self.x.index(x)]*(self.dv_util_dsch[x,i,j] - self.dv_util_chrg[x,i,j]) for u, x, i, j in self.dv_util_stor),"con_sto_net")
    
    # Peak
    self.con_sto_net_pk = self.model.addConstrs((self.pk_util_stor[u,x,i,j] == self.dv_util_stor[u,x,i,j]/self.dt_j(j) for u, x, i, j in self.pk_util_stor),"con_sto_net_pk")    
    
    
    # "Either/or" constraints for charge/discharge, and maximum charge/discharge rates of storage technologies.
    print("About to do the either/or constraints:")
    self.con_sto_chrg     = self.model.addConstrs((self.dv_util_chrg[x,i,j] <= (self.dv_dir_store[x,i,j]*(self.bigM)) for x, i, j in self.dv_util_chrg), "con_sto_chrg")

    self.con_sto_chrg_amt = self.model.addConstrs((self.dv_util_chrg[x,i,j] <= (self.dv_max_store[x]*self.d_stor[x]['MAX_CH_RATE']) for x, i, j in self.dv_util_chrg), "con_sto_chrg_amt")
        
    self.con_sto_dsch     = self.model.addConstrs((self.dv_util_dsch[x,i,j] <= (1 - self.dv_dir_store[x,i,j])*(self.bigM) for x, i, j in self.dv_util_dsch), "con_sto_dsch")

    self.con_sto_dsch_amt = self.model.addConstrs((self.dv_util_dsch[x,i,j] <= (self.dv_max_store[x]*self.d_stor[x]['MAX_DS_RATE']) for x, i, j in self.dv_util_dsch), "con_sto_dsch_amt")
    
    
    # Minimum allowable level of stored energy
    self.con_sto_min_store = self.model.addConstrs((self.dv_min_store[x] == self.dv_max_store[x]*self.d_stor[x]['MIN_LIMIT'] for x in self.dv_min_store), "con_sto_min_store")
    
    # Bound the amount of stored energy between min and max settings
    self.con_sto_bound_min = self.model.addConstrs((self.dv_lvl_store[x,i,j] >= self.dv_min_store[x] for x, i, j in self.dv_lvl_store), "con_sto_bound_min")
    self.con_sto_bound_max = self.model.addConstrs((self.dv_lvl_store[x,i,j] <= self.dv_max_store[x] for x, i, j in self.dv_lvl_store), "con_sto_bound_max")
    
        
    # Stored Energy Levels
    for x in self.x:
      for i in self.i:
        for j in self.j:
          if ((i==1) and (j==1)):
            self.model.addConstr(self.dv_lvl_store[x,i,j] == self.dv_util_chrg[x,i,j]*self.d_stor[x]['EFF_CH'] - self.dv_util_dsch[x,i,j]/self.d_stor[x]['EFF_DS'], "con_sto_flow_init")
          elif (j==1):
            self.model.addConstr(self.dv_lvl_store[x,i,j] == self.dv_lvl_store[x,(i-1),len(self.j)]*(1-self.d_stor[x]['DECAY']) + self.dv_util_chrg[x,i,j]*self.d_stor[x]['EFF_CH'] - self.dv_util_dsch[x,i,j]/self.d_stor[x]['EFF_DS'], "con_sto_flow_wrap")
          else:
            self.model.addConstr(self.dv_lvl_store[x,i,j] == self.dv_lvl_store[x,i,(j-1)]*(1-self.d_stor[x]['DECAY']) + self.dv_util_chrg[x,i,j]*self.d_stor[x]['EFF_CH'] - self.dv_util_dsch[x,i,j]/self.d_stor[x]['EFF_DS'], "con_sto_flow_else")
            

    # Energy Balance (Horizontal) Across Storage Technologies
    # Nominal
    self.con_util_sto = self.model.addConstrs((self.dv_util_flow_sto[u,i,j] == self.dv_util_stor.sum(u,'*',i,j) for u, i ,j in self.dv_util_flow_sto),"con_util_sto")  
    
    # Peak
    self.con_util_sto_pk = self.model.addConstrs((self.pk_util_flow_sto[u,i,j] == self.pk_util_stor.sum(u,'*',i,j) for u, i ,j in self.pk_util_flow_sto),"con_util_sto_pk")  
        
    


    
    
    
    #--------------------------
    #----NET ENERGY BALANCE----
    #--------------------------
    
    # Nominal Utility Energy Balance Constraints Across All Technology Types, Imports, Exports, Wastes, Demands
    self.con_util_bal = self.model.addConstrs((  self.dv_util_flow_imp[u,i,j]*self.d_util[u]["IMPORTABLE"]
                                               + self.dv_util_flow_tec[u,i,j] 
                                               + self.dv_util_flow_wnd[u,i,j]
                                               + self.dv_util_flow_sto[u,i,j]
                                               - self.dv_util_flow_exp[u,i,j]*self.d_util[u]["EXPORTABLE"]
                                               - self.dv_util_flow_was[u,i,j]*self.d_util[u]["WASTEABLE"]
                                               - self.dv_util_flow_dem[u,i,j]*self.d_util[u]["DEMANDED"]
                                               == 0 for u, i, j in self.dv_util_flow_dem), "con_util_bal")
    
    
    self.con_imp_rescale = self.model.addConstrs((self.dv_util_flow_imp[u,i,j] == self.dv_kutil_flow_imp[u,i,j]*1000 for u, i, j in self.dv_util_flow_imp),"con_prod_rescale")
    

    # Peak Utility Energy Balance Constraints Across All Technology Types, Imports, Exports, Wastes, Demands
    self.con_util_bal_pk = self.model.addConstrs((   self.pk_util_flow_imp[u,i,j]*self.d_util[u]["IMPORTABLE"]  # limited only by external capacity
                                                   + self.pk_util_flow_tec[u,i,j]                               
                                                   + self.pk_util_flow_wnd[u,i,j]                               
                                                   + self.pk_util_flow_sto[u,i,j]      # how much storage is nominally being discharged that day / 24 hours
                                                   - self.pk_util_flow_exp[u,i,j]*self.d_util[u]["EXPORTABLE"]
                                                   - self.pk_util_flow_was[u,i,j]*self.d_util[u]["WASTEABLE"]
                                                   - self.pk_util_flow_dem[u,i,j]*self.d_util[u]["DEMANDED"]# user-defined peak

                                                   == 0 for u, i, j in self.pk_util_flow_dem), "con_util_bal_pk")




    
    #---------------
    #----DEMANDS---- 
    #---------------
    self.con_util_dem = self.model.addConstrs((self.dv_util_flow_dem[u,i,j] == self.d_demands[u,i,j]*self.d_util[u]['MWh/LOGICAL_UNIT'] for u, i, j in self.dv_util_flow_dem),"con_util_dem")
       
    self.con_util_dem_pk = self.model.addConstrs((self.pk_util_flow_dem[u,i,j] == self.d_peaks[u,i,j]*self.d_util[u]['MWh/LOGICAL_UNIT'] for u, i, j in self.pk_util_flow_dem),"con_util_dem_pk")





    #-----------------------
    #----GRID CAPACITIES----
    #-----------------------
    
    # External Connections:
    
    # Nominal - make sure nominal value is below the external capacity
    self.con_util_cap_imp = self.model.addConstrs((self.dv_util_flow_imp[u,i,j] <= self.d_util[u]["IMPORTABLE"]*self.dv_util_cpct_ext[u]*self.dt_j(j) for u, i, j in self.dv_util_flow_imp), "con_util_cap_imp")
    self.con_util_cap_exp = self.model.addConstrs((self.dv_util_flow_exp[u,i,j] <= self.d_util[u]["EXPORTABLE"]*self.dv_util_cpct_ext[u]*self.dt_j(j) for u, i, j in self.dv_util_flow_exp), "con_util_cap_exp")
      
    # Peak - make sure peak value is below the external capacity
    self.con_util_cap_imp_pk = self.model.addConstrs((self.pk_util_flow_imp[u,i,j] <= self.d_util[u]["IMPORTABLE"]*self.dv_util_cpct_ext[u] for u, i, j in self.pk_util_flow_imp), "con_util_cap_imp_pk")
    self.con_util_cap_exp_pk = self.model.addConstrs((self.pk_util_flow_exp[u,i,j] <= self.d_util[u]["EXPORTABLE"]*self.dv_util_cpct_ext[u] for u, i, j in self.pk_util_flow_exp), "con_util_cap_exp_pk")
    
    
    self.con_util_ext_ind = self.model.addConstrs(self.dv_util_cpct_ext[u] <= self.dv_util_ext_conx[u]*self.bigM for u in self.dv_util_cpct_ext)
    
    self.con_util_int_ind = self.model.addConstrs(self.dv_util_cpct_int[u] <= self.dv_util_int_conx[u]*self.bigM for u in self.dv_util_cpct_int)
    
    
    
    # Internal Connections:
    # Peaks - Make sure nominal value is below the internal capacity
    
    self.con_util_cap_int_tec_ub = self.model.addConstrs((self.pk_util_tech[u,v,i,j] <= self.dv_util_cpct_int[u] for u, v, i, j in self.pk_util_tech), "con_util_cap_tec_ub")
    self.con_util_cap_int_tec_lb = self.model.addConstrs((self.pk_util_tech[u,v,i,j] >= -self.dv_util_cpct_int[u] for u, v, i, j in self.pk_util_tech), "con_util_cap_tec_lb")    
        
    self.con_util_cap_int_wnd = self.model.addConstrs((self.pk_util_wind[u,w,i,j] <= self.dv_util_cpct_int[u] for u, w, i, j in self.pk_util_wind), "con_util_cap_wnd")
    
    self.con_util_cap_int_sto_ub = self.model.addConstrs((self.pk_util_stor[u,x,i,j] <= self.dv_util_cpct_int[u] for u, x, i, j in self.dv_util_stor), "con_util_cap_sto")
    
    self.con_util_cap_int_sto_lb = self.model.addConstrs((self.pk_util_stor[u,x,i,j] >= -self.dv_util_cpct_int[u] for u, x, i, j in self.dv_util_stor), "con_util_cap_sto")
    
    self.con_util_cap_int_was = self.model.addConstrs((self.pk_util_flow_was[u,i,j] <= self.dv_util_cpct_int[u] for u, i, j in self.pk_util_flow_was), "con_util_cap_was")
    self.con_util_cap_int_dem = self.model.addConstrs((self.pk_util_flow_dem[u,i,j] <= self.dv_util_cpct_int[u] for u, i, j in self.pk_util_flow_dem), "con_util_cap_dem")

    # Peaks - Make sure the sum of all positive flows is below the grid capacity. This is the true upper bound on grid capacity.
    
    self.con_tech_pos = self.model.addConstrs((self.pk_util_tech_p[u,v,i,j] == self.pk_util_tech[u,v,i,j]*self.mat_supr_pro[self.u.index(u),self.v.index(v)] for u, v, i, j in self.pk_util_tech_p),"con_tech_pos")
    self.con_wind_pos = self.model.addConstrs((self.pk_util_wind_p[u,w,i,j] == self.pk_util_wind[u,w,i,j]*self.mat_wind_supr_pro[self.u.index(u),self.w.index(w)] for u, w, i, j in self.pk_util_wind_p),"con_wind_pos")
    self.con_stor_pos = self.model.addConstrs((self.pk_util_stor_p[u,x,i,j] == self.pk_util_stor[u,x,i,j]*self.mat_stor_supr_pro[self.u.index(u),self.x.index(x)] for u, x, i, j in self.pk_util_stor_p),"con_stor_pos")
    
       
    
    self.con_tech_sum_pos = self.model.addConstrs((self.pk_util_flow_tec_p[u,i,j] == self.pk_util_tech_p.sum(u,'*',i,j) for u, v, i, j in self.pk_util_tech_p), "con_tech_sum_pos")
    
    self.con_wind_sum_pos = self.model.addConstrs((self.pk_util_flow_wnd_p[u,i,j] == self.pk_util_wind_p.sum(u,'*',i,j) for u, w, i, j in self.pk_util_wind_p), "con_wind_sum_pos")
    
    self.con_stor_sum_pos = self.model.addConstrs((self.pk_util_flow_sto_p[u,i,j] == self.pk_util_stor_p.sum(u,'*',i,j) for u, x, i, j in self.pk_util_stor_p), "con_stor_sum_pos")

    self.con_util_cap_int_flow = self.model.addConstrs((self.pk_util_flow_imp[u,i,j] + self.pk_util_flow_tec_p[u,i,j] + self.pk_util_flow_wnd_p[u,i,j] + self.pk_util_flow_sto_p[u,i,j] <= self.dv_util_cpct_int[u] for u, i, j in self.pk_util_flow_imp), "con_util_cap_int_flow")


    
    print("Constraints have been defined.\n")
    
    #Model Update
    self.model.update()
    

    
  
    
    


#%% Define Objective Statements
  def define_Obj(self):
    
    print("Defining Objective Statements...")
    
    
    # The current strategy is to minimize the present value of all costs.
    ##   All cash flows are assumed to be disbursed at the end of periods i, other than capital costs which are disbursed at time zero. The discount rate input is compounded at the end of each period i. 
    
    #-------------------
    #----FIXED COSTS----
    #-------------------
    
    
    self.c_fix_techCAPEX = quicksum(self.dv_num_techs[v]*self.d_tech[v]['CAPEX']/1000 for v in self.dv_num_techs) 
    self.c_fix_windCAPEX = quicksum(self.dv_num_winds[w]*self.d_windtech[w]['CAPEX']/1000 for w in self.dv_num_winds) 
      
    self.c_fix_storCAPEX = quicksum(self.dv_max_store[x]*self.d_stor[x]['SP_CAPEX']/1000 for x in self.dv_max_store)

    self.c_fix_utilCAPEX_ext =  quicksum(self.dv_util_cpct_ext[u]*(self.d_util[u]['CAPEX_a_EXT']/1000+self.d_util[u]['CAPEX_b_EXT']/1000*self.d_util[u]['D_TRANS_EXT']) + self.d_util[u]['CAPEX_c_EXT']/1000*self.d_util[u]['D_TRANS_EXT']*self.dv_util_ext_conx[u] for u in self.dv_util_cpct_ext)
    self.c_fix_utilCAPEX_int =  quicksum(self.dv_util_cpct_int[u]*(self.d_util[u]['CAPEX_a_INT']/1000+self.d_util[u]['CAPEX_b_INT']/1000*self.d_util[u]['D_TRANS_INT']) + self.d_util[u]['CAPEX_c_INT']/1000*self.d_util[u]['D_TRANS_INT']*self.dv_util_int_conx[u] for u in self.dv_util_cpct_int)
    
    

    self.c_fix_techOPEX_fix = quicksum(self.dv_num_techs[v]*self.d_tech[v]['OPEX_FIX']/1000*self.econ_PAIN_at[self.ni] for v in self.dv_num_techs)
    self.c_fix_windOPEX_fix = quicksum(self.dv_num_winds[w]*self.d_windtech[w]['OPEX_FIX']/1000*self.econ_PAIN_at[self.ni] for w in self.dv_num_winds)

    self.c_fix = self.c_fix_techCAPEX + self.c_fix_windCAPEX + self.c_fix_storCAPEX + self.c_fix_utilCAPEX_ext + self.c_fix_utilCAPEX_int + self.c_fix_techOPEX_fix + self.c_fix_windOPEX_fix
    
    #----------------------
    #----VARIABLE COSTS----
    #----------------------
    
    self.c_var_imports  = quicksum(self.dv_kutil_flow_imp[u,i,j]*self.d_import_prices[u,i,j]*1000/1000/self.d_util[u]['MWh/LOGICAL_UNIT']*self.econ_PFIN_at[i] for u,i,j in self.dv_util_flow_imp)
    self.c_var_exports  = quicksum(self.dv_util_flow_exp[u,i,j]*self.d_export_prices[u,i,j]/1000/self.d_util[u]['MWh/LOGICAL_UNIT']*self.econ_PFIN_at[i] for u,i,j in self.dv_util_flow_exp)
    self.c_var_techOPEX_var = quicksum(self.dv_kprod_tech[v,i,j]*self.d_tech[v]['OPEX_VAR']*1000*1000/1000*self.econ_PFIN_at[i] for v,i,j in self.dv_prod_tech)
    self.c_var_windOPEX_var = quicksum(self.dv_prod_wind[w,i,j]*self.d_windtech[w]['OPEX_VAR']/1000*self.econ_PFIN_at[i] for w,i,j in self.dv_prod_wind)
    
    self.c_var = (self.c_var_imports - self.c_var_exports) + self.c_var_techOPEX_var + self.c_var_windOPEX_var
    
    
    #-------------------
    #----TOTAL COSTS----
    #-------------------
    
    self.c_tot = (self.c_fix + self.c_var)
    
    self.model.setObjective(self.c_tot, GRB.MINIMIZE)
    
    
    print("Objective statements created.\n")
  
  #%% Optimize the Model
  def run_optimization(self):
    print('--------------------------------')
    print("Commencing model optimization...")
    self.model.update()
    self.model.optimize()

    self.gen_supr_at_all()
    self.gen_supr_at_i()
    self.maxflow = self.get_max_dv_util_flow()
    self.maxflow_i = self.get_max_dv_util_flow_i()
    print('--------------------------------')
    print("-----Optimization complete------")
    print('--------------------------------')
    print('')
    if not self.bigM_is_valid():
      print("------------------------PROBLEM---------------------------------")
      print("-      bigM value has been reached - increase value of bigM    -")
      print("----------------------------------------------------------------")      
    else:
      print("Notice: bigM values were adequately large and did not overconstrain the model.")
    

#%% Run everything
 
  def run_model(self):
    self.import_Data()
    self.define_DVs()
    self.define_Constrs()
    self.define_Obj()
    self.run_optimization()
        
#%% Run everything and export excel results
 
  def run_all(self):
    self.run_model()
    self.get_cashflows()
    self.get_results()
    
#%% Result Array Generating Functions 
  
  def get_imports_array(self):
    array = np.zeros((len(self.u),self.ij))
    
    for u in self.u:
      for i in self.i:
        for j in self.j:
          array[self.u.index(u)][(i-1)*len(self.j)+(j-1)] = self.dv_util_flow_imp[u,i,j].X
          
    return array
  
  def get_exports_array(self):
    array = np.zeros((len(self.u),self.ij))
    
    for u in self.u:
      for i in self.i:
        for j in self.j:
          array[self.u.index(u)][(i-1)*len(self.j)+(j-1)] = self.dv_util_flow_exp[u,i,j].X
          
    return array
    
  
  def get_wastes_array(self):
    array = np.zeros((len(self.u),self.ij))
    
    for u in self.u:
      for i in self.i:
        for j in self.j:
          array[self.u.index(u)][(i-1)*len(self.j)+(j-1)] = self.dv_util_flow_was[u,i,j].X
          
    return array        
  
  def get_demands_array(self):
    array = np.zeros((len(self.u),self.ij))
    
    for u in self.u:
      for i in self.i:
        for j in self.j:
          array[self.u.index(u)][(i-1)*len(self.j)+(j-1)] = self.dv_util_flow_dem[u,i,j].X
          
    return array         
        
  def get_production_array(self):
    array = np.zeros((len(self.v)+len(self.w),self.ij))
    
    for v in self.v:
      for i in self.i:
        for j in self.j:
          array[self.v.index(v)][(i-1)*len(self.j)+(j-1)] = self.dv_prod_tech[v,i,j].X
    for w in self.w:
      for i in self.i:
        for j in self.j:
          array[len(self.v)+self.w.index(w)][(i-1)*len(self.j)+(j-1)] = self.dv_prod_wind[w,i,j].X
  
    return array             

  def get_storage_array(self):
    array = np.zeros((len(self.x),self.ij))
  
    for x in self.x:
      for i in self.i:
        for j in self.j:
          array[self.x.index(x)][(i-1)*len(self.j)+(j-1)] = self.dv_lvl_store[x,i,j].X  
       
    return array       
  
  def get_storage_ch_array(self,x):
    array = np.zeros((4,self.ij))
    
    for i in self.i:
      for j in self.j:
        array[0,(i-1)*len(self.j)+(j-1)] = self.dv_lvl_store[x,i,j].X   #storage level
        array[1,(i-1)*len(self.j)+(j-1)] = self.dv_util_chrg[x,i,j].X   #charge value
        array[2,(i-1)*len(self.j)+(j-1)] = -self.dv_util_dsch[x,i,j].X  #discharge value
        array[3,(i-1)*len(self.j)+(j-1)] = self.dv_dir_store[x,i,j].X  #direction binary
    
    return array

  def get_util_array(self,u):
    util = self.u.index(u)
    rows = len(self.v+self.w+self.x)+4
    cols = self.ij
    array = np.zeros((rows, cols))
    for r in range(rows):
      for i in self.i:
        for j in self.j:
          array[r][(i-1)*len(self.j)+(j-1)] = self.d_suprs[tuple((i,j))][util][r]     
          
    
          
    return array
          
#%% Get cashflows of v, w, x, u, imports, exports, for period j


  def get_cf_techCAPEX(self,v):
    cf = np.zeros(len(self.i)) 
    cap = self.dv_num_techs[v].X*self.d_tech[v]['CAPEX']
    cf = np.insert(cf,0,cap) 
    return cf
    
  def get_cf_windCAPEX(self,w):
    cf = np.zeros(len(self.i)) 
    cap = self.dv_num_winds[w].X*self.d_windtech[w]['CAPEX']
    cf = np.insert(cf,0,cap)     
    return cf
    
  def get_cf_storCAPEX(self,x):
    cf = np.zeros(len(self.i)) 
    cap = self.dv_max_store[x].X*self.d_stor[x]['SP_CAPEX']
    cf = np.insert(cf,0,cap)     
    return cf
    
  def get_cf_utilCAPEX_ext(self,u):
    cf = np.zeros(len(self.i)) 
    cap = self.dv_util_cpct_ext[u].X*(self.d_util[u]['CAPEX_a_EXT']+self.d_util[u]['CAPEX_b_EXT']*self.d_util[u]['D_TRANS_EXT']) + self.d_util[u]['CAPEX_c_EXT']*self.d_util[u]['D_TRANS_EXT']*self.dv_util_ext_conx[u].X
    cf = np.insert(cf,0,cap) 
    return cf
  
  def get_cf_utilCAPEX_int(self,u):
    cf = np.zeros(len(self.i)) 
    cap = self.dv_util_cpct_int[u].X*(self.d_util[u]['CAPEX_a_INT']+self.d_util[u]['CAPEX_b_INT']*self.d_util[u]['D_TRANS_INT']) + self.d_util[u]['CAPEX_c_INT']*self.d_util[u]['D_TRANS_INT']*self.dv_util_int_conx[u].X
    cf = np.insert(cf,0,cap)
    return cf    
    
  def get_cf_techOPEX_fix(self,v):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = self.dv_num_techs[v].X*self.d_tech[v]['OPEX_FIX']    
    return cf        
    
  def get_cf_windOPEX_fix(self,w):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = self.dv_num_winds[w].X*self.d_windtech[w]['OPEX_FIX']   
    return cf
  
  def get_cf_imports(self,u):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = 0
      for j in self.j:
          cf[i] = cf[i] + self.dv_util_flow_imp[u,i,j].X*self.d_import_prices[u,i,j]/self.d_util[u]['MWh/LOGICAL_UNIT']
    
    return cf 
    
  def get_cf_exports(self,u):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = 0
      for j in self.j:
          cf[i] = cf[i] - self.dv_util_flow_exp[u,i,j].X*self.d_export_prices[u,i,j]/self.d_util[u]['MWh/LOGICAL_UNIT']
    
    return cf 
  
  def get_cf_techOPEX_var(self,v):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = 0
      for j in self.j:
          cf[i] = cf[i] + self.dv_prod_tech[v,i,j].X*self.d_tech[v]['OPEX_VAR']
    
    return cf 
    
  def get_cf_windOPEX_var(self,w):
    cf = np.zeros(len(self.i)+1) 
    for i in self.i:
      cf[i] = 0
      for j in self.j:
          cf[i] = cf[i] + self.dv_prod_wind[w,i,j].X*self.d_windtech[w]['OPEX_VAR']
    
    return cf
  
  def get_cashflows(self):
    cashflows = {}
    
    #Utilities
    for u in self.u:
      cashflows[u+"_CAPEX_INT"] = self.get_cf_utilCAPEX_int(u) 
      cashflows[u+"_CAPEX_EXT"] = self.get_cf_utilCAPEX_ext(u)
      cashflows[u+"_IMPORTS"] = self.get_cf_imports(u)
      cashflows[u+"_EXPORTS"] = self.get_cf_exports(u)  
      
    #Technologies
    for v in self.v:
      cashflows[v+"_CAPEX"] = self.get_cf_techCAPEX(v)
      cashflows[v+"_OPEX_FIX"] = self.get_cf_techOPEX_fix(v)
      cashflows[v+"_OPEX_VAR"] = self.get_cf_techOPEX_var(v)
      
    
    #Wind Technologies
    for w in self.w:
      cashflows[w+"_CAPEX"] = self.get_cf_windCAPEX(w)
      cashflows[w+"_OPEX_FIX"] = self.get_cf_windOPEX_fix(w)
      cashflows[w+"_OPEX_VAR"] = self.get_cf_windOPEX_var(w)    
    
    #Storage Technologies
    for x in self.x:
      cashflows[x+"_CAPEX"] = self.get_cf_storCAPEX(x)
      
    cashflows_data = pd.DataFrame.from_dict(cashflows, orient = 'index')
    cashflows_data.to_excel("_cashflows_" + self.modelname+ ".xlsx",sheet_name = "cashflows")
    print(self.modelname +" cashflows have been exported")
  
  def get_results(self):
    results = {}
      
    #Utilities
    for u in self.u:
      results["CAPACITY_INT"+u] = self.dv_util_cpct_int[u].X
      results["CAPACITY_EXT"+u] = self.dv_util_cpct_ext[u].X
      
    #Technologies
    for v in self.v:
      results["num_techs_"+v] = self.dv_num_techs[v].X
    
    #Wind Technologies
    for w in self.w:
      results["num_winds_"+w] = self.dv_num_winds[w].X
  
    #Storage Technologies
    for x in self.x:
      results["max_store_"+x] = self.dv_max_store[x].X
    
    costs = self.get_costs()
    
    results.update(costs)
      
    results_data = pd.DataFrame.from_dict(results, orient = 'index')
    results_data.to_excel("_results_"+ self.modelname +".xlsx", sheet_name = "results")
    
    storageFlows = {}
    
    for x in self.x:
      frame = pd.DataFrame(np.transpose(self.get_storage_ch_array(x)))
      frame.to_excel("_storageflow_" + self.modelname + "_" + x + ".xlsx",sheet_name = "results")    
    
    print(self.modelname + " results have been exported.")













