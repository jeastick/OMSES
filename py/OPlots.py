# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 19:17:43 2020

@author: Jeff Eastick
"""

import numpy as np
import matplotlib as mpl
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

from matplotlib import cm
from collections import OrderedDict

cmaps = OrderedDict()

mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "CMU Serif"
mpl.rcParams['axes.unicode_minus']=False

lw = 0.3

def plt_imports(model):
  y = model.get_imports_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  for row in range(len(model.u)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label= model.u[row], linewidth = lw)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Imports')
  plt.show()

def plt_exports(model):
  y = model.get_exports_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  for row in range(len(model.u)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label= model.u[row], linewidth = lw)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Exports')
  plt.show()
  
def plt_wastes(model):
  y = model.get_wastes_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  for row in range(len(model.u)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label= model.u[row], linewidth = lw)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Wastage')
  plt.show()  
  
def plt_demands(model):
  y = model.get_demands_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  for row in range(len(model.u)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label= model.u[row], linewidth = lw)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Demands')
  plt.show()  

def plt_production(model):
  y = model.get_production_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  labels = model.v+model.w
  for row in range(len(model.v+model.w)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label = labels[row], linewidth = lw)
  # for row in range(len(model.w)):
  #   ax.plot(x, y[row], label = model.w[row])
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Production of Characteristic Utility by Technology')
  plt.show()

def plt_storage(model):
  y = model.get_storage_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  labels = model.x
  for row in range(len(model.x)):
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row], label = labels[row], linewidth = lw)

  ax.legend(loc='right')
  plt.ylabel('MWh')
  plt.title('Storage Level of Characteristic Utility by Storage Technology')
  plt.show()

def plt_storage_ch_x(model,x):
  y = model.get_storage_ch_array(x)
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  labels = ['Storage Level', 'MWh Charged', 'MWh Discharged']
  for row in range(3):
    ax.plot(x, y[row], label = labels[row], linewidth = lw)

  ax.legend(loc='right')
  plt.ylabel('MWh')
  plt.title('Storage Level of Characteristic Utility by Storage Technology')
  plt.show()


def plt_imports_stacked(model):
  y = model.get_imports_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  ax.stackplot(x,y,labels=model.u)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Imports (Stacked)')
  plt.show()

def plt_exports_stacked(model):
  y = model.get_exports_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  ax.stackplot(x,y,labels=model.u)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Exports (Stacked)')
  plt.show()
  
def plt_wastes_stacked(model):
  y = model.get_wastes_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  ax.stackplot(x,y,labels=model.u)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Wastage (Stacked)')
  plt.show()  

def plt_demands_stacked(model):
  y = model.get_demands_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  ax.stackplot(x,y,labels=model.u)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Utility Demands (Stacked)')
  plt.show()  

def plt_production_stacked(model):
  y = model.get_production_array()
  x = np.arange(model.ij)
  fix, ax= plt.subplots()
  ax.stackplot(x,y,labels=model.v+model.w)
  ax.legend(loc='right')
  plt.ylabel('MWh/period')
  plt.title('Production of Characteristic Utility by Technology (Stacked)')
  plt.show()
  
  
def plt_util_flow(model,u):
  y = model.get_util_array(u)
  x = np.arange(model.ij)
  labels = model.v.copy() + model.w.copy() + model.x.copy()
  labels.insert(0,"Imports")
  labels.append("Exports")
  labels.append("Waste")
  labels.append("Demand")
  
  colour = iter(cm.rainbow(np.linspace(0,1,len(labels))))
  
  fig, ax= plt.subplots()
  fig.set_size_inches(9, 6)

  for row in range(len(labels)):
    c = next(colour)
    if (np.count_nonzero(y[row]) > 0):
      ax.plot(x, y[row],c=c, label = labels[row], linewidth = lw)
  ax.legend(loc='right')
  
  def ijformat(x,y):
    j = int((x%model.nj)+1)
    i = int(math.ceil((x+1)/model.nj))
    return "({i},{j})".format(i=i,j=j)
  

  ax.xaxis.set_ticks(np.arange(0, model.ij, model.nj))  
  
  
  ax.xaxis.set_major_formatter(tick.FuncFormatter(ijformat))
  
  plt.xlim(0,model.ij)
  plt.xlabel('Period (i,j)')
  plt.ylabel('MWh/period')
  plt.title('Utility Balance - ' + u)
  fig.tight_layout()

  plt.show()
  
  
### MATPLOTLIB Function that takes the results from an OMSES model , and prints a heatmap of the superstructure for period i,j

def plot_structure_at(model,i,j):
    #get the import/tech/export/waste/demand data for a period (i,j) into a dataframe or numpy array
    #plot that numpy array using a heatmap
  ss = model.get_supr_at(i,j).astype(int)
    
  list_head = ["Imports"]
  list_tail = ["Exports","Waste","Demands"]
  
  utils = model.u
  techs = list_head + model.v + model.w + model.x + list_tail
  
  fig, ax = plt.subplots()
  im = ax.imshow(ss , cmap = plt.cm.RdYlGn, vmin = -model.maxflow, vmax = model.maxflow)
  ax.set_xticks(np.arange(len(techs)))
  ax.set_yticks(np.arange(len(utils)))
  ax.set_xticklabels(techs)
  ax.set_yticklabels(utils)  

  ax.tick_params(top=True, bottom=False,
               labeltop=True, labelbottom=False)
  plt.setp(ax.get_xticklabels(), rotation=60, ha="left",
         rotation_mode="anchor")
  
  for u in range(len(utils)):
      for v in range(len(techs)):
          text = ax.text(v, u, ss[u, v],
                         ha="center", va="center", color="k")          

  ax.set_ylim(len(utils)-0.5, -0.5)
  ax.set_title("Energy Balance in MWh at i = " + str(i) + ", j = " + str(j))
  fig.tight_layout()
  plt.show()
  
  
def plot_structure_at_i(model,i):
  ss = model.d_suprs_i[i].astype(int)
  
  list_head = ["Imports"]
  list_tail = ["Exports","Waste","Demands"]
  
  utils = model.u
  techs = list_head + model.v + model.w +model.x + list_tail
  
  fig, ax = plt.subplots()
  fig.set_size_inches(12,9)
  im = ax.imshow(ss , cmap = plt.cm.RdYlGn, vmin = -model.maxflow_i, vmax = model.maxflow_i)
  ax.set_xticks(np.arange(len(techs)))
  ax.set_yticks(np.arange(len(utils)))
  ax.set_xticklabels(techs)
  ax.set_yticklabels(utils)  

  ax.tick_params(top=True, bottom=False,
               labeltop=True, labelbottom=False)
  plt.setp(ax.get_xticklabels(), rotation=60, ha="left",
         rotation_mode="anchor")
  
  for u in range(len(utils)):
      for v in range(len(techs)):
          text = ax.text(v, u, ss[u, v],
                         ha="center", va="center", color="k")          

  ax.set_ylim(len(utils)-0.5, -0.5)
  ax.set_title("Annual nominal energy balance in MWh in year 10, with wind specific capital cost = $3,300/kW and a 70% reduction in SMR CAPEX")
  fig.tight_layout()
  # fig.colorbar(im)

  plt.show()  
    
  


def plot_structure_animation(model):
  print("starting animation function")
  fig, ax = plt.subplots()
  list_head = ["Imports"]
  list_tail = ["Exports","Waste","Demands"]
 
  utils = model.u
  techs = list_head + model.v + list_tail  
  ax.set_xticks(np.arange(len(techs)))
  ax.set_yticks(np.arange(len(utils)))
  ax.set_xticklabels(techs)
  ax.set_yticklabels(utils)
  ax.tick_params(top=True, bottom=False,
               labeltop=True, labelbottom=False)
  plt.setp(ax.get_xticklabels(), rotation=60, ha="left",
         rotation_mode="anchor")    
  fig.tight_layout()
  ax.set_title('System Superstructure Energy Flows [MWh] for period')
  ims = []
  print("starting for loop")
  for frame in range(model.ij):
    i = math.ceil((frame+1)/model.nj)
    j = frame % model.nj+1
    ss = model.d_suprs[(i,j)].astype(int)
    im = ax.imshow(ss, cmap = plt.cm.RdYlGn, animated = True , vmin = -model.maxflow, vmax = model.maxflow)
    ti = ax.text(0.5,1.45,("i = " + str(i) + ", j = " + str(j)),
                 size = plt.rcParams["axes.titlesize"],
                 ha = "center", transform = ax.transAxes)
    art = [im,ti]

    for u in range(len(utils)):
      for v in range(len(techs)):
          art.append(ax.text(v, u, ss[u, v],
                          ha="center", va="center", color="k"))
    

    ims.append(art)
    print("did i = " + str(i) + " j = " + str(j))

  
  ani = animation.ArtistAnimation(fig, ims, interval = 50, blit = False, repeat_delay = 1000)
  print("animation created, saving to video")
  ani.save("movie.mp4")
  print("running plt.show()")

  plt.show()
  
  
  
  