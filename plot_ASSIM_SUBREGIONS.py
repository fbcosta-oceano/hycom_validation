from netCDF4 import Dataset
from scipy import interpolate
import scipy.io
import remo
import datetime
import sys
import os
import time 
import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature

pd.options.display.float_format = '{:.2f}'.format
np.set_printoptions(threshold=sys.maxsize)

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
#data_inicial = str(content_list[3].rstrip('\n'))
#data_final = str(content_list[4].rstrip('\n'))
#inc_tempo = float(content_list[5].rstrip('\n'))

output_dir = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/resultados/figs'
dir_hycom = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/previsao/hycom_2_2_18/proc'
dir_assim = 'assim/slasst_gridded_map_error_or_argoZ_stats/Input'
sub_file = 'subregions_008i_0256'

if '004j' in rodada[0]:
   IDM = 502
   JDM = 564
   dir_ATL = dir_hycom+'/ATLj0.04/'

   lat_step_plot = 4
   lon_step_plot = 4
   resolution_plot = 300
elif '008d' in rodada[0]:
   IDM = 1717
   JDM = 2345
   dir_ATL = dir_hycom+'/ATLd0.08/'

   lat_step_plot = 5
   lon_step_plot = 5
   resolution_plot = 300
elif '008i' in rodada[0]:
   IDM = 628
   JDM = 780
   dir_ATL = dir_hycom+'/ATLi0.08/'

   lat_step_plot = 5
   lon_step_plot = 5
   resolution_plot = 300
else:
   print('CASE '+rodada[0]+' NOT DEFINED')
   exit()

dir_expt = dir_ATL+expt[0]
IJDM=IDM*JDM
npad=4096-(IJDM%4096)
counter = -1

f = open(dir_ATL+'/topo/regional.grid.a','rb')
f.seek(0*4*(npad+IJDM))
field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
lon_hycom = np.reshape(field,(JDM,IDM))
lon_hycom = lon_hycom[0,:]
f.seek(1*4*(npad+IJDM))
field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
lat_hycom = np.reshape(field,(JDM,IDM))
lat_hycom = lat_hycom[:,0]
f.close()

max_lon_hycom = max(lon_hycom)
min_lon_hycom = min(lon_hycom)
max_lat_hycom = max(lat_hycom)
min_lat_hycom = min(lat_hycom)

opt=open(dir_expt+'/'+dir_assim+'/'+sub_file+'.ascii', "r")
content_list = opt.readlines()
#print(len(test))
#print(test[0],int(test[4]))
#print(content_list[0].rstrip('\n').split(','),len(content_list))
#exit()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([min_lon_hycom, max_lon_hycom, min_lat_hycom, max_lat_hycom], crs=ccrs.PlateCarree())

for ind in range(0,len(content_list),1):
#for ind in range(0,2,1):
    sub = str(content_list[ind].rstrip('\n').split(','))
    sub = sub.replace('[\'','').replace('\']','').split(' ')
    x0 = lon_hycom[int(sub[1])-1]
    x1 = lon_hycom[int(sub[2])-1]
    y0 = lat_hycom[int(sub[3])-1]
    y1 = lat_hycom[int(sub[4])-1]

    plt.plot([x0, x1], [y0, y0], color='black', linewidth=0.7)
    plt.plot([x0, x1], [y1, y1], color='black', linewidth=0.7)
    plt.plot([x0, x0], [y0, y1], color='black', linewidth=0.7)
    plt.plot([x1, x1], [y0, y1], color='black', linewidth=0.7)

states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

parallels = np.arange(int(min_lat_hycom),int(max_lat_hycom)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon_hycom)+2,int(max_lon_hycom)-1,lon_step_plot)
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

plt.savefig(output_dir+'/'+sub_file+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
exit()


