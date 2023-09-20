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
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

output_dir = '/home/filipe/resultados/figs'
dir_hycom = '/disco1/remo/data/hycom_ufba'
dir_obs = '/disco1/remo/data'

obs_types = ['CTD','MRB','XBT','GLD','APB','OSD','DRB']
#obs_types = ['MRB','XBT','GLD']
#obs_types = ['XBT']

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

if '004j' in rodada[0]:
   IDM = 502
   JDM = 564
   dir_ATL = dir_hycom+'/ATLj0.04/'

   lat_step_plot = 4
   lon_step_plot = 4
   point_color = ['dodgerblue','red','green','blue','yellow']
   resolution_plot = 300
elif '008d' in rodada[0]:
   IDM = 1717
   JDM = 2345
   dir_ATL = dir_hycom+'/ATLd0.08/'

   lat_step_plot = 14
   lon_step_plot = 14
   point_color = ['dodgerblue','red','green','black','yellow']
   point_color = ['dodgerblue','red','green','yellow','blue','magenta','tan','gray']
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

#max_lon_hycom = 179
#min_lon_hycom = -179
#max_lat_hycom = max(lat_hycom)
#min_lat_hycom = min(lat_hycom)

print(len(obs_types))
for obs in range(0,len(obs_types),1):
    current_data = data_inicial
    print(obs_types[obs])
    arq_obs = dir_obs+'/'+obs_types[obs]+'/'+obs_types[obs]+'S'+current_data.strftime("%Y")+'_TEMP_FORMATTED.ascii'
    if not (os.path.isfile(arq_obs)):
       arq_obs = dir_obs+'/'+obs_types[obs]+'/'+obs_types[obs]+'O'+current_data.strftime("%Y")+'_TEMP_FORMATTED.ascii'
    if not (os.path.isfile(arq_obs)):
        print('FILE FOR OBS: '+obs_types[obs]+' YEAR: '+current_data.strftime("%Y")+' NOT FOUND')
        continue

    opt=open(arq_obs, "r")
    content_list = opt.readlines()
    temp = np.ma.array(np.zeros((len(content_list),7)), mask=True)
    for i in range(0,len(content_list),1):
        cont = content_list[i].rstrip('\n').split(' ')
        cont = list(filter(None,cont))
        cont = np.array(cont)
        temp[i] = cont[:]
    
    del opt, content_list, cont

    cond = (temp[:,0] < max_lon_hycom) & (temp[:,0] > min_lon_hycom) & \
           (temp[:,1] < max_lat_hycom) & (temp[:,1] > min_lat_hycom)
    temp_work = temp[cond,:]
    if not cond.any():
       print('NO '+obs_types[obs]+' FOUND INSIDE DOMAIN')
       continue
    del cond

    while(current_data <= data_final):

        cond = (temp_work[:,5] == float(current_data.strftime("%d"))) & (temp_work[:,3] == float(current_data.strftime("%m"))) & \
               (temp_work[:,2] == float(current_data.strftime("%Y")))
        if not cond.any():
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        if 'temp_plot' in locals():
           temp_plot = np.ma.concatenate((temp_plot,temp_work[cond,:]),axis=0)
        else:
           temp_plot = temp_work[cond,:]
        del cond

        current_data = current_data + datetime.timedelta(days=inc_tempo)

    if not ('temp_plot' in locals()):
       print('NO '+obs_types[obs]+' FOUND IN PERIOD')
       continue
    dif_lon = np.diff(temp_plot[:,0])
    dif_lat = np.diff(temp_plot[:,1])
    cond = (abs(dif_lon)>0.1) | (abs(dif_lat)>0.1)
    if not cond.any():
       ind=0
    #   continue
    else:
       ind = [i for i, val in enumerate(cond) if val]
       ind.insert(0,-1)

    if not ('ax' in locals()):
       ax = plt.axes(projection=ccrs.PlateCarree())
       ax.set_extent([min_lon_hycom, max_lon_hycom, min_lat_hycom, max_lat_hycom], crs=ccrs.PlateCarree())

    print(obs_types[obs])
    plt.scatter(x=temp_plot[ind,0],y=temp_plot[ind,1],color=point_color[obs], s=1, alpha=0.5, transform=ccrs.PlateCarree(),label=obs_types[obs])
    del temp_work, temp_plot, temp

ax.add_feature(cfeature.GSHHSFeature(scale='high'))
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)

#states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
#ax.add_feature(states_provinces, edgecolor='gray')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.right_labels = False

parallels = np.arange(int(min_lat_hycom),int(max_lat_hycom)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon_hycom)+2,int(max_lon_hycom)-1,lon_step_plot)
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

plt.legend()
plt.savefig('OBS_POINT_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
exit()


