import datetime
import time
import remo
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
import sys
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings("ignore")

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

#IDM=1717
#JDM=2345

IDM=628
JDM=780
path='/home/filipe.costa/previsao/hycom_2_2_18/proc/ATLi0.08'
year=2018
lat_step_plot = 10
lon_step_plot = 10

#IDM=502
#JDM=564

IJDM=IDM*JDM
npad=4096-(IJDM%4096)

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial

record = 0
while(current_data <= data_final):
    if (current_data==data_inicial):
       f = open(path+'/topo/regional.grid.a','rb')
       f.seek(0*4*(npad+IJDM))
       field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
       lon = np.reshape(field,(JDM,IDM))
       lon = lon[0,:]
       f.seek(1*4*(npad+IJDM))
       field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
       lat = np.reshape(field,(JDM,IDM))
       lat = lat[:,0]
       f.close()
       
       max_lat=max(lat)
       min_lat=min(lat)
       min_lon=min(lon)
       max_lon=max(lon)

       parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
       meridians = np.arange(int(min_lon),int(max_lon)+0.1,lon_step_plot)
       under='#000033'
       over='#330000'
       bad='#FFFFFF'
       cmap=matplotlib.cm.get_cmap('jet')
       cmap.set_under(under)
       cmap.set_over(over)
       cmap.set_bad(bad)
    
    for rec in range(0,4,1):
        f = open(path+'/force/cfsr/'+str(year)+'/airtmp.a','rb')
        f.seek(record*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        field = np.reshape(field,(JDM,IDM))
        f.close()
        
        f = open(path+'/force/cfsr/'+str(year)+'_PERT_0/airtmp.a','rb')
        f.seek(record*4*(IJDM+npad))
        field2 = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        cond = abs(field2) > 999999
        field2 = np.ma.masked_where(cond,field2)
        field2 = np.reshape(field2,(JDM,IDM))
        f.close()
        
        fig, ax = plt.subplots(nrows=1,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        ax = ax.flatten()
        
        ax[0].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
        im = ax[0].contourf(lon, lat, field, np.arange(10,30.5,.5),cmap=cmap, \
                         vmin=10,vmax=30,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
        
        ax[0].add_feature(cfeature.LAND)
        ax[0].add_feature(cfeature.COASTLINE, linewidth=0.5)
        
        states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
        ax[0].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
        
        gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(meridians)
        gl.ylocator = mticker.FixedLocator(parallels)
        gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        ax[0].text(-65, -8, 'NO PERT', horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())
        
        
        ax[1].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
        im = ax[1].contourf(lon, lat, field2, np.arange(10,30.5,.5),cmap=cmap, \
                         vmin=10,vmax=30,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
        
        ax[1].add_feature(cfeature.LAND)
        ax[1].add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax[1].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
        
        gl = ax[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.right_labels = False
        gl.left_labels = False
        gl.xlocator = mticker.FixedLocator(meridians)
        gl.ylocator = mticker.FixedLocator(parallels)
        gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        ax[1].text(-65, -8, 'PERT', horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

        ax[2].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
        im = ax[2].contourf(lon, lat, field2-field, np.arange(-3.0,3.25,.25),cmap=cmap, \
                         vmin=-3,vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')

        ax[2].add_feature(cfeature.LAND)
        ax[2].add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax[2].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

        gl = ax[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.right_labels = False
        gl.left_labels = False
        gl.xlocator = mticker.FixedLocator(meridians)
        gl.ylocator = mticker.FixedLocator(parallels)
        gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        ax[2].text(-65, -8, 'DIFF', horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())
        
        
        #cbar_ticks = np.arange(10,32,2.0)
        cbar_ticks = np.arange(-3.0,3.25,0.25)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
        cbar_ax = fig.add_axes([0.99, 0.15, 0.02, 0.7])
        cbar=fig.colorbar(im, ticks = cbar_ticks, cax=cbar_ax,orientation='vertical')
        plt.suptitle('AIRTMP ('+current_data.strftime("%Y%m%d")+' RECORD: '+str(record)+')', y=0.90)
        
        fig.savefig('AIRTMP_'+current_data.strftime("%Y%m%d")+'_RECORD-'+str(record),dpi=300,transparent=False,bbox_inches='tight',pad_inches=0.05)
        record = record+1
    current_data = current_data + datetime.timedelta(days=inc_tempo)

