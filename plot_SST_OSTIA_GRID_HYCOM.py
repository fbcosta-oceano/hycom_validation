import remo
import datetime
import time
import sys
import os
import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings("ignore")

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

input_dir = '/home/filipe.costa/resultados'
fig_output_dir = '/home/filipe.costa/resultados/figs'

month_label_interval = 1

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

under='#000033'
over='#330000'
bad='#FFFFFF'
cmap=matplotlib.cm.get_cmap('jet')

if '004j' in rodada[0]:
   min_lat = -33.74+2
   max_lat = -12.32-2
   min_lon = -53.33
   max_lon = -32.46-2

   lat_step_plot = 4
   lon_step_plot = 4
   font_size = 10
   line_color = ['red','blue','green','yellow']
   resolution_plot = 300

elif '008d' in rodada[0]:
   min_lat = -77.54
   max_lat = 48.27
   min_lon = -98.00
   max_lon = 43.00

   lat_step_plot = 10
   lon_step_plot = 10
   font_size = 8
   line_color = ['red','blue','green','yellow']
   resolution_plot = 300

elif '008i' in rodada[0]:
   min_lat = -47.93
   max_lat = 10.11
   min_lon = -70.00
   max_lon = -17.75

   lat_step_plot = 5
   lon_step_plot = 5
   font_size = 8
   line_color = ['red','blue','green','yellow']
   resolution_plot = 300

if min_lat<0:
   s1 = "S"
else:
   s1 = "N"
if max_lat<0:
   s2 = "S"
else:
   s2 = "N"
if min_lon<0:
   s3 = "W"
else:
   s3 = "E"
if max_lon<0:
   s4 = "W"
else:
   s4 = "E"
dominio_str = f"{abs(min_lat):.2f}"+s1+"-"+f"{abs(max_lat):.2f}"+s2+"_" + f"{abs(min_lon):.2f}"+s3+"-"+f"{abs(max_lon):.2f}"+s4

for rod in range(0,len(rodada),1):
    filename = "SST_OSTIA_GRID_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str+".mat"
    mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[rod]+"/"+filename)
    if rod==0:
       locals()['lat_ostia'] = mat_contents['lat_ostia'].flatten()
       locals()['lon_ostia'] = mat_contents['lon_ostia'].flatten()
    locals()['sst_hycom_mean_'+rodada[rod]] = mat_contents['sst_hycom_mean_'+rodada[rod]]
    locals()['sst_ostia_mean'] = mat_contents['sst_ostia_mean']
    locals()['sst_hycom_std_'+rodada[rod]] = mat_contents['sst_hycom_std_'+rodada[rod]]
    locals()['sst_ostia_std'] = mat_contents['sst_ostia_std']
    locals()['sst_rmsd_map_'+rodada[rod]] = mat_contents['sst_rmsd_map_'+rodada[rod]]
    locals()['sst_rmsd_time_'+rodada[rod]] = mat_contents['sst_rmsd_time_'+rodada[rod]]

    cond = abs(mat_contents['sst_hycom_mean_'+rodada[rod]]) > 99999
    locals()['sst_hycom_mean_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['sst_hycom_mean_'+rodada[rod]])
    cond = abs(mat_contents['sst_ostia_mean']) > 99999
    locals()['sst_ostia_mean'] = np.ma.masked_where(cond,mat_contents['sst_ostia_mean'])

    cond = abs(mat_contents['sst_hycom_std_'+rodada[rod]]) > 99999
    locals()['sst_hycom_std_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['sst_hycom_std_'+rodada[rod]])
    cond = abs(mat_contents['sst_ostia_std']) > 99999
    locals()['sst_ostia_std'] = np.ma.masked_where(cond,mat_contents['sst_ostia_std'])

    cond = abs(mat_contents['sst_rmsd_map_'+rodada[rod]]) > 99999
    locals()['sst_rmsd_map_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['sst_rmsd_map_'+rodada[rod]])

    cond = abs(mat_contents['sst_rmsd_time_'+rodada[rod]]) > 99999
    locals()['sst_rmsd_time_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['sst_rmsd_time_'+rodada[rod]]).flatten()
    del mat_contents, cond

####### SST MEAN MAP ###########

sst_ostia_mean = np.ma.masked_where(np.ma.getmask(locals()['sst_hycom_mean_'+rodada[rod]]),sst_ostia_mean)
sst_ostia_std = np.ma.masked_where(np.ma.getmask(locals()['sst_hycom_std_'+rodada[rod]]),sst_ostia_std)

if len(rodada)<=2:
   fig = plt.figure(figsize=(6.4, 4.8), constrained_layout=True)
   gs = gridspec.GridSpec(1, len(rodada)+1)
   ax = fig.add_subplot(gs[0, 0:1],projection=ccrs.PlateCarree())
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==4:
   fig = plt.figure(figsize=(6.4, 4.8), constrained_layout=True)
   gs = gridspec.GridSpec(3, 4)
   ax = fig.add_subplot(gs[0, 1:3],projection=ccrs.PlateCarree())

all_rod_save = rodada[0]
for rod in range(1,len(rodada),1):
    all_rod_save = all_rod_save+"_"+rodada[rod]
out_dir = fig_output_dir+"/"+all_rod_save+"/SST"
if not (os.path.isdir(out_dir)):
   os.system("mkdir -p "+out_dir)

ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon)+2,int(max_lon)-1,lon_step_plot)

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax.contourf(lon_ostia, lat_ostia, locals()['sst_ostia_mean'], np.arange(0,31,1),cmap=cmap,vmin=0,vmax=30,linestyles='dashed', \
                     transform=ccrs.PlateCarree(),extend='both')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
if len(rodada)>=3:
   gl.bottom_labels = False

gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}

if '004j' in rodada[0]:
   ax.text(min_lon+1, max_lat-3, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())
elif '008d' in rodada[0]:
   ax.text(-96, -10, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())
elif '008i' in rodada[0]:
   ax.text(-68, -7, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())

all_rod = 'OSTIA'
for rod in range(0,len(rodada),1):

    ####### PLOT MEAN MAP ###########

    if len(rodada)<=2:
       if rod==0:
          ax = fig.add_subplot(gs[0, 1:2],projection=ccrs.PlateCarree())
       elif rod==1:
          ax = fig.add_subplot(gs[0, 2:3],projection=ccrs.PlateCarree())
    elif len(rodada)==4:
       if rod==0:
          ax = fig.add_subplot(gs[1, 0:2],projection=ccrs.PlateCarree())
       elif rod==1:
          ax = fig.add_subplot(gs[1, 2:4],projection=ccrs.PlateCarree())
       elif rod==2:
          ax = fig.add_subplot(gs[2, 0:2],projection=ccrs.PlateCarree())
       elif rod==3:
          ax = fig.add_subplot(gs[2, 2:4],projection=ccrs.PlateCarree())

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    im = ax.contourf(lon_ostia, lat_ostia, locals()['sst_hycom_mean_'+rodada[rod]], np.arange(0,31,1),cmap=cmap, \
                         vmin=0,vmax=30,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.right_labels = False
    if len(rodada)<=2:
       gl.left_labels = False
    elif len(rodada)==3:
       if rod==0:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==2:
          gl.left_labels = False
    elif len(rodada)==4:
       if rod==0:
          gl.bottom_labels = False
       elif rod==1:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==3:
          gl.left_labels = False          
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}

    ax.tick_params(axis='both', labelsize=font_size)
    if '004j' in rodada[0]:
       ax.text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                      transform=ccrs.PlateCarree())
    elif '008d' in rodada[rod]:
       ax.text(-96, -10, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())
    elif '008i' in rodada[rod]:
       ax.text(-68, -7, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())

    all_rod = all_rod+"_"+rodada[rod]
    if rod==len(rodada)-1:
       fig.tight_layout()
       
       ## Add a colorbar axis at the bottom of the graph
       if '004j' in rodada[rod]:
          if len(rodada)==1:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
          elif len(rodada)==2:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
          elif len(rodada)==3:
             xmax = ax[0].get_position().xmax
             xmin = ax[1].get_position().xmin
             if round(min_lat,2) == -33.74 and max_lat == -12.32 and min_lon == -53.33 and max_lon == -32.46:
                fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-4*(xmin-xmax)+0.03, hspace=0.02)
             elif round(min_lat,2) == -31.74 and max_lat == -14.32 and min_lon == -53.33 and max_lon == -34.46:
                fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-2*(xmin-xmax)-0.12, hspace=0.02)
             xmax = ax[rod+1].get_position().xmax
             ymax = ax[0].get_position().ymax
             ymin = ax[rod+1].get_position().ymin
             cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       elif '008d' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.310, 0.02, 0.37])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.715)
       elif '008i' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.985, 0.16, 0.02, 0.68])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=0.880)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.310, 0.02, 0.37])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.715)
          elif len(rodada)==4:
             fig.subplots_adjust(bottom=0.0, top=0.98, left=0.0, wspace=-0.82, hspace= 0.02)
             cbar_ax = fig.add_axes([0.72, 0.005, 0.02, 0.95])
             plt.suptitle('SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=1.02)

       cbar=fig.colorbar(im, ticks = np.arange(0,32,2), cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('\N{DEGREE SIGN}C')
       cbar.ax.tick_params(labelsize=9)

       fig.savefig(out_dir+'/SST_MEAN_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

del all_rod,fig,ax

####### DIFF SST MEAN ##############
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

if len(rodada)<=3:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==4:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
for rod in range(0,len(rodada),1):
    if len(rodada)>1:
       if rod==0: 
          ax = ax.flatten()
       ax_curr = ax[rod]
    else:
       ax_curr = ax

    ax_curr.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    im = ax_curr.contourf(lon_ostia, lat_ostia, locals()['sst_hycom_mean_'+rodada[rod]]-sst_ostia_mean, np.arange(-1.0,1.0001,.1), \
                     cmap='bwr',vmin=-1.0,vmax=1.0001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax_curr.contour(lon_ostia, lat_ostia, locals()['sst_hycom_mean_'+rodada[rod]]-sst_ostia_mean,levels = [0.0],colors=('k',), \
           linestyles=('--',),linewidths=(1))
    ax_curr.add_feature(cfeature.LAND)
    ax_curr.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax_curr.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
    gl = ax_curr.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    mean = np.mean(abs(locals()['sst_hycom_mean_'+rodada[rod]]-sst_ostia_mean).compressed())
    gl.xlabels_top = False
    gl.right_labels = False
    if (len(rodada)<=3) and (rod>0):
       gl.left_labels = False
    elif (len(rodada)==4):
       if rod==0:
          gl.bottom_labels = False
       elif rod==1:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==3:
          gl.left_labels = False       
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}

    if len(rodada)<=3:
       ax_curr.tick_params(axis='both', labelsize=font_size)
       if '004j' in rodada[rod]:
          ax_curr.text(min_lon+1, max_lat-3, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008d' in rodada[rod]:
          ax_curr.text(-96, -10, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
             transform=ccrs.PlateCarree())
          ax_curr.text(-92, -14, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
             transform=ccrs.PlateCarree())
       elif '008i' in rodada[rod]:
          ax_curr.text(-68, -7, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
             transform=ccrs.PlateCarree())
          ax_curr.text(-68, -8.5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
             transform=ccrs.PlateCarree())
    else:
       ax_curr.tick_params(axis='both', labelsize=font_size)
       if '004j' in rodada[rod]:
          ax_curr.text(min_lon+1, max_lat-3, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008d' in rodada[rod]:
          ax_curr.text(-96, -10, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-92, -14, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008i' in rodada[rod]:
          ax_curr.text(-68, -7, legenda[rod]+' - OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-68, -10.0, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]

    if rod==len(rodada)-1:
       fig.tight_layout()

       ## Add a colorbar axis at the bottom of the graph
       if '004j' in rodada[rod]:
          if len(rodada)==1:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
          elif len(rodada)==3:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             xmax = ax_curr.get_position().xmax
             ymax = ax[0].get_position().ymax
             ymin = ax_curr.get_position().ymin
             cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       elif '008d' in rodada[rod]:
          if len(rodada)==1:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.205, 0.02, 0.57])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.815)
       elif '008i' in rodada[rod]:
          if len(rodada)==1:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.83, 0.012, 0.02, 0.975])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.020)
          elif len(rodada)==2:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.205, 0.02, 0.57])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.815)
          elif len(rodada)==4:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.47,hspace= 0.02)
             cbar_ax = fig.add_axes([0.835, 0.005, 0.02, 0.985])
             plt.suptitle('DIFF SST MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=1.025)

       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('\N{DEGREE SIGN}C')

       fig.savefig(out_dir+'/DIFF_MAP_'+all_rod+'-OSTIA_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()

######### DESVIO PADRAO ###########
if len(rodada)<=2:
   fig = plt.figure(figsize=(6.4, 4.8), constrained_layout=True)
   gs = gridspec.GridSpec(1, len(rodada)+1)
   ax = fig.add_subplot(gs[0, 0:1],projection=ccrs.PlateCarree())
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==4:
   fig = plt.figure(figsize=(6.4, 4.8), constrained_layout=True)
   gs = gridspec.GridSpec(3, 4)

ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon)+2,int(max_lon)-1,lon_step_plot)

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax.contourf(lon_ostia, lat_ostia, locals()['sst_ostia_std'], np.arange(.0,3.5+.25,.25),cmap=cmap,vmin=0.0,vmax=3.5,linestyles='dashed', \
                     transform=ccrs.PlateCarree(),extend='both')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
if len(rodada)>=3:
   gl.bottom_labels = False
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}

if '004j' in rodada[0]:
   ax.text(min_lon+1, max_lat-3, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())
elif '008d' in rodada[0]:
   ax.text(-96, -10, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
           transform=ccrs.PlateCarree())
elif '008i' in rodada[0]:
   ax.text(-68, -7, 'OSTIA', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
           transform=ccrs.PlateCarree())

all_rod = 'OSTIA'
for rod in range(0,len(rodada),1):

    ####### PLOT STD MAP ###########

    if len(rodada)<=2:
       if rod==0:
          ax = fig.add_subplot(gs[0, 1:2],projection=ccrs.PlateCarree())
       elif rod==1:
          ax = fig.add_subplot(gs[0, 2:3],projection=ccrs.PlateCarree())
    elif len(rodada)==4:
       if rod==0:
          ax = fig.add_subplot(gs[1, 0:2],projection=ccrs.PlateCarree())
       elif rod==1:
          ax = fig.add_subplot(gs[1, 2:4],projection=ccrs.PlateCarree())
       elif rod==2:
          ax = fig.add_subplot(gs[2, 0:2],projection=ccrs.PlateCarree())
       elif rod==3:
          ax = fig.add_subplot(gs[2, 2:4],projection=ccrs.PlateCarree())

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    im = ax.contourf(lon_ostia, lat_ostia, locals()['sst_hycom_std_'+rodada[rod]], np.arange(.0,3.5+.25,.25),cmap=cmap, \
                         vmin=0.0,vmax=3.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if len(rodada)<=2:
       gl.ylabels_left = False
    elif len(rodada)==3:
       if rod==0:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==2:
          gl.left_labels = False
    elif len(rodada)==4:
       if rod==0:
          gl.bottom_labels = False
       elif rod==1:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==3:
          gl.left_labels = False          
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}

    ax.tick_params(axis='both', labelsize=font_size)
    if '004j' in rodada[rod]:
       ax.text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
               transform=ccrs.PlateCarree())
    elif '008d' in rodada[rod]:
       ax.text(-96, -10, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
               transform=ccrs.PlateCarree())
    elif '008i' in rodada[rod]:
       ax.text(-68, -7, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
               transform=ccrs.PlateCarree())

    all_rod = all_rod+"_"+rodada[rod]
    if rod==len(rodada)-1:
       fig.tight_layout()
       
       ## Add a colorbar axis at the bottom of the graph
       if '004j' in rodada[rod]:
          if len(rodada)==1:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
          elif len(rodada)==2:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
             cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
          elif len(rodada)==3:
             xmax = ax[0].get_position().xmax
             xmin = ax[1].get_position().xmin
             if round(min_lat,2) == -33.74 and max_lat == -12.32 and min_lon == -53.33 and max_lon == -32.46:
                fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-4*(xmin-xmax)+0.03, hspace=0.02)
             elif round(min_lat,2) == -31.74 and max_lat == -14.32 and min_lon == -53.33 and max_lon == -34.46:
                fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-2*(xmin-xmax)-0.12, hspace=0.02)
             xmax = ax[rod+1].get_position().xmax
             ymax = ax[0].get_position().ymax
             ymin = ax[rod+1].get_position().ymin
             cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       elif '008d' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.310, 0.02, 0.37])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.715)
       elif '008i' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.985, 0.16, 0.02, 0.68])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=0.880)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.310, 0.02, 0.37])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.715)
          elif len(rodada)==4:
             fig.subplots_adjust(bottom=0.0, top=0.98, left=0.0, wspace=-0.82, hspace= 0.02)
             cbar_ax = fig.add_axes([0.72, 0.005, 0.02, 0.95])
             plt.suptitle('SST STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=1.02)

       cbar=fig.colorbar(im, ticks = np.arange(.0,3.65,.5), cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('\N{DEGREE SIGN}C')

       fig.savefig(out_dir+'/SST_STD_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

del all_rod,fig,ax

####### RMSD MAP ###################
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

if len(rodada)<=3:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==4:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
for rod in range(0,len(rodada),1):
    if len(rodada)>1:
       if rod==0: 
          ax = ax.flatten()
       ax_curr = ax[rod]
    else:
       ax_curr = ax 

    ax_curr.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    im = ax_curr.contourf(lon_ostia, lat_ostia, locals()['sst_rmsd_map_'+rodada[rod]], np.arange(0,3.001,.125),cmap=cmap, \
                         vmin=0,vmax=3.001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax_curr.add_feature(cfeature.LAND)
    ax_curr.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax_curr.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax_curr.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl.xlabels_top = False
    gl.ylabels_right = False
    if (len(rodada)<=3) and (rod>0):
       gl.left_labels = False
    elif (len(rodada)==4):
       if rod==0:
          gl.bottom_labels = False
       elif rod==1:
          gl.left_labels = False
          gl.bottom_labels = False
       elif rod==3:
          gl.left_labels = False       
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    mean = np.mean(locals()['sst_rmsd_map_'+rodada[rod]].compressed())

    if len(rodada)<=3:
       ax_curr.tick_params(axis='both', labelsize=font_size)
       if '004j' in rodada:
          ax_curr.text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())
          ax_curr.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008d' in rodada[0]:
          ax_curr.text(-96, -10, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-92, -14, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008i' in rodada[0]:
          ax_curr.text(-68, -7, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-68, -8.5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())

    else:
       ax_curr.tick_params(axis='both', labelsize=font_size)
       if '004j' in rodada[rod]:
          ax_curr.text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008d' in rodada[rod]:
          ax_curr.text(-96, -10, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-92, -14, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
       elif '008i' in rodada[rod]:
          ax_curr.text(-68, -7, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())
          ax_curr.text(-68, -10.0, str(round(mean,2)), horizontalalignment='left', fontsize=font_size, fontweight='bold', \
                  transform=ccrs.PlateCarree())

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]

    if rod==len(rodada)-1:
       fig.tight_layout()
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)

       if '004j' in rodada[rod]:
          ## Add a colorbar axis at the bottom of the graph
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
          elif len(rodada)==3:
             xmax = ax_curr.get_position().xmax
             ymax = ax[0].get_position().ymax
             ymin = ax_curr.get_position().ymin
             cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       elif '008d' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.205, 0.02, 0.57])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.815)
       elif '008i' in rodada[rod]:
          if len(rodada)==1:
             cbar_ax = fig.add_axes([0.83, 0.013, 0.02, 0.985])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.025)
          elif len(rodada)==2:
             cbar_ax = fig.add_axes([0.99, 0.205, 0.02, 0.57])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.815)
          elif len(rodada)==4:
             fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.47,hspace= 0.02)
             cbar_ax = fig.add_axes([0.835, 0.005, 0.02, 0.985])
             plt.suptitle('SST RMSD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=1.025)             

       cbar=fig.colorbar(im, ticks = np.arange(0,3.001,.25), cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('\N{DEGREE SIGN}C')

       fig.savefig(out_dir+'/RMSD_MAP_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
del ax

####### PLOT RMSD TIME ###########
plt.figure(figsize=(6.4, 4.8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    data_plot = pd.Series(data=locals()['sst_rmsd_time_'+rodada[rod]], index = pd.date_range(data_inicial,data_final,freq='d'))
    ax = data_plot.plot(label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0,1.51,.2)) 
       ax.set_yticks(np.arange(0,1.51,.1),minor=True) 
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.1f'))
       plt.ylabel('\N{DEGREE SIGN}C')
       plt.ylim(0,1.51)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('SST RMSD')
       plt.savefig(out_dir+"/SST_HYCOM_RMSD_TIME_"+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                  '_'+dominio_str+"_"+all_rod+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
    del date_vec, data_plot
exit()

