import datetime
import remo
import time
import numpy as np
import scipy.io
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
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

month_label_interval = 1
input_dir = '/home/filipe/resultados'
fig_output_dir = '/home/filipe/resultados/figs'

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
   line_color = ['red','blue','green']
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
    filename = "CORR_ADT_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+ \
               dominio_str+".mat"
    mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[rod]+"/"+filename)
    if rod==0:
       locals()['lat_roms'] = mat_contents['lat_roms'].flatten()
       locals()['lon_roms'] = mat_contents['lon_roms'].flatten()
    locals()['ssh_hycom_mean_'+rodada[rod]] = mat_contents['ssh_hycom_mean_'+rodada[rod]]
    locals()['ssh_roms_mean'] = mat_contents['ssh_roms_mean']
    locals()['ssh_hycom_std_'+rodada[rod]] = mat_contents['ssh_hycom_std_'+rodada[rod]]
    locals()['ssh_roms_std'] = mat_contents['ssh_roms_std']
    locals()['ssh_corr_map_'+rodada[rod]] = mat_contents['ssh_corr_map_'+rodada[rod]]
    locals()['ssh_rmsd_map_'+rodada[rod]] = mat_contents['ssh_rmsd_map_'+rodada[rod]]
    locals()['ssh_corr_time_'+rodada[rod]] = mat_contents['ssh_corr_time_'+rodada[rod]]
    locals()['ssh_rmsd_time_'+rodada[rod]] = mat_contents['ssh_rmsd_time_'+rodada[rod]]

    cond = abs(mat_contents['ssh_hycom_mean_'+rodada[rod]]) > 99999
    locals()['ssh_hycom_mean_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_mean_'+rodada[rod]])
    cond = abs(mat_contents['ssh_roms_mean']) > 99999
    locals()['ssh_roms_mean'] = np.ma.masked_where(cond,mat_contents['ssh_roms_mean'])

    cond = abs(mat_contents['ssh_hycom_std_'+rodada[rod]]) > 99999
    locals()['ssh_hycom_std_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_std_'+rodada[rod]])
    cond = abs(mat_contents['ssh_roms_std']) > 99999
    locals()['ssh_roms_std'] = np.ma.masked_where(cond,mat_contents['ssh_roms_std'])

    cond = abs(mat_contents['ssh_corr_map_'+rodada[rod]]) > 99999
    locals()['ssh_corr_map_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_corr_map_'+rodada[rod]])
    cond = abs(mat_contents['ssh_rmsd_map_'+rodada[rod]]) > 99999
    locals()['ssh_rmsd_map_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_map_'+rodada[rod]])

    cond = abs(mat_contents['ssh_corr_time_'+rodada[rod]]) > 99999
    locals()['ssh_corr_time_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_corr_time_'+rodada[rod]]).flatten()
    cond = abs(mat_contents['ssh_rmsd_time_'+rodada[rod]]) > 99999
    locals()['ssh_rmsd_time_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_time_'+rodada[rod]])
    del mat_contents, cond

######## ADT MEAN ###############
if len(rodada)<=2:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada)+1,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))

ax = ax.flatten()
all_rod_save = rodada[0]
for rod in range(1,len(rodada),1):
    all_rod_save = all_rod_save+"_"+rodada[rod]
out_dir = fig_output_dir+"/"+all_rod_save+"/ADT"
if not (os.path.isdir(out_dir)):
   os.system("mkdir -p "+out_dir)

ax[0].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon)+2,int(max_lon)-1,lon_step_plot)

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax[0].contourf(lon_roms, lat_roms, locals()['ssh_roms_mean'], np.arange(-0.10,.4501,.0200),cmap=cmap,vmin=-0.10,vmax=0.4501, \
        linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
#ax[0].add_feature(cfeature.GSHHSFeature(scale='high'))
ax[0].add_feature(cfeature.LAND)
ax[0].add_feature(cfeature.COASTLINE)

states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
ax[0].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
if len(rodada)==3:
   gl.bottom_labels = False

gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

ax[0].text(min_lon+1, max_lat-3, 'Nature Run', horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())

all_rod = 'ROMS'
for rod in range(0,len(rodada),1):

    ####### PLOT MEAN MAP ###########

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax[rod+1].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    offset = remo.nanmean(locals()['ssh_hycom_mean_'+rodada[rod]].compressed()) - remo.nanmean(locals()['ssh_roms_mean'].compressed())
    im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['ssh_hycom_mean_'+rodada[rod]]-offset, np.arange(-0.10,.4501,.0200),cmap=cmap, \
                         vmin=-0.10,vmax=0.4501,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #ax[rod+1].add_feature(cfeature.GSHHSFeature(scale='intermediate'))
    #ax[rod+1].add_feature(cfeature.GSHHSFeature(scale='high'))
    ax[rod+1].add_feature(cfeature.LAND)
    ax[rod+1].add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax[rod+1].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax[rod+1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
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
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    ax[rod+1].tick_params(axis='both', labelsize=8)
    ax[rod+1].text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())

    all_rod = all_rod+"_"+rodada[rod]
    if rod==len(rodada)-1:
       fig.tight_layout()
       
       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       elif len(rodada)==2:
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
          plt.suptitle('ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
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
          plt.suptitle('ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('m')

       fig.savefig(out_dir+'/ADT_MEAN_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

del all_rod,fig,ax,offset

####### DIFF ADT MEAN ##############
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

if len(rodada)<=3:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
for rod in range(0,len(rodada),1):
    offset = remo.nanmean(locals()['ssh_hycom_mean_'+rodada[rod]].compressed()) - remo.nanmean(locals()['ssh_roms_mean'].compressed())
    if len(rodada)==1:
       ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
       im = ax.contourf(lon_roms, lat_roms, locals()['ssh_hycom_mean_'+rodada[rod]]-offset-ssh_roms_mean, np.arange(-0.3,0.3001,.025), \
                        cmap='bwr',vmin=-0.3,vmax=0.3001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       ax.contour(lon_roms, lat_roms, locals()['ssh_hycom_mean_'+rodada[rod]]-offset-ssh_roms_mean,levels = [0.0],colors=('k',), \
              linestyles=('--',),linewidths=(1))
       #ax.add_feature(cfeature.GSHHSFeature(scale='high'))
       ax.add_feature(cfeature.LAND)
       ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
       ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
    else:
       ax[rod].set_extent([min_lon, max_lon, min_lat, max_lat])

       im = ax[rod].contourf(lon_roms, lat_roms, locals()['ssh_hycom_mean_'+rodada[rod]]-offset-ssh_roms_mean, np.arange(-0.3,0.3001,.025), \
                        cmap='bwr',vmin=-0.3,vmax=0.3001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       ax[rod].contour(lon_roms, lat_roms, locals()['ssh_hycom_mean_'+rodada[rod]]-offset-ssh_roms_mean,levels = [0.0],colors=('k',), \
                   linestyles=('--',),linewidths=(1))
       #ax[rod].add_feature(cfeature.GSHHSFeature(scale='high'))
       ax[rod].add_feature(cfeature.LAND)
       ax[rod].add_feature(cfeature.COASTLINE)
       ax[rod].add_feature(states_provinces, edgecolor='gray')

       gl = ax[rod].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    mean = np.mean(abs(locals()['ssh_hycom_mean_'+rodada[rod]]-offset-ssh_roms_mean).compressed())
    gl.xlabels_top = False
    gl.right_labels = False
    if (len(rodada)<=3) and (rod>0):
       gl.left_labels = False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    if len(rodada)==1:
       ax.tick_params(axis='both', labelsize=8)
       ax.text(min_lon+1, max_lat-3, legenda[rod]+' - Nature Run', horizontalalignment='left', fontsize=10, fontweight='bold', \
               transform=ccrs.PlateCarree())
       ax.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=10, fontweight='bold', \
               transform=ccrs.PlateCarree())
    else:
       ax[rod].tick_params(axis='both', labelsize=8)
       ax[rod].text(min_lon+1, max_lat-3, legenda[rod]+' - Nature Run', horizontalalignment='left', fontsize=10, fontweight='bold', \
               transform=ccrs.PlateCarree())
       ax[rod].text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=10, fontweight='bold', \
               transform=ccrs.PlateCarree())

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]

    if rod==len(rodada)-1:
       fig.tight_layout()
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)

       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
          plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
       elif len(rodada)==2:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       elif len(rodada)==3:
          xmax = ax[rod].get_position().xmax
          ymax = ax[0].get_position().ymax
          ymin = ax[rod].get_position().ymin
          cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
          plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)

       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('m')

       fig.savefig(out_dir+'/DIFF_MAP_'+all_rod+'-ROMS_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()

######### DESVIO PADRAO ###########
if len(rodada)<=2:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada)+1,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))

ax = ax.flatten()
ax[0].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon)+2,int(max_lon)-1,lon_step_plot)

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax[0].contourf(lon_roms, lat_roms, locals()['ssh_roms_std'], np.arange(.0,.2,.0100),cmap=cmap,vmin=0.0,vmax=0.2,linestyles='dashed', \
                     transform=ccrs.PlateCarree(),extend='both')
#ax[0].add_feature(cfeature.GSHHSFeature(scale='high'))
ax[0].add_feature(cfeature.LAND)
ax[0].add_feature(cfeature.COASTLINE, linewidth=0.5)

states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
ax[0].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
if len(rodada)==3:
   gl.bottom_labels = False
gl.xlocator = mticker.FixedLocator(meridians)
gl.ylocator = mticker.FixedLocator(parallels)
gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

ax[0].text(min_lon+1, max_lat-3, 'Nature Run', horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())

all_rod = 'ROMS'
for rod in range(0,len(rodada),1):

    ####### PLOT STD MAP ###########

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax[rod+1].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['ssh_hycom_std_'+rodada[rod]], np.arange(.0,.2,.0100),cmap=cmap, \
                         vmin=0.0,vmax=0.2,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #ax[rod+1].add_feature(cfeature.GSHHSFeature(scale='high'))
    ax[rod+1].add_feature(cfeature.LAND)
    ax[rod+1].add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax[rod+1].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax[rod+1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
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
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    ax[rod+1].tick_params(axis='both', labelsize=8)
    ax[rod+1].text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())

    all_rod = all_rod+"_"+rodada[rod]
    if rod==len(rodada)-1:
       fig.tight_layout()
       
       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('ADT STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       elif len(rodada)==2:
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
          plt.suptitle('ADT STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
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
          plt.suptitle('ADT STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('m')

       fig.savefig(out_dir+'/ADT_STD_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

del all_rod,fig,ax

####### CORR MAP ###################
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

if len(rodada)<=3:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
for rod in range(0,len(rodada),1):
    if len(rodada)==1:
       ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

       im = ax.contourf(lon_roms, lat_roms, locals()['ssh_corr_map_'+rodada[rod]], np.arange(-1,1.001,.025),cmap='bwr', \
                            vmin=-1,vmax=1.001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       ax.contour(lon_roms, lat_roms, locals()['ssh_corr_map_'+rodada[rod]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       #ax.add_feature(cfeature.GSHHSFeature(scale='high'))
       ax.add_feature(cfeature.LAND)
       ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
       ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

       gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    else:
       ax[rod].set_extent([min_lon, max_lon, min_lat, max_lat])

       im = ax[rod].contourf(lon_roms, lat_roms, locals()['ssh_corr_map_'+rodada[rod]], np.arange(-1,1.001,.025),cmap='bwr', \
                            vmin=-1,vmax=1.001,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       ax[rod].contour(lon_roms, lat_roms, locals()['ssh_corr_map_'+rodada[rod]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       #ax[rod].add_feature(cfeature.GSHHSFeature(scale='high'))
       ax[rod].add_feature(cfeature.LAND)
       ax[rod].add_feature(cfeature.COASTLINE, linewidth=0.5)
       ax[rod].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

       gl = ax[rod].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if (len(rodada)<=3) and (rod>0):
       gl.left_labels = False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    mean = np.mean(locals()['ssh_corr_map_'+rodada[rod]].compressed())

    if len(rodada)==1:
       ax.tick_params(axis='both', labelsize=8)
       ax.text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
       ax.text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
    else:
       ax[rod].tick_params(axis='both', labelsize=8)
       ax[rod].text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
       ax[rod].text(min_lon+1, max_lat-5, str(round(mean,2)), horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]

    if rod==len(rodada)-1:
       fig.tight_layout()
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)

       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
       elif len(rodada)==2:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       elif len(rodada)==3:
          xmax = ax[rod].get_position().xmax
          ymax = ax[0].get_position().ymax
          ymin = ax[rod].get_position().ymin
          cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)

       cbar=fig.colorbar(im, ticks = np.arange(-1,1.001,.2), cax=cbar_ax,orientation='vertical')

       fig.savefig(out_dir+'/CORR_MAP_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
del ax

####### PLOT CORR TIME ###########
plt.figure(figsize=(6.4, 4.8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):

    if rod==0:
       all_rod = rodada[rod]
    else:
       all_rod = all_rod+"_"+rodada[rod]
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    data_plot = pd.Series(data=locals()['ssh_corr_time_'+rodada[rod]], index = pd.date_range(data_inicial,data_final,freq='d'))
    ax = data_plot.plot(label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0,1.01,.2)) 
       ax.set_yticks(np.arange(0,1.01,.1),minor=True) 
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       #ax.grid(False, which='minor', axis = 'y')
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.1f'))
       plt.ylim(0,1)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('ADT CORRELATION')
       plt.savefig(out_dir+"/ADT_HYCOM_CORR_TIME_"+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                  '_'+dominio_str+"_"+all_rod+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
    del date_vec, data_plot
exit()

