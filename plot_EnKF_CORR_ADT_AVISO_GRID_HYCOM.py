import datetime
import remo
import time
import numpy as np
import scipy.io
import os
import sys
import ast
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
expt = ast.literal_eval(content_list[0])
rodada = ast.literal_eval(content_list[1])
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

input_dir = '/home/filipe.costa/resultados'
fig_output_dir = '/home/filipe.costa/resultados/figs'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

under='#000033'
over='#330000'
bad='#FFFFFF'
cmap=matplotlib.cm.get_cmap('jet')

if '004j' in rodada[0][0]:
   min_lat = -33.74
   max_lat = -12.32
   min_lon = -53.33
   max_lon = -32.46

   lat_step_plot = 4
   lon_step_plot = 4
   line_color = ['red','blue']
   resolution_plot = 300

elif '008d' in rodada[0][0]:
   min_lat = -79.54
   max_lat = 50.27
   min_lon = -98.00
   max_lon = 45.00

   lon_leg = -95
   lat_leg = -5
   lat_leg_mean = -10

   min_mean = -1
   max_mean = 1.2
   inc_mean = 0.05
   min_diff = -0.4
   max_diff = 0.4
   inc_diff = 0.02
   min_std = 0
   max_std = 0.4
   inc_std = 0.01
   lat_step_plot = 15
   lon_step_plot = 15
   line_color = ['red','blue']
   resolution_plot = 300

elif '008i' in rodada[0][0]:
   min_lat = -47.93
   max_lat = 10.11
   min_lon = -70.00
   max_lon = -17.75

   lon_leg = -68
   lat_leg = -5
   lat_leg_mean = -7.5

   min_mean = -.5
   max_mean = 1.0
   inc_mean = 0.05
   min_diff = -0.5
   max_diff = 0.5
   inc_diff = 0.025
   min_std = 0
   max_std = 0.4
   inc_std = 0.01
   lat_step_plot = 15
   lon_step_plot = 15
   line_color = ['red','blue']
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
EnKF_runs = 0

for ii in range(0,len(rodada),1):
    all_rod = rodada[ii][0]
    for rod in range(0,len(rodada[ii]),1):
        if len(rodada[ii])>0 and rod>0:
           all_rod = all_rod + '_' + rodada[ii][rod][5:]

        if len(rodada[ii])>0 and rod==len(rodada[ii])-1:
           EnKF_runs = EnKF_runs+1
           filename = "CORR_ADT_AVISO_GRID_HYCOM_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+ \
                      dominio_str+".mat"
           mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[ii][0]+"/"+filename)

           locals()['ssh_hycom_mean_all_'+rodada[ii][0]] = mat_contents['ssh_hycom_mean_'+rodada[ii][0]]
           locals()['ssh_hycom_std_all_'+rodada[ii][0]] = mat_contents['ssh_hycom_std_'+rodada[ii][0]]
           locals()['ssh_corr_map_all_'+rodada[ii][0]] = mat_contents['ssh_corr_map_'+rodada[ii][0]]
           locals()['ssh_rmsd_map_all_'+rodada[ii][0]] = mat_contents['ssh_rmsd_map_'+rodada[ii][0]]
           locals()['ssh_corr_time_all_'+rodada[ii][0]] = mat_contents['ssh_corr_time_'+rodada[ii][0]]
           locals()['ssh_rmsd_time_all_'+rodada[ii][0]] = mat_contents['ssh_rmsd_time_'+rodada[ii][0]]
    
           cond = abs(mat_contents['ssh_hycom_mean_'+rodada[ii][0]]) > 99999
           locals()['ssh_hycom_mean_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_mean_'+rodada[ii][0]])

           cond = abs(mat_contents['ssh_hycom_std_'+rodada[ii][0]]) > 99999
           locals()['ssh_hycom_std_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_std_'+rodada[ii][0]])

           cond = abs(mat_contents['ssh_corr_map_'+rodada[ii][0]]) > 99999
           locals()['ssh_corr_map_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_corr_map_'+rodada[ii][0]])
           cond = abs(mat_contents['ssh_rmsd_map_'+rodada[ii][0]]) > 99999
           locals()['ssh_rmsd_map_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_map_'+rodada[ii][0]])

           cond = abs(mat_contents['ssh_corr_time_'+rodada[ii][0]]) > 99999
           locals()['ssh_corr_time_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_corr_time_'+rodada[ii][0]]).flatten()
           cond = abs(mat_contents['ssh_rmsd_time_'+rodada[ii][0]]) > 99999
           locals()['ssh_rmsd_time_all_'+rodada[ii][0]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_time_'+rodada[ii][0]])
           del mat_contents, cond

        filename = "CORR_ADT_AVISO_GRID_HYCOM_"+rodada[ii][rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+ \
                   dominio_str+".mat"
        mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[ii][rod]+"/"+filename)
        if ii==0 and rod==0:
           locals()['lat_aviso'] = mat_contents['lat_aviso'].flatten()
           locals()['lon_aviso'] = mat_contents['lon_aviso'].flatten()
           locals()['ssh_aviso_mean'] = mat_contents['ssh_aviso_mean']
           locals()['ssh_aviso_std'] = mat_contents['ssh_aviso_std']

           cond = abs(mat_contents['ssh_aviso_mean']) > 99999
           locals()['ssh_aviso_mean'] = np.ma.masked_where(cond,mat_contents['ssh_aviso_mean'])
           cond = abs(mat_contents['ssh_aviso_std']) > 99999
           locals()['ssh_aviso_std'] = np.ma.masked_where(cond,mat_contents['ssh_aviso_std'])

        locals()['ssh_hycom_mean_'+rodada[ii][rod]] = mat_contents['ssh_hycom_mean_'+rodada[ii][rod]]
        locals()['ssh_hycom_std_'+rodada[ii][rod]] = mat_contents['ssh_hycom_std_'+rodada[ii][rod]]
        locals()['ssh_corr_map_'+rodada[ii][rod]] = mat_contents['ssh_corr_map_'+rodada[ii][rod]]
        locals()['ssh_rmsd_map_'+rodada[ii][rod]] = mat_contents['ssh_rmsd_map_'+rodada[ii][rod]]
        locals()['ssh_corr_time_'+rodada[ii][rod]] = mat_contents['ssh_corr_time_'+rodada[ii][rod]]
        locals()['ssh_rmsd_time_'+rodada[ii][rod]] = mat_contents['ssh_rmsd_time_'+rodada[ii][rod]]
    
        cond = abs(mat_contents['ssh_hycom_mean_'+rodada[ii][rod]]) > 99999
        locals()['ssh_hycom_mean_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_mean_'+rodada[ii][rod]])
    
        cond = abs(mat_contents['ssh_hycom_std_'+rodada[ii][rod]]) > 99999
        locals()['ssh_hycom_std_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_hycom_std_'+rodada[ii][rod]])
    
        cond = abs(mat_contents['ssh_corr_map_'+rodada[ii][rod]]) > 99999
        locals()['ssh_corr_map_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_corr_map_'+rodada[ii][rod]])
        cond = abs(mat_contents['ssh_rmsd_map_'+rodada[ii][rod]]) > 99999
        locals()['ssh_rmsd_map_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_map_'+rodada[ii][rod]])
    
        cond = abs(mat_contents['ssh_corr_time_'+rodada[ii][rod]]) > 99999
        locals()['ssh_corr_time_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_corr_time_'+rodada[ii][rod]]).flatten()
        cond = abs(mat_contents['ssh_rmsd_time_'+rodada[ii][rod]]) > 99999
        locals()['ssh_rmsd_time_'+rodada[ii][rod]] = np.ma.masked_where(cond,mat_contents['ssh_rmsd_time_'+rodada[ii][rod]])
        del mat_contents, cond

######## ADT MEAN ###############
if len(rodada)<=2:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada)+1,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))

ax = ax.flatten()
all_rod_save = rodada[0][0]
for ii in range(1,len(rodada),1):
    all_rod_save = all_rod_save+"_"+rodada[ii][0]
out_dir = fig_output_dir+"/"+all_rod_save+"/ADT"
if not (os.path.isdir(out_dir)):
   os.system("mkdir -p "+out_dir)

ax[0].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
meridians = np.arange(int(min_lon),int(max_lon)+0.1,lon_step_plot)

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax[0].contourf(lon_aviso, lat_aviso, locals()['ssh_aviso_mean'], np.arange(min_mean,max_mean+inc_mean,inc_mean),cmap=cmap, \
           vmin=min_mean,vmax=max_mean,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
ax[0].add_feature(cfeature.LAND)
ax[0].add_feature(cfeature.COASTLINE, linewidth=0.5)

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

ax[0].text(lon_leg, lat_leg, 'AVISO', horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

all_rod = 'AVISO'
for ii in range(0,len(rodada),1):

    ####### PLOT MEAN MAP ###########

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax[ii+1].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    if len(rodada[ii]) == 1:
       offset = remo.nanmean(locals()['ssh_hycom_mean_'+rodada[ii][0]].compressed()) - remo.nanmean(locals()['ssh_aviso_mean'].compressed())
       im = ax[ii+1].contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_'+rodada[ii][0]]-offset, \
            np.arange(min_mean,max_mean+inc_mean,inc_mean),cmap=cmap,vmin=min_mean,vmax=max_mean,linestyles='dashed', \
            transform=ccrs.PlateCarree(),extend='both')
    else:
       offset = remo.nanmean(locals()['ssh_hycom_mean_all_'+rodada[ii][0]].compressed()) - remo.nanmean(locals()['ssh_aviso_mean'].compressed())
       im = ax[ii+1].contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_all_'+rodada[ii][0]]-offset, \
            np.arange(min_mean,max_mean+inc_mean,inc_mean),cmap=cmap,vmin=min_mean,vmax=max_mean,linestyles='dashed', \
            transform=ccrs.PlateCarree(),extend='both')
    ax[ii+1].add_feature(cfeature.LAND)
    ax[ii+1].add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax[ii+1].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax[ii+1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
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

    ax[ii+1].tick_params(axis='both', labelsize=7)
    ax[ii+1].text(lon_leg, lat_leg, legenda[ii], horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

    if len(rodada[ii]) == 1:
       all_rod = all_rod+"_"+rodada[ii][0]
    else:
       all_rod = all_rod+"_"+rodada[ii][0]+'AllMembers'
    if ii==len(rodada)-1:
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
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
          plt.suptitle('ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
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
ax = ax.flatten()
for ii in range(0,len(rodada),1):
    if len(rodada)==1:
       new_ax = ax
    else:
       new_ax = ax[ii]

    new_ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    if len(rodada[ii]) == 1:
       offset = remo.nanmean(locals()['ssh_hycom_mean_'+rodada[ii][0]].compressed()) - remo.nanmean(locals()['ssh_aviso_mean'].compressed())
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_'+rodada[ii][0]]-offset-ssh_aviso_mean, \
                           np.arange(min_diff,max_diff+inc_diff,inc_diff), \
                           cmap='seismic',vmin=min_diff,vmax=max_diff,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_'+rodada[ii][0]]-offset-ssh_aviso_mean,levels = [0.0],colors=('k',), \
                     linestyles=('--',),linewidths=(1))
       mean = np.mean(abs(locals()['ssh_hycom_mean_'+rodada[ii][0]]-offset-ssh_aviso_mean).compressed())
    else:
       offset = remo.nanmean(locals()['ssh_hycom_mean_all_'+rodada[ii][0]].compressed()) - remo.nanmean(locals()['ssh_aviso_mean'].compressed())
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_all_'+rodada[ii][0]]-offset-ssh_aviso_mean, \
                           np.arange(min_diff,max_diff+inc_diff,inc_diff), \
                           cmap='seismic',vmin=min_diff,vmax=max_diff,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_hycom_mean_all_'+rodada[ii][0]]-offset-ssh_aviso_mean,levels = [0.0],colors=('k',), \
                     linestyles=('--',),linewidths=(1))
       mean = np.mean(abs(locals()['ssh_hycom_mean_all_'+rodada[ii][0]]-offset-ssh_aviso_mean).compressed())

    new_ax.add_feature(cfeature.LAND)
    new_ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    new_ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
    gl = new_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl.xlabels_top = False
    gl.right_labels = False
    if (len(rodada)<=3) and (rod>0):
       gl.left_labels = False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    new_ax.tick_params(axis='both', labelsize=7)
    new_ax.text(lon_leg, lat_leg, legenda[ii]+' - AVISO', horizontalalignment='left', fontsize=7, fontweight='bold', \
                transform=ccrs.PlateCarree())
    new_ax.text(lon_leg, lat_leg_mean, str(round(mean,2)), horizontalalignment='left', fontsize=7, fontweight='bold', \
               transform=ccrs.PlateCarree())

    if ii==0 and len(rodada[ii]) == 1:
       all_rod = rodada[ii][0] 
    elif ii==0 and len(rodada[ii]) > 1:
       all_rod = rodada[ii][0]+'AllMembers'
    elif len(rodada[ii]) == 1:
       all_rod = all_rod+"_"+rodada[ii][0]
    else:
       all_rod = all_rod+"_"+rodada[ii][0]+'AllMembers'

    if ii==len(rodada)-1:
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.00)
       fig.tight_layout()

       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
          plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
       elif len(rodada)==2:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       #elif len(rodada)==3:
       #   xmax = ax[rod].get_position().xmax
       #   ymax = ax[0].get_position().ymax
       #   ymin = ax[rod].get_position().ymin
       #   cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
       #   plt.suptitle('DIFF ADT MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=ymax+0.035)

       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
       cbar.ax.set_ylabel('m')
       fig.savefig(out_dir+'/DIFF_MAP_'+all_rod+'-AVISO_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
del all_rod,fig,ax,offset,new_ax

######### DESVIO PADRAO ###########
if len(rodada)<=2:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada)+1,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
elif len(rodada)==3:
   fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))

ax = ax.flatten()
ax[0].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

im = ax[0].contourf(lon_aviso, lat_aviso, locals()['ssh_aviso_std'], np.arange(min_std,max_std+inc_std,inc_std),cmap=cmap,vmin=min_std, \
                    vmax=max_std,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
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

ax[0].text(lon_leg, lat_leg, 'AVISO', horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

all_rod = 'AVISO'
for ii in range(0,len(rodada),1):

    ####### PLOT STD MAP ###########

    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    ax[ii+1].set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    if len(rodada[ii]) == 1:
       im = ax[ii+1].contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_std_'+rodada[ii][0]], np.arange(min_std,max_std+inc_std,inc_std),cmap=cmap, \
                             vmin=min_std,vmax=max_std,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    else:
       im = ax[ii+1].contourf(lon_aviso, lat_aviso, locals()['ssh_hycom_std_all_'+rodada[ii][0]], np.arange(min_std,max_std+inc_std,inc_std),cmap=cmap, \
                             vmin=min_std,vmax=max_std,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax[ii+1].add_feature(cfeature.LAND)
    ax[ii+1].add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax[ii+1].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = ax[ii+1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
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

    ax[ii+1].tick_params(axis='both', labelsize=7)
    ax[ii+1].text(lon_leg, lat_leg, legenda[ii], horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

    if len(rodada[ii]) == 1:
       all_rod = all_rod+"_"+rodada[ii][0]
    else:
       all_rod = all_rod+"_"+rodada[ii][0]+'AllMembers'
    if ii==len(rodada)-1:
       fig.tight_layout()
       
       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
          plt.suptitle('ADT STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
       elif len(rodada)==2:
          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
          cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
          plt.suptitle('ADT STD ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
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
ax = ax.flatten()

for ii in range(0,len(rodada),1):
    if len(rodada)==1:
       new_ax = ax
    else:
       new_ax = ax[ii]
    new_ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    if len(rodada[ii])==1:
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_corr_map_'+rodada[ii][0]], np.arange(-1,1,.025),cmap='seismic', \
                            vmin=-1,vmax=1,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_corr_map_'+rodada[ii][0]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       mean = np.mean(locals()['ssh_corr_map_'+rodada[ii][0]].compressed())
    else:
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_corr_map_all_'+rodada[ii][0]], np.arange(-1,1,.025),cmap='seismic', \
                            vmin=-1,vmax=1,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_corr_map_all_'+rodada[ii][0]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       mean = np.mean(locals()['ssh_corr_map_all_'+rodada[ii][0]].compressed())
    new_ax.add_feature(cfeature.LAND)
    new_ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    new_ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = new_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if (len(rodada)<=3) and (ii>0):
       gl.left_labels = False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    new_ax.tick_params(axis='both', labelsize=7)
    new_ax.text(lon_leg, lat_leg, legenda[ii], horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())
    new_ax.text(lon_leg, lat_leg_mean, str(round(mean,2)), horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

    if ii==0 and len(rodada[ii]) == 1:
       all_rod = rodada[ii][0]
    elif ii==0 and len(rodada[ii]) > 1:
       all_rod = rodada[ii][0]+'AllMembers'
    elif len(rodada[ii]) == 1:
       all_rod = all_rod+"_"+rodada[ii][0]
    else:
       all_rod = all_rod+"_"+rodada[ii][0]+'AllMembers'

    if ii==len(rodada)-1:
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.00)
       fig.tight_layout()

       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
       elif len(rodada)==2:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)

       cbar=fig.colorbar(im, ticks = np.arange(-1,1.001,.2), cax=cbar_ax, orientation='vertical')
       fig.savefig(out_dir+'/CORR_MAP_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

       plt.close()
del ax, new_ax, all_rod

####### RMSD MAP ###################
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(bad)

if len(rodada)+EnKF_runs<=3:
   fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
ax = ax.flatten()

plot_mean_rmsd = False
for ii in range(0,len(rodada),1):
    if len(rodada)==1:
       new_ax = ax
    else:
       new_ax = ax[ii]
    new_ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    if len(rodada[ii])==1:
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_corr_map_'+rodada[ii][0]], np.arange(-1,1,.025),cmap='seismic', \
                            vmin=-1,vmax=1,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_corr_map_'+rodada[ii][0]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       mean = np.mean(locals()['ssh_corr_map_'+rodada[ii][0]].compressed())
    elif not (plot_mean_rmsd):
       im = new_ax.contourf(lon_aviso, lat_aviso, locals()['ssh_corr_map_all_'+rodada[ii][0]], np.arange(-1,1,.025),cmap='seismic', \
                            vmin=-1,vmax=1,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
       new_ax.contour(lon_aviso, lat_aviso, locals()['ssh_corr_map_all_'+rodada[ii][0]],levels = [0.7],colors=('k',),linestyles=('--',),linewidths=(1))
       mean = np.mean(locals()['ssh_corr_map_all_'+rodada[ii][0]].compressed())
       plot_mean_rmsd = True
    else:
       plot_mean_rmsd = False
    new_ax.add_feature(cfeature.LAND)
    new_ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    new_ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

    gl = new_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if (len(rodada)<=3) and (ii>0):
       gl.left_labels = False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

    new_ax.tick_params(axis='both', labelsize=7)
    new_ax.text(lon_leg, lat_leg, legenda[ii], horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())
    new_ax.text(lon_leg, lat_leg_mean, str(round(mean,2)), horizontalalignment='left', fontsize=7, fontweight='bold', transform=ccrs.PlateCarree())

    if ii==0 and len(rodada[ii]) == 1:
       all_rod = rodada[ii][0]
    elif ii==0 and len(rodada[ii]) > 1:
       all_rod = rodada[ii][0]+'AllMembers'
    elif len(rodada[ii]) == 1:
       all_rod = all_rod+"_"+rodada[ii][0]
    else:
       all_rod = all_rod+"_"+rodada[ii][0]+'AllMembers'

    if plot_mean_rmsd:
       ii = ii-1
    if ii==len(rodada)-1:
       fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.00)
       fig.tight_layout()

       ## Add a colorbar axis at the bottom of the graph
       if len(rodada)==1:
          cbar_ax = fig.add_axes([0.86, 0.02, 0.02, 0.93])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')',y=1.000)
       elif len(rodada)==2:
          cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
          plt.suptitle('ADT CORR ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)

       cbar=fig.colorbar(im, ticks = np.arange(-1,1.001,.2), cax=cbar_ax, orientation='vertical')
       fig.savefig(out_dir+'/CORR_MAP_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+dominio_str+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)

       plt.close()
del ax, new_ax, all_rod

exit()
####### PLOT CORR TIME ###########
plt.figure(figsize=(6.4, 4.8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):

    if rod==0:
       all_rod = rodada[ii][rod]
    else:
       all_rod = all_rod+"_"+rodada[ii][rod]
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    data_plot = pd.Series(data=locals()['ssh_corr_time_'+rodada[ii][rod]], index = pd.date_range(data_inicial,data_final,freq='d'))
    ax = data_plot.plot(label=rodada[ii][rod], color=line_color[rod])
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0,1.01,.2)) 
       ax.set_yticks(np.arange(0,1.01,.1),minor=True) 
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       #ax.grid(False, which='minor', axis = 'y')
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=2))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.1f'))
       plt.ylim(0,1)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.savefig(out_dir+"/ADT_HYCOM_CORR_TIME_"+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                  '_'+dominio_str+"_"+all_rod+".png",dpi=resolution_plot)
    del date_vec, data_plot
exit()

