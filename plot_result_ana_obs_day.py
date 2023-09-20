from netCDF4 import Dataset
from scipy import interpolate
from scipy.stats import pearsonr
import scipy.io
import remo
import datetime
import sys
import os
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
from mpl_toolkits.basemap import Basemap
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
inc_assim = float(content_list[6].rstrip('\n'))

fig_output_dir = '/home/filipe.costa/resultados/figs'
dir_hycom = '/home/filipe.costa/previsao/hycom_2_2_18/proc'
dir_hycom_assim = 'restart_assim'
dir_obs_adt = '/home/filipe.costa/dados_obs/get_AVISO/adt'
dir_obs_sst = '/home/filipe.costa/dados_obs/get_OSTIA_SST/ostia'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

under='#000033'
over='#330000'
bad='#FFFFFF'
cmap=matplotlib.cm.get_cmap('jet')
states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       dir_ATL = dir_hycom+'/ATLj0.04/'
       record_sst = 193  ## RESTART
       lat_step_plot = 4
       lon_step_plot = 5
       resolution_plot = 150
    elif '008d' in rodada[rod]:
       IDM = 1717
       JDM = 2345
       dir_ATL = dir_hycom+'/ATLd0.08/'
       record_sst = 193  ## RESTART
       lat_step_plot = 15
       lon_step_plot = 15
       resolution_plot = 150
       plt_min_ssh = -1.0
       plt_max_ssh = 1.4
       plt_inc_ssh = 0.10
       plt_min_diff_ssh = -0.5
       plt_max_diff_ssh = 0.5
       plt_inc_diff_ssh = 0.025
    elif '008i' in rodada[rod]:
       IDM = 628
       JDM = 780
       dir_ATL = dir_hycom+'/ATLi0.08/'
       record_sst = 193  ## RESTART
       lat_step_plot = 5
       lon_step_plot = 5
       resolution_plot = 150
       plt_min_ssh = -1.0
       plt_max_ssh = 1.4
       plt_inc_ssh = 0.10
       plt_min_diff_ssh = -0.5
       plt_max_diff_ssh = 0.5
       plt_inc_diff_ssh = 0.025
       plt_min_sst = 10
       plt_max_sst = 30
       plt_inc_sst = 1
       plt_min_diff_sst = -5
       plt_max_diff_sst = 5
       plt_inc_diff_sst = 0.5
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1

    offset = np.ma.array(np.zeros((data_final - data_inicial).days+1), mask=True)

    while(current_data <= data_final):
        print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))
        if (current_data==data_inicial):
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
           if min_lat_hycom<0:
              s1 = "S"
           else:
              s1 = "N"
           if max_lat_hycom<0:
              s2 = "S"
           else:
              s2 = "N"
           if min_lon_hycom<0:
              s3 = "W"
           else:
              s3 = "E"
           if max_lon_hycom<0:
              s4 = "W"
           else:
              s4 = "E"
           dominio_str = f"{abs(min_lat_hycom):.2f}"+s1+"-"+f"{abs(max_lat_hycom):.2f}"+s2+"_" + \
                         f"{abs(min_lon_hycom):.2f}"+s3+"-"+f"{abs(max_lon_hycom):.2f}"+s4

        arq_obs_adt = dir_obs_adt+'/'+current_data.strftime("%Y")+'/dt_global_allsat_phy_l4_'+current_data.strftime("%Y%m%d")+'.nc'
        if not (os.path.isfile(arq_obs_adt)):
           print('ADT OBS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_obs_adt)
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        arq_obs_sst = dir_obs_sst+'/'+current_data.strftime("%Y")+'/'+current_data.strftime("%Y%m%d")+ \
                      '120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB_REP-v02.0-fv02.0.nc'
        if not (os.path.isfile(arq_obs_sst)):
           print('SST OBS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_obs_sst)
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        restart_hycom_back = dir_expt+'/data/restart_temp/restart_' \
                             +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
        #print(restart_hycom_back)
        if not (os.path.isfile(restart_hycom_back)):
           #print('HYCOM BACKGROUND FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' NOT FOUND IN RESTART_TEMP')
           restart_hycom_back = dir_expt+'/data/restart_files/restart_' \
                                +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                                +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
           print('RESTART BACK: ',restart_hycom_back)
           if not (os.path.isfile(restart_hycom_back)):
              print('HYCOM BACKGROUND FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' NOT FOUND IN RESTART_FILES')
              current_data = current_data + datetime.timedelta(days=inc_tempo)
              continue

        restart_hycom_assim = dir_expt+'/data/'+dir_hycom_assim+'/restart_' \
                              +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                              +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
        print('RESTART ASSIM: ',restart_hycom_assim)
        if not (os.path.isfile(restart_hycom_assim)):
           print('HYCOM ANALYSIS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(restart_hycom_assim)

        archv_hycom_back = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d") \
                             +'/archv.'+str((current_data).timetuple().tm_year).zfill(4)+'_' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        if not (os.path.isfile(archv_hycom_back)):
           archv_hycom_back = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-inc_assim)).strftime("%Y%m%d") \
                                +'/archv.'+str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'

        archv_hycom_assim = dir_expt+'/output/ab/'+current_data.strftime("%Y%m%d") \
                             +'/archv.'+str((current_data).timetuple().tm_year).zfill(4)+'_' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'

        print('ARCHV BACK: ',archv_hycom_back)
        print('ARCHV ASSIM: ',archv_hycom_assim)
        counter = counter+1

        nc = Dataset(arq_obs_adt, 'r')
        if (current_data==data_inicial):
           lat_aviso = nc.variables['latitude'][:]
           lon_aviso = nc.variables['longitude'][:]
           cond_lat_aviso = (lat_aviso<=max_lat_hycom) & (lat_aviso>=min_lat_hycom)
           cond_lon_aviso = (lon_aviso<=max_lon_hycom) & (lon_aviso>=min_lon_hycom)
           lat_aviso = lat_aviso[cond_lat_aviso]
           lon_aviso = lon_aviso[cond_lon_aviso]
           min_lat_aviso = min(lat_aviso)
           max_lat_aviso = max(lat_aviso)
           min_lon_aviso = min(lon_aviso)
           max_lon_aviso = max(lon_aviso)
        ssh = nc.variables['adt'][:]
        ssh = np.moveaxis(ssh,0,-1)
        cond = abs(ssh) > 99999
        ssh = np.ma.masked_where(cond,ssh)
        ssh = ssh[cond_lat_aviso,:,0]
        ssh_aviso = ssh[:,cond_lon_aviso]
        del ssh


        nc = Dataset(arq_obs_sst, 'r')
        if (current_data==data_inicial):
           lat_ostia = nc.variables['lat'][:]
           lon_ostia = nc.variables['lon'][:]
           cond_lat_ostia = (lat_ostia<=max_lat_hycom) & (lat_ostia>=min_lat_hycom)
           cond_lon_ostia = (lon_ostia<=max_lon_hycom) & (lon_ostia>=min_lon_hycom)
           lat_ostia = lat_ostia[cond_lat_ostia]
           lon_ostia = lon_ostia[cond_lon_ostia]
           min_lat_ostia = min(lat_ostia)
           max_lat_ostia = max(lat_ostia)
           min_lon_ostia = min(lon_ostia)
           max_lon_ostia = max(lon_ostia)
        sst = nc.variables['analysed_sst'][:]
        sst = np.moveaxis(sst,0,-1)-273.15
        cond = abs(sst) > 99999
        sst = np.ma.masked_where(cond,sst)
        sst = sst[cond_lat_ostia,:,0]
        sst_ostia = sst[:,cond_lon_ostia]
        del sst

        f = open(restart_hycom_back,'rb')
        f.seek(record_sst*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_ostia,lat_ostia)
        cond = abs(field) > 999999
        sst_hycom_back = np.ma.masked_where(cond,field)
        sst_hycom_back = np.ma.masked_where(np.ma.getmask(sst_ostia), sst_hycom_back)
        f.close()
        del field, cond

        f = open(restart_hycom_assim,'rb')
        f.seek(record_sst*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_ostia,lat_ostia)
        cond = abs(field) > 999999
        sst_hycom_assim = np.ma.masked_where(cond,field)
        sst_hycom_assim = np.ma.masked_where(np.ma.getmask(sst_ostia), sst_hycom_assim)
        f.close()
        del field, cond

        out_dir = fig_output_dir+"/"+rodada[rod]+"/BACK_ANA/DAILY"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)

        cmap.set_under(under)
        cmap.set_over(over)
        cmap.set_bad(bad)

        ###### SST PLOTS #########
        print('PLOT SST')

        parallels = np.arange(int(min_lat_ostia),int(max_lat_ostia)+0.1,lat_step_plot)
        meridians = np.arange(int(min_lon_ostia),int(max_lon_ostia)-1,lon_step_plot)
        
        fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        ax = ax.flatten()

        for i in range(0,6,1):
            ax[i].set_extent([min_lon_ostia, max_lon_ostia, min_lat_ostia, max_lat_ostia], crs=ccrs.PlateCarree())
            gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = mticker.FixedLocator(meridians)
            gl.ylocator = mticker.FixedLocator(parallels)
            gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

            if i==0:
               im = ax[i].contourf(lon_ostia, lat_ostia, sst_ostia, np.arange(plt_min_sst,plt_max_sst+plt_inc_sst,plt_inc_sst), \
                    cmap=cmap,vmin=plt_min_sst,vmax=plt_max_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
            
            if i==1:
               im = ax[i].contourf(lon_ostia, lat_ostia, sst_hycom_back, np.arange(plt_min_sst,plt_max_sst+plt_inc_sst,plt_inc_sst), \
                    cmap=cmap,vmin=plt_min_sst,vmax=plt_max_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
               gl.left_labels = False
            
            if i==2:
               im_mean = ax[i].contourf(lon_ostia, lat_ostia, sst_hycom_assim, np.arange(plt_min_sst,plt_max_sst+plt_inc_sst,plt_inc_sst), \
                         cmap=cmap,vmin=plt_min_sst,vmax=plt_max_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
               gl.left_labels = False
            
            if i==3:
               im_diff = ax[i].contourf(lon_ostia, lat_ostia, sst_hycom_assim-sst_hycom_back, \
                         np.arange(plt_min_diff_sst,plt_max_diff_sst+plt_inc_diff_sst,plt_inc_diff_sst),cmap='bwr', \
                         vmin=plt_min_diff_sst,vmax=plt_max_diff_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               mean = abs(sst_hycom_assim-sst_hycom_back)
               cond = mean<=0.001
               mean = np.ma.masked_where(cond,mean)
               mean = np.mean(mean.compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008i' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               del mean

            if i==4:
               im_diff = ax[i].contourf(lon_ostia, lat_ostia, sst_hycom_back-sst_ostia, \
                         np.arange(plt_min_diff_sst,plt_max_diff_sst+plt_inc_diff_sst,plt_inc_diff_sst), \
                         cmap='bwr',vmin=plt_min_diff_sst,vmax=plt_max_diff_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.left_labels = False
               mean = np.mean(abs(sst_hycom_back-sst_ostia).compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008d' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
            if i==5:
               im = ax[i].contourf(lon_ostia, lat_ostia, sst_hycom_assim-sst_ostia,  \
                    np.arange(plt_min_diff_sst,plt_max_diff_sst+plt_inc_diff_sst,plt_inc_diff_sst), \
                    cmap='bwr',vmin=plt_min_diff_sst,vmax=plt_max_diff_sst,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.left_labels = False
               mean = np.mean(abs(sst_hycom_assim-sst_ostia).compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008d' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())

            ax[i].add_feature(cfeature.LAND)
            ax[i].add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax[i].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)

        if '008i' in rodada[rod]:
           ax[0].text(-68, -7, 'OSTIA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[1].text(-68, -7, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[2].text(-68, -7, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[3].text(-68, -7, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[4].text(-68, -7, 'BACK - OSTIA', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[5].text(-68, -7, 'ANA - OSTIA', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())

        elif '008d' in rodada[rod]:
           ax[0].text(-95, -5, 'OSTIA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[1].text(-95, -5, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[2].text(-95, -5, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[3].text(-95, -5, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[4].text(-95, -5, 'BACK - OSTIA', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[5].text(-95, -5, 'ANA - OSTIA', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())

        fig.tight_layout()
        fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=0.0)

        ## Add a colorbar axis at the bottom of the graph
        ymax = ax[5].get_position().ymax
        ymin = ax[5].get_position().ymin
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('m',fontsize=7)

        ymin = ax[2].get_position().ymin
        ymax = ax[2].get_position().ymax
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('m',fontsize=7)

        plt.suptitle('SST - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
        fig.savefig(out_dir+'/SST_'+current_data.strftime("%Y%m%d")+'.png',dpi=resolution_plot,transparent=False, \
                    bbox_inches='tight',pad_inches=0.05)
        plt.close()

        f = open(archv_hycom_back,'rb')
        f.seek(1*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        field = field/9.806
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_aviso,lat_aviso)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        f.close()
        offset = remo.nanmean(field.compressed()) - remo.nanmean(ssh_aviso.compressed())
        ssh_hycom_back = field-offset
        del field

        f = open(archv_hycom_assim,'rb')
        f.seek(1*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        field = field/9.806
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_aviso,lat_aviso)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        f.close()
        offset = remo.nanmean(field.compressed()) - remo.nanmean(ssh_aviso.compressed())
        ssh_hycom_assim = field-offset
        del field
        print('PLOT ADT')

        parallels = np.arange(int(min_lat_aviso),int(max_lat_aviso)+0.1,lat_step_plot)
        meridians = np.arange(int(min_lon_aviso),int(max_lon_aviso)-1,lon_step_plot)
        
        fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        ax = ax.flatten()

        for i in range(0,6,1):
            ax[i].set_extent([min_lon_aviso, max_lon_aviso, min_lat_aviso, max_lat_aviso], crs=ccrs.PlateCarree())
            gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = mticker.FixedLocator(meridians)
            gl.ylocator = mticker.FixedLocator(parallels)
            gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}

            if i==0:
               im = ax[i].contourf(lon_aviso, lat_aviso, ssh_aviso, np.arange(plt_min_ssh,plt_max_ssh+plt_inc_ssh,plt_inc_ssh), \
                    cmap=cmap,vmin=plt_min_ssh,vmax=plt_max_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
            
            if i==1:
               im = ax[i].contourf(lon_aviso, lat_aviso, ssh_hycom_back, np.arange(plt_min_ssh,plt_max_ssh+plt_inc_ssh,plt_inc_ssh), \
                    cmap=cmap,vmin=plt_min_ssh,vmax=plt_max_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
               gl.left_labels = False
            
            if i==2:
               im_mean = ax[i].contourf(lon_aviso, lat_aviso, ssh_hycom_assim, np.arange(plt_min_ssh,plt_max_ssh+plt_inc_ssh,plt_inc_ssh), \
                         cmap=cmap,vmin=plt_min_ssh,vmax=plt_max_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.bottom_labels = False
               gl.left_labels = False
            
            if i==3:
               im_diff = ax[i].contourf(lon_aviso, lat_aviso, ssh_hycom_assim-ssh_hycom_back, \
                         np.arange(plt_min_diff_ssh,plt_max_diff_ssh+plt_inc_diff_ssh,plt_inc_diff_ssh),cmap='bwr', \
                         vmin=plt_min_diff_ssh,vmax=plt_max_diff_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               mean = abs(ssh_hycom_assim-ssh_hycom_back)
               cond = mean<=0.001
               mean = np.ma.masked_where(cond,mean)
               mean = np.mean(mean.compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008d' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())

            if i==4:
               im_diff = ax[i].contourf(lon_aviso, lat_aviso, ssh_hycom_back-ssh_aviso, \
                         np.arange(plt_min_diff_ssh,plt_max_diff_ssh+plt_inc_diff_ssh,plt_inc_diff_ssh), \
                         cmap='bwr',vmin=plt_min_diff_ssh,vmax=plt_max_diff_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.left_labels = False
               mean = np.mean(abs(ssh_hycom_back-ssh_aviso).compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008d' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
            if i==5:
               im = ax[i].contourf(lon_aviso, lat_aviso, ssh_hycom_assim-ssh_aviso,  \
                    np.arange(plt_min_diff_ssh,plt_max_diff_ssh+plt_inc_diff_ssh,plt_inc_diff_ssh), \
                    cmap='bwr',vmin=plt_min_diff_ssh,vmax=plt_max_diff_ssh,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               gl.left_labels = False
               mean = np.mean(abs(ssh_hycom_assim-ssh_aviso).compressed())
               if '008i' in rodada[rod]:
                  ax[i].text(-68, -9, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               elif '008d' in rodada[rod]:
                  ax[i].text(-95, -12, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())

            ax[i].add_feature(cfeature.LAND)
            ax[i].add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax[i].add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
               

        if '008i' in rodada[rod]:
           ax[0].text(-68, -7, 'AVISO', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[1].text(-68, -7, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[2].text(-68, -7, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[3].text(-68, -7, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[4].text(-68, -7, 'BACK - AVISO', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[5].text(-68, -7, 'ANA - AVISO', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())

        elif '008d' in rodada[rod]:
           ax[0].text(-95, -5, 'AVISO', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[1].text(-95, -5, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[2].text(-95, -5, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                      transform=ccrs.PlateCarree())
           ax[3].text(-95, -5, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[4].text(-95, -5, 'BACK - AVISO', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())
           ax[5].text(-95, -5, 'ANA - AVISO', horizontalalignment='left', fontsize=7, \
                      fontweight='bold',transform=ccrs.PlateCarree())


        fig.tight_layout()
        fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=0.0)

        ## Add a colorbar axis at the bottom of the graph
        ymax = ax[5].get_position().ymax
        ymin = ax[5].get_position().ymin
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('m',fontsize=7)

        ymin = ax[2].get_position().ymin
        ymax = ax[2].get_position().ymax
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('m',fontsize=7)

        plt.suptitle('ADT - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
        fig.savefig(out_dir+'/ADT_'+current_data.strftime("%Y%m%d")+'.png',dpi=resolution_plot,transparent=False, \
                    bbox_inches='tight',pad_inches=0.05)
        plt.close()

        current_data = current_data + datetime.timedelta(days=inc_assim)

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

