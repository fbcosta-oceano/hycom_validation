from netCDF4 import Dataset
from scipy import interpolate
from scipy.stats import pearsonr
import scipy.io
import remo
import datetime
import sys
import os
import time
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import pandas as pd
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import warnings
#jet = matplotlib.cm.get_cmap('jet')
cmap = copy.copy(matplotlib.cm.get_cmap('jet'))
#cmap_diff = copy.copy(matplotlib.cm.get_cmap('bwr'))

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

fig_output_dir = '/prj/rodas/filipe.costa/resultados/figs'
dir_hycom = '/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc'
dir_hycom_assim = 'restart_assim'
dir_obs_roms = '/prj/rodas/leonardo.pires/ROMS'
dir_obs_roms_daily = '/prj/rodas/filipe.costa/ROMS/mma'

#prof = np.array([[0,0], [100,100], [400,400], [800,800], [2000,2000]])
prof = np.array([[0,0], [100,100], [400,400], [800,800]])
k_inic = 0

lat_step_plot = 4
lon_step_plot = 5
resolution_plot = 150

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       dir_ATL = dir_hycom+'/ATLj0.04/'
       record_sst = 193  ## RESTART
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

        restart_hycom_back = dir_expt+'/data/restart_temp/restart_' \
                             +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
        print(restart_hycom_back)
        if not (os.path.isfile(restart_hycom_back)):
           print('HYCOM BACKGROUND FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' NOT FOUND IN RESTART_TEMP')
           restart_hycom_back = dir_expt+'/data/restart_files/restart_' \
                                +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                                +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
           print(restart_hycom_back)
           if not (os.path.isfile(restart_hycom_back)):
              print('HYCOM BACKGROUND FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' NOT FOUND IN RESTART_FILES')
              current_data = current_data + datetime.timedelta(days=inc_tempo)
              continue

        restart_hycom_assim = dir_expt+'/data/'+dir_hycom_assim+'/restart_' \
                              +str((current_data).timetuple().tm_year).zfill(4)+'d' \
                              +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'
        print(restart_hycom_assim)

        archive_hycom_back = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
                             +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        print(archive_hycom_back)
        archive_hycom_assim = dir_expt+'/output/ab/'+current_data.strftime("%Y%m%d")+'/archv.' \
                              +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                              +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        print(archive_hycom_assim)

        if not (os.path.isfile(restart_hycom_assim)):
           print('HYCOM ANALYSIS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(restart_hycom_assim)

        counter = counter+1

        if (current_data==data_inicial):
           arq_roms1 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_1de2.nc4'
           arq_roms2 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_2de2.nc4'
           if not (os.path.isfile(arq_roms1)) or not (os.path.isfile(arq_roms2)):
              print('ROMS FILE FOR YEAR: '+current_data.strftime("%Y")+' DOES NOT EXIST')
              print(arq_roms1)
              print(arq_roms2)
              exit()
           nc1 = Dataset(arq_roms1, 'r')
           nc2 = Dataset(arq_roms2, 'r')
           lat_roms = nc1.variables['lat_rho'][:,0] #564
           lon_roms = nc1.variables['lon_rho'][0,:] #502
           cond_lat_roms = (lat_roms<=max_lat_hycom) & (lat_roms>=min_lat_hycom)
           cond_lon_roms = (lon_roms<=max_lon_hycom) & (lon_roms>=min_lon_hycom)
           lat_roms = lat_roms[cond_lat_roms]
           lon_roms = lon_roms[cond_lon_roms]
           max_lat_roms = max(lat_roms)
           min_lat_roms = min(lat_roms)
           max_lon_roms = max(lon_roms)
           min_lon_roms = min(lon_roms)
           time_roms1 = nc1.variables['ocean_time'][:]
           time_roms2 = nc2.variables['ocean_time'][:]
           time_roms = np.concatenate((time_roms1,time_roms2),axis=0)
           time_roms = time_roms/(60*60*24) #AJUSTANDO PARA DIAS
           time_roms = time_roms-.5 #SUBTRAINDO 12 h PARA FICAR EM 00Z
           ref = datetime.datetime(*time.strptime('01/01/1970',"%d/%m/%Y")[0:3])
           date_roms = [ref + datetime.timedelta(days=time_roms[i]) for i in range(0,len(time_roms),1)]

           ssh1 = nc1.variables['zeta'][:] #(238, 564, 502)
           ssh2 = nc2.variables['zeta'][:]
           ssh = np.ma.concatenate((ssh1,ssh2),axis=0)
           ssh = np.moveaxis(ssh,0,-1)
           ssh = ssh[cond_lat_roms,:,:]
           ssh_roms = ssh[:,cond_lon_roms,:]
           del ssh1, ssh2, ssh

           sst1 = nc1.variables['temp'][:] #(238, 32, 564, 502)
           sst2 = nc2.variables['temp'][:]
           sst = np.ma.concatenate((sst1,sst2),axis=0)
           
           sst = np.moveaxis(sst,1,-1)
           sst = np.moveaxis(sst,0,-1)
           sst = sst[cond_lat_roms,:,-1,:]
           sst_roms = sst[:,cond_lon_roms,:]
           del sst1, sst2, sst
           nc1.close()
           nc2.close()

        f = open(restart_hycom_back,'rb')
        f.seek(record_sst*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_roms,lat_roms)
        cond = abs(field) > 999999
        sst_hycom_back = np.ma.masked_where(cond,field)
        f.close()
        del field, cond

        f = open(restart_hycom_assim,'rb')
        f.seek(record_sst*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_roms,lat_roms)
        cond = abs(field) > 999999
        sst_hycom_assim = np.ma.masked_where(cond,field)
        f.close()
        del field, cond

        if not (os.path.isfile(archive_hycom_assim)):
           print('HYCOM ARCHIVE ANALYSIS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(archive_hycom_assim)
        else:
           f = open(archive_hycom_back,'rb')
           f.seek(1*4*(IJDM+npad))
           field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
           field = np.reshape(field,(JDM,IDM))/9.806
           I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
           field = I(lon_roms,lat_roms)
           cond = abs(field) > 999999
           ssh_hycom_back = np.ma.masked_where(cond,field)
           f.close()
           del field, cond

           f = open(archive_hycom_assim,'rb')
           f.seek(1*4*(IJDM+npad))
           field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
           field = np.reshape(field,(JDM,IDM))/9.806
           I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
           field = I(lon_roms,lat_roms)
           cond = abs(field) > 999999
           ssh_hycom_assim = np.ma.masked_where(cond,field)
           f.close()
           del field, cond


        ind = np.where(date_roms == np.datetime64(datetime.datetime(current_data.year, current_data.month, current_data.day, 0, 0)))
        ind = int(ind[0])
        cond = abs(sst_roms[:,:,ind]) > 999999
        sst_hycom_back = np.ma.masked_where(cond,sst_hycom_back)
        sst_hycom_assim = np.ma.masked_where(cond,sst_hycom_assim)

        out_dir = fig_output_dir+"/"+rodada[rod]+"/BACK_ANA/DAILY"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)

        parallels = np.arange(int(min_lat_roms),int(max_lat_roms)+0.1,lat_step_plot)
        meridians = np.arange(int(min_lon_roms)+2,int(max_lon_roms)-1,lon_step_plot)
        
        under='#000033'
        over='#330000'
        bad='#FFFFFF'
        cmap.set_under(under)
        cmap.set_over(over)
        cmap.set_bad(bad)
        states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')

        ############### 2X3 SST ####################################
        fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        ax = ax.flatten()

        t_min = 16
        t_max = 30
        delta = 0.5
        t_range = np.arange(t_min,t_max+delta,delta)
        for i in range(0,6,1):
            ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
            gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = mticker.FixedLocator(meridians)
            gl.ylocator = mticker.FixedLocator(parallels)
            gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            
            if i==0:
               im = ax[i].contourf(lon_roms, lat_roms, sst_roms[:,:,ind], t_range,cmap=cmap,vmin=t_min, \
                    vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
                          transform=ccrs.PlateCarree())
               gl.bottom_labels = False
            
            if i==1:
               im = ax[i].contourf(lon_roms, lat_roms, sst_hycom_back, t_range,cmap=cmap,vmin=t_min, \
                    vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                          transform=ccrs.PlateCarree())
               gl.bottom_labels = False
               gl.left_labels = False
            
            if i==2:
               im_mean = ax[i].contourf(lon_roms, lat_roms, sst_hycom_assim, t_range,cmap=cmap,vmin=t_min, \
                         vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                          transform=ccrs.PlateCarree())
               gl.left_labels = False
               gl.bottom_labels = False
            
            if i==3:
               im_diff = ax[i].contourf(lon_roms, lat_roms, sst_hycom_assim-sst_hycom_back, np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                         vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                          fontweight='bold',transform=ccrs.PlateCarree())
               mean = abs(sst_hycom_assim-sst_hycom_back)
               cond = mean<=0.01
               mean = np.ma.masked_where(cond,mean)
               mean = np.mean(mean.compressed())
               ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                          fontweight='bold',transform=ccrs.PlateCarree())

            if i==4:
               im_diff = ax[i].contourf(lon_roms, lat_roms, sst_hycom_back-sst_roms[:,:,ind], np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                         vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
                          fontweight='bold',transform=ccrs.PlateCarree())
               gl.left_labels = False
               mean = np.mean(abs(sst_hycom_back-sst_roms[:,:,ind]).compressed())
               ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                          fontweight='bold',transform=ccrs.PlateCarree())

            if i==5:
               im = ax[i].contourf(lon_roms, lat_roms, sst_hycom_assim-sst_roms[:,:,ind], np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                    vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
               ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
                          transform=ccrs.PlateCarree())
               gl.left_labels = False
               mean = np.mean(abs(sst_hycom_assim-sst_roms[:,:,ind]).compressed())
               ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                          fontweight='bold',transform=ccrs.PlateCarree())

            ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
            ax[i].add_feature(cfeature.LAND)
            ax[i].add_feature(cfeature.COASTLINE)

            ax[i].add_feature(states_provinces, edgecolor='gray')
            
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

        ## Add a colorbar axis at the bottom of the graph
        ymax = ax[5].get_position().ymax
        ymin = ax[5].get_position().ymin
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('\N{DEGREE SIGN}C',fontsize=8)

        ymin = ax[2].get_position().ymin
        ymax = ax[2].get_position().ymax
        cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
        cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel('\N{DEGREE SIGN}C',fontsize=8)

        plt.suptitle('SST - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
        fig.savefig(out_dir+'/SST_'+current_data.strftime("%Y%m%d")+'.png',dpi=resolution_plot,transparent=False, \
                    bbox_inches='tight',pad_inches=0.05)
        plt.close()
        ################## 2X3 SST ###############

        ################## 2X3 SSH ###############
        if (os.path.isfile(archive_hycom_assim)):
           fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
           ax = ax.flatten()
           offset_back = remo.nanmean(ssh_roms[:,:,ind])-remo.nanmean(ssh_hycom_back)

           for i in range(0,6,1):
               ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
               gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
               gl.top_labels = False
               gl.right_labels = False
               gl.xlocator = mticker.FixedLocator(meridians)
               gl.ylocator = mticker.FixedLocator(parallels)
               gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
               gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
               
               if i==0:
                  im = ax[i].contourf(lon_roms, lat_roms, ssh_roms[:,:,ind], np.arange(-0.2,0.51,.025),cmap=cmap,vmin=-0.2, \
                       vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
                             transform=ccrs.PlateCarree())
                  gl.bottom_labels = False
               
               if i==1:
                  im = ax[i].contourf(lon_roms, lat_roms, ssh_hycom_back+offset_back, np.arange(-0.2,0.51,.025),cmap=cmap,vmin=-0.2, \
                       vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                             transform=ccrs.PlateCarree())
                  gl.bottom_labels = False
                  gl.left_labels = False
               
               if i==2:
                  im_mean = ax[i].contourf(lon_roms, lat_roms, ssh_hycom_assim+offset_back, np.arange(-0.2,0.51,.025),cmap=cmap,vmin=-0.2, \
                            vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                             transform=ccrs.PlateCarree())
                  gl.bottom_labels = False
                  gl.left_labels = False
               
               if i==3:
                  im_diff = ax[i].contourf(lon_roms, lat_roms, ssh_hycom_assim-ssh_hycom_back, np.arange(-0.5,0.51,.05),cmap='bwr', \
                            vmin=-0.5,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
                  mean = abs(ssh_hycom_assim-ssh_hycom_back)
                  cond = mean<=0.001
                  mean = np.ma.masked_where(cond,mean)
                  mean = np.mean(mean.compressed())
                  ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
                  del mean

               if i==4:
                  im_diff = ax[i].contourf(lon_roms, lat_roms, ssh_hycom_back+offset_back-ssh_roms[:,:,ind], np.arange(-0.5,0.51,.05), \
                            cmap='bwr',vmin=-0.5,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
                  gl.left_labels = False
                  mean = np.mean(abs(ssh_hycom_back+offset_back-ssh_roms[:,:,ind]).compressed())
                  ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
               if i==5:
                  im = ax[i].contourf(lon_roms, lat_roms, ssh_hycom_assim+offset_back-ssh_roms[:,:,ind], np.arange(-0.5,0.51,.05), \
                       cmap='bwr',vmin=-0.5,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())
                  gl.left_labels = False
                  mean = np.mean(abs(ssh_hycom_assim+offset_back-ssh_roms[:,:,ind]).compressed())
                  ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                             fontweight='bold',transform=ccrs.PlateCarree())

               ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
               ax[i].add_feature(cfeature.LAND)
               ax[i].add_feature(cfeature.COASTLINE)

               ax[i].add_feature(states_provinces, edgecolor='gray')
               

           fig.tight_layout()
           fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

           ## Add a colorbar axis at the bottom of the graph
           ymax = ax[5].get_position().ymax
           ymin = ax[5].get_position().ymin
           cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
           cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
           cbar.ax.tick_params(labelsize=8)
           cbar.ax.set_ylabel('m',fontsize=8)

           ymin = ax[2].get_position().ymin
           ymax = ax[2].get_position().ymax
           cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
           cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
           cbar.ax.tick_params(labelsize=8)
           cbar.ax.set_ylabel('m',fontsize=8)

           plt.suptitle('ADT - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
           fig.savefig(out_dir+'/ADT_'+current_data.strftime("%Y%m%d")+'.png',dpi=resolution_plot,transparent=False, \
                       bbox_inches='tight',pad_inches=0.05)
           plt.close()

        for k in range(k_inic,len(prof[:,0]),1):
            arq_roms = dir_obs_roms_daily+'/ROMS_'+current_data.strftime("%Y%m%d")+'.nc'
            if not (os.path.isfile(arq_roms)):
               print('ROMS FILE: '+arq_roms+' DOES NOT EXIST')
               current_data = current_data + datetime.timedelta(days=inc_tempo)
               continue

            arq_dia_hycom_back_ts = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                                 +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                 +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'

            arq_dia_hycom_back_uv = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                                 +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                 +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'

            arq_dia_hycom_assim_ts = dir_expt+'/output/Ncdf/'+current_data.strftime("%Y%m%d")+'00/archv.' \
                                  +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                  +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'

            arq_dia_hycom_assim_uv = dir_expt+'/output/Ncdf/'+current_data.strftime("%Y%m%d")+'00/archv.' \
                                  +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                  +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'

            if not (os.path.isfile(arq_dia_hycom_back_ts)) or not (os.path.isfile(arq_dia_hycom_assim_ts)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom_back_ts,arq_dia_hycom_assim_ts)
               continue

            if not (os.path.isfile(arq_dia_hycom_back_uv)) or not (os.path.isfile(arq_dia_hycom_assim_uv)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom_back_uv,arq_dia_hycom_assim_uv)
               continue

            nc = Dataset(arq_roms, 'r')
            z_roms = nc.variables['depth'][:] #102
            ind1 = [x for x in enumerate(np.flip(z_roms)) if x[1] <= prof[k,0]][0]
            ind2 = [x for x in enumerate(z_roms) if x[1] >= prof[k,1]][0]
            ind1 = abs(ind1[0]+1-len(z_roms))
            indz_roms = np.arange(ind1,ind2[0]+1,1)
            t_roms = nc.variables['temp'][:] #(564, 502, 102)
            t_roms = t_roms[cond_lat_roms,:,:]
            t_roms = t_roms[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               t_roms = np.squeeze(t_roms[:,:,indz_roms])
            else:
               t_roms = np.mean(t_roms[:,:,indz_roms],axis=-1)
            cond = abs(t_roms) > 999999
            t_roms = np.ma.masked_where(cond,t_roms)

            s_roms = nc.variables['salt'][:] #(564, 502, 102)
            s_roms = s_roms[cond_lat_roms,:,:]
            s_roms = s_roms[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               s_roms = np.squeeze(s_roms[:,:,indz_roms])
            else:
               s_roms = np.mean(s_roms[:,:,indz_roms],axis=-1)
            cond = abs(s_roms) > 999999
            s_roms = np.ma.masked_where(cond,s_roms)

            v_roms = nc.variables['v'][:] #(564, 502, 102)
            v_roms = v_roms[cond_lat_roms,:,:]
            v_roms = v_roms[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               v_roms = np.squeeze(v_roms[:,:,indz_roms])
            else:
               v_roms = np.mean(v_roms[:,:,indz_roms],axis=-1)
            cond = abs(v_roms) > 999999
            v_roms = np.ma.masked_where(cond,v_roms)

            u_roms = nc.variables['u'][:] #(564, 502, 102)
            u_roms = u_roms[cond_lat_roms,:,:]
            u_roms = u_roms[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               u_roms = np.squeeze(u_roms[:,:,indz_roms])
            else:
               u_roms = np.mean(u_roms[:,:,indz_roms],axis=-1)
            cond = abs(u_roms) > 999999
            u_roms = np.ma.masked_where(cond,u_roms)

            vel_roms = np.sqrt(v_roms**2+u_roms**2)
            nc.close()
            del nc

            nc = Dataset(arq_dia_hycom_back_ts, 'r')
            depth_hycom = nc.variables['Depth'][:]
            t_back = nc.variables['temperature'][0,:,:,:] #(33, 564, 502)
            t_back = np.moveaxis(t_back,0,-1)
            s_back = nc.variables['salinity'][0,:,:,:] #(33, 564, 502)
            s_back = np.moveaxis(s_back,0,-1)
            nc.close()
            del nc

            nc = Dataset(arq_dia_hycom_assim_ts, 'r')
            t_assim = nc.variables['temperature'][0,:,:,:] #(33, 564, 502)
            t_assim = np.moveaxis(t_assim,0,-1)
            s_assim = nc.variables['salinity'][0,:,:,:] #(33, 564, 502)
            s_assim = np.moveaxis(s_assim,0,-1)
            nc.close()
            del nc

            nc = Dataset(arq_dia_hycom_back_uv, 'r')
            v_back = nc.variables['v'][0,:,:,:] #(33, 564, 502)
            v_back = np.moveaxis(v_back,0,-1)
            u_back = nc.variables['u'][0,:,:,:] #(33, 564, 502)
            u_back = np.moveaxis(u_back,0,-1)
            nc.close()
            del nc

            nc = Dataset(arq_dia_hycom_assim_uv, 'r')
            v_assim = nc.variables['v'][0,:,:,:] #(33, 564, 502)
            v_assim = np.moveaxis(v_assim,0,-1)
            u_assim = nc.variables['u'][0,:,:,:] #(33, 564, 502)
            u_assim = np.moveaxis(u_assim,0,-1)
            nc.close()
            del nc

            I = interpolate.interp1d(depth_hycom,t_back,axis=-1)
            t_back = I(z_roms[indz_roms].data)
            cond = abs(t_back) > 999999
            t_back = np.ma.masked_where(cond,t_back)
            I = interpolate.interp1d(depth_hycom,t_assim,axis=-1)
            t_assim = I(z_roms[indz_roms].data)
            cond = abs(t_assim) > 999999
            t_assim = np.ma.masked_where(cond,t_assim)
            t_back = np.mean(t_back,axis=2)
            t_assim = np.mean(t_assim,axis=2)

            I = interpolate.interp1d(depth_hycom,s_back,axis=-1)
            s_back = I(z_roms[indz_roms].data)
            cond = abs(s_back) > 999999
            s_back = np.ma.masked_where(cond,s_back)
            I = interpolate.interp1d(depth_hycom,s_assim,axis=-1)
            s_assim = I(z_roms[indz_roms].data)
            cond = abs(s_assim) > 999999
            s_assim = np.ma.masked_where(cond,s_assim)
            s_back = np.mean(s_back,axis=2)
            s_assim = np.mean(s_assim,axis=2)

            I = interpolate.interp1d(depth_hycom,v_back,axis=-1)
            v_back = I(z_roms[indz_roms].data)
            cond = abs(v_back) > 999999
            v_back = np.ma.masked_where(cond,v_back)
            I = interpolate.interp1d(depth_hycom,v_assim,axis=-1)
            v_assim = I(z_roms[indz_roms].data)
            cond = abs(v_assim) > 999999
            v_assim = np.ma.masked_where(cond,v_assim)
            v_back = np.mean(v_back,axis=2)
            v_assim = np.mean(v_assim,axis=2)

            I = interpolate.interp1d(depth_hycom,u_back,axis=-1)
            u_back = I(z_roms[indz_roms].data)
            cond = abs(u_back) > 999999
            u_back = np.ma.masked_where(cond,u_back)
            I = interpolate.interp1d(depth_hycom,u_assim,axis=-1)
            u_assim = I(z_roms[indz_roms].data)
            cond = abs(u_assim) > 999999
            u_assim = np.ma.masked_where(cond,u_assim)
            u_back = np.mean(u_back,axis=2)
            u_assim = np.mean(u_assim,axis=2)

            vel_back = np.sqrt(v_back**2+u_back**2)
            vel_assim = np.sqrt(v_assim**2+u_assim**2)

            ############### 2X3 TEMP ####################################
            if prof[k,1]!=0:
               fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
               ax = ax.flatten()

               if prof[k,1]<=100:
                  t_min = 16
                  t_max = 30
                  delta = 0.5
                  t_range = np.arange(t_min,t_max+delta,delta)
               elif prof[k,1]<=400:
                  t_min = 5
                  t_max = 16
                  delta = 0.5
                  t_range = np.arange(t_min,t_max+delta,delta)
               elif prof[k,1]<=800:
                  t_min = 2
                  t_max = 9
                  delta = 0.25
                  t_range = np.arange(t_min,t_max+delta,delta)
               else:
                  t_min = 1.5
                  t_max = 4.5
                  delta = 0.25
                  t_range = np.arange(t_min,t_max+delta,delta)

               for i in range(0,6,1):
                   ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
                   gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
                   gl.top_labels = False
                   gl.right_labels = False
                   gl.xlocator = mticker.FixedLocator(meridians)
                   gl.ylocator = mticker.FixedLocator(parallels)
                   gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
                   gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
                   
                   if i==0:
                      im = ax[i].contourf(lon_roms, lat_roms, t_roms, t_range,cmap=cmap,vmin=t_min, \
                           vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                      gl.bottom_labels = False
                   
                   if i==1:
                      im = ax[i].contourf(lon_roms, lat_roms, t_back, t_range,cmap=cmap,vmin=t_min, \
                           vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                      gl.bottom_labels = False
                      gl.left_labels = False
                   
                   if i==2:
                      im_mean = ax[i].contourf(lon_roms, lat_roms, t_assim, t_range,cmap=cmap,vmin=t_min, \
                                vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                      gl.left_labels = False
                      gl.bottom_labels = False
                   
                   if i==3:
                      im_diff = ax[i].contourf(lon_roms, lat_roms, t_assim-t_back, np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                                vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                                 fontweight='bold',transform=ccrs.PlateCarree())
                      mean = abs(t_assim-t_back)
                      cond = mean<=0.01
                      mean = np.ma.masked_where(cond,mean)
                      mean = np.mean(mean.compressed())
                      ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                                 fontweight='bold',transform=ccrs.PlateCarree())

                   if i==4:
                      im_diff = ax[i].contourf(lon_roms, lat_roms, t_back-t_roms, np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                                vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
                                 fontweight='bold',transform=ccrs.PlateCarree())
                      gl.left_labels = False
                      mean = np.mean(abs(t_back-t_roms).compressed())
                      ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                                 fontweight='bold',transform=ccrs.PlateCarree())

                   if i==5:
                      im = ax[i].contourf(lon_roms, lat_roms, t_assim-t_roms, np.arange(-3,3.1,.25),cmap='bwr',vmin=-3, \
                           vmax=3,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                      ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', \
                                 fontsize=7, fontweight='bold',transform=ccrs.PlateCarree())
                      gl.left_labels = False
                      mean = np.mean(abs(t_assim-t_roms).compressed())
                      ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
                                 fontweight='bold',transform=ccrs.PlateCarree())

                   ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
                   ax[i].add_feature(cfeature.LAND)
                   ax[i].add_feature(cfeature.COASTLINE)

                   ax[i].add_feature(states_provinces, edgecolor='gray')
                   
               fig.tight_layout()
               fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

               ## Add a colorbar axis at the bottom of the graph
               ymax = ax[5].get_position().ymax
               ymin = ax[5].get_position().ymin
               cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
               cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
               cbar.ax.tick_params(labelsize=8)
               cbar.ax.set_ylabel('\N{DEGREE SIGN}C',fontsize=8)

               ymin = ax[2].get_position().ymin
               ymax = ax[2].get_position().ymax
               cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
               cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
               cbar.ax.tick_params(labelsize=8)
               cbar.ax.set_ylabel('\N{DEGREE SIGN}C',fontsize=8)

               plt.suptitle('TEMP DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
               fig.savefig(out_dir+'/TEMP_'+current_data.strftime("%Y%m%d")+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png', \
                           dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
               plt.close()
            ################## 2X3 TEMP ###############

            ################# 2X3 SALT ####################################
            ##fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
            ##ax = ax.flatten()

            ##if prof[k,1]<=100:
            ##   t_min = 35
            ##   t_max = 37.5
            ##   delta = 0.050
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##elif prof[k,1]<=400:
            ##   t_min = 34
            ##   t_max = 36.0
            ##   delta = 0.050
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##elif prof[k,1]<=800:
            ##   t_min = 34
            ##   t_max = 35.0
            ##   delta = 0.025
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##else:
            ##   t_min = 34.5
            ##   t_max = 35.5
            ##   delta = 0.025
            ##   t_range = np.arange(t_min,t_max+delta,delta)

            ##for i in range(0,6,1):
            ##    ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
            ##    gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
            ##    gl.top_labels = False
            ##    gl.right_labels = False
            ##    gl.xlocator = mticker.FixedLocator(meridians)
            ##    gl.ylocator = mticker.FixedLocator(parallels)
            ##    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            ##    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            ##    
            ##    if i==0:
            ##       im = ax[i].contourf(lon_roms, lat_roms, s_roms, t_range,cmap=cmap,vmin=t_min, \
            ##            vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.bottom_labels = False
            ##    
            ##    if i==1:
            ##       im = ax[i].contourf(lon_roms, lat_roms, s_back, t_range,cmap=cmap,vmin=t_min, \
            ##            vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.bottom_labels = False
            ##       gl.left_labels = False
            ##    
            ##    if i==2:
            ##       im_mean = ax[i].contourf(lon_roms, lat_roms, s_assim, t_range,cmap=cmap,vmin=t_min, \
            ##                 vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       gl.bottom_labels = False
            ##    
            ##    if i==3:
            ##       im_diff = ax[i].contourf(lon_roms, lat_roms, s_assim-s_back, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##                 vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())
            ##       mean = np.mean(abs(s_assim-s_back).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    if i==4:
            ##       im_diff = ax[i].contourf(lon_roms, lat_roms, s_back-s_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##                 vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       mean = np.mean(abs(s_back-s_roms).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    if i==5:
            ##       im = ax[i].contourf(lon_roms, lat_roms, s_assim-s_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##            vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', \
            ##                  fontsize=7, fontweight='bold',transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       mean = np.mean(abs(s_assim-s_roms).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,2)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
            ##    ax[i].add_feature(cfeature.LAND)
            ##    ax[i].add_feature(cfeature.COASTLINE)

            ##    ax[i].add_feature(states_provinces, edgecolor='gray')
            ##    
            ##fig.tight_layout()
            ##fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

            #### Add a colorbar axis at the bottom of the graph
            ##ymax = ax[5].get_position().ymax
            ##ymin = ax[5].get_position().ymin
            ##cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            ##cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
            ##cbar.ax.tick_params(labelsize=8)
            ##cbar.ax.set_ylabel('PSU',fontsize=8)

            ##ymin = ax[2].get_position().ymin
            ##ymax = ax[2].get_position().ymax
            ##cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            ##cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
            ##cbar.ax.tick_params(labelsize=8)
            ##cbar.ax.set_ylabel('PSU',fontsize=8)

            ##plt.suptitle('SALT DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
            ##fig.savefig(out_dir+'/SALT_'+current_data.strftime("%Y%m%d")+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png', \
            ##            dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
            ##plt.close()
            #################### 2X3 SALT ###############

            ################# 2X3 V ####################################
            ##fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
            ##ax = ax.flatten()

            ##if prof[k,1]<=100:
            ##   t_min = -1
            ##   t_max = 1
            ##   delta = 0.100
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##elif prof[k,1]<=400:
            ##   t_min = -1
            ##   t_max = 1
            ##   delta = 0.100
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##elif prof[k,1]<=800:
            ##   t_min = -0.8
            ##   t_max = 0.8
            ##   delta = 0.100
            ##   t_range = np.arange(t_min,t_max+delta,delta)
            ##else:
            ##   t_min = -0.6
            ##   t_max = 0.6
            ##   delta = 0.050
            ##   t_range = np.arange(t_min,t_max+delta,delta)

            ##for i in range(0,6,1):
            ##    ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
            ##    gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
            ##    gl.top_labels = False
            ##    gl.right_labels = False
            ##    gl.xlocator = mticker.FixedLocator(meridians)
            ##    gl.ylocator = mticker.FixedLocator(parallels)
            ##    gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            ##    gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            ##    
            ##    if i==0:
            ##       im = ax[i].contourf(lon_roms, lat_roms, v_roms, t_range,cmap=cmap,vmin=t_min, \
            ##            vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.bottom_labels = False
            ##    
            ##    if i==1:
            ##       im = ax[i].contourf(lon_roms, lat_roms, v_back, t_range,cmap=cmap,vmin=t_min, \
            ##            vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.bottom_labels = False
            ##       gl.left_labels = False
            ##    
            ##    if i==2:
            ##       im_mean = ax[i].contourf(lon_roms, lat_roms, v_assim, t_range,cmap=cmap,vmin=t_min, \
            ##                 vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
            ##                  transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       gl.bottom_labels = False
            ##    
            ##    if i==3:
            ##       im_diff = ax[i].contourf(lon_roms, lat_roms, v_assim-v_back, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##                 vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())
            ##       mean = np.mean(abs(v_assim-v_back).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    if i==4:
            ##       im_diff = ax[i].contourf(lon_roms, lat_roms, v_back-v_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##                 vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       mean = np.mean(abs(v_back-v_roms).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    if i==5:
            ##       im = ax[i].contourf(lon_roms, lat_roms, v_assim-v_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
            ##            vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', \
            ##                  fontsize=7, fontweight='bold',transform=ccrs.PlateCarree())
            ##       gl.left_labels = False
            ##       mean = np.mean(abs(v_assim-v_roms).compressed())
            ##       ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
            ##                  fontweight='bold',transform=ccrs.PlateCarree())

            ##    ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
            ##    ax[i].add_feature(cfeature.LAND)
            ##    ax[i].add_feature(cfeature.COASTLINE)

            ##    ax[i].add_feature(states_provinces, edgecolor='gray')
            ##    
            ##fig.tight_layout()
            ##fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

            #### Add a colorbar axis at the bottom of the graph
            ##ymax = ax[5].get_position().ymax
            ##ymin = ax[5].get_position().ymin
            ##cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            ##cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
            ##cbar.ax.tick_params(labelsize=8)
            ##cbar.ax.set_ylabel('m/s',fontsize=8)

            ##ymin = ax[2].get_position().ymin
            ##ymax = ax[2].get_position().ymax
            ##cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            ##cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
            ##cbar.ax.tick_params(labelsize=8)
            ##cbar.ax.set_ylabel('m/s',fontsize=8)

            ##plt.suptitle('V VEL DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
            ##fig.savefig(out_dir+'/V_VEL_'+current_data.strftime("%Y%m%d")+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png', \
            ##            dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
            ##plt.close()
            #################### 2X3 V ###############

            ############### 2X3 VEL RES ####################################
            under='#FFFFFF'
            cmap.set_under(under)
            fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
            ax = ax.flatten()

            inc_latlon = 17
            if prof[k,1]<=100:
               t_min = 0.2
               t_max = 1.3
               delta = 0.050
               scale = 12
               t_range = np.arange(t_min,t_max+delta,delta)
            elif prof[k,1]<=400:
               t_min = 0.15
               t_max = 1
               delta = 0.050
               scale = 9
               t_range = np.arange(t_min,t_max+delta,delta)
            elif prof[k,1]<=800:
               t_min = 0.1
               t_max = 0.6
               delta = 0.025
               scale = 9
               t_range = np.arange(t_min,t_max+delta,delta)
            else:
               t_min = 0.1
               t_max = 0.5
               delta = 0.025
               scale = 5
               t_range = np.arange(t_min,t_max+delta,delta)

            for i in range(0,6,1):
                ax[i].set_extent([min_lon_roms, max_lon_roms, min_lat_roms, max_lat_roms])
                gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='gray', alpha=0.6, linestyle='--')
                gl.top_labels = False
                gl.right_labels = False
                gl.xlocator = mticker.FixedLocator(meridians)
                gl.ylocator = mticker.FixedLocator(parallels)
                gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
                gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
                
                if i==0:
                   im = ax[i].contourf(lon_roms, lat_roms, vel_roms, t_range,cmap=cmap,vmin=t_min, \
                        vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'Nature Run', horizontalalignment='left', fontsize=7, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   x,y = np.meshgrid(lon_roms[0:-1:inc_latlon],lat_roms[0:-1:inc_latlon])
                   ax[i].quiver(x, y, u_roms[0:-1:inc_latlon,0:-1:inc_latlon], v_roms[0:-1:inc_latlon,0:-1:inc_latlon], \
                                transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5, alpha=0.7)
                   x,y = np.meshgrid(min_lon_roms+1,max_lat_roms-3.4)
                   ax[i].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5)
                   ax[i].text(min_lon_roms+1, max_lat_roms-4.8, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   gl.bottom_labels = False
                
                if i==1:
                   im = ax[i].contourf(lon_roms, lat_roms, vel_back, t_range,cmap=cmap,vmin=t_min, \
                        vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK', horizontalalignment='left', fontsize=7, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   x,y = np.meshgrid(lon_roms[0:-1:inc_latlon],lat_roms[0:-1:inc_latlon])
                   ax[i].quiver(x, y, u_back[0:-1:inc_latlon,0:-1:inc_latlon], v_back[0:-1:inc_latlon,0:-1:inc_latlon], \
                                transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5, alpha=0.7)
                   x,y = np.meshgrid(min_lon_roms+1,max_lat_roms-3.4)
                   ax[i].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5)
                   ax[i].text(min_lon_roms+1, max_lat_roms-4.8, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   gl.bottom_labels = False
                   gl.left_labels = False
                
                if i==2:
                   im_mean = ax[i].contourf(lon_roms, lat_roms, vel_assim, t_range,cmap=cmap,vmin=t_min, \
                             vmax=t_max,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA', horizontalalignment='left', fontsize=7, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   x,y = np.meshgrid(lon_roms[0:-1:inc_latlon],lat_roms[0:-1:inc_latlon])
                   ax[i].quiver(x, y, u_assim[0:-1:inc_latlon,0:-1:inc_latlon], v_assim[0:-1:inc_latlon,0:-1:inc_latlon], \
                                transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5, alpha=0.7)
                   x,y = np.meshgrid(min_lon_roms+1,max_lat_roms-3.4)
                   ax[i].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=scale,width=0.007,headaxislength=5)
                   ax[i].text(min_lon_roms+1, max_lat_roms-4.8, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                              transform=ccrs.PlateCarree())

                   gl.left_labels = False
                   gl.bottom_labels = False
                
                if i==3:
                   im_diff = ax[i].contourf(lon_roms, lat_roms, vel_assim-vel_back, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
                             vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - BACK', horizontalalignment='left', fontsize=7, \
                              fontweight='bold',transform=ccrs.PlateCarree())
                   mean = np.mean(abs(vel_assim-vel_back).compressed())
                   ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                              fontweight='bold',transform=ccrs.PlateCarree())

                if i==4:
                   im_diff = ax[i].contourf(lon_roms, lat_roms, vel_back-vel_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
                             vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'BACK - Nature Run', horizontalalignment='left', fontsize=7, \
                              fontweight='bold',transform=ccrs.PlateCarree())
                   gl.left_labels = False
                   mean = np.mean(abs(vel_back-vel_roms).compressed())
                   ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                              fontweight='bold',transform=ccrs.PlateCarree())

                if i==5:
                   im = ax[i].contourf(lon_roms, lat_roms, vel_assim-vel_roms, np.arange(-0.5,0.55,0.05),cmap='bwr',vmin=-0.5, \
                        vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                   ax[i].text(min_lon_roms+1, max_lat_roms-3, 'ANA - Nature Run', horizontalalignment='left', \
                              fontsize=7, fontweight='bold',transform=ccrs.PlateCarree())
                   gl.left_labels = False
                   mean = np.mean(abs(vel_assim-vel_roms).compressed())
                   ax[i].text(min_lon_roms+1, max_lat_roms-5, str(round(mean,3)), horizontalalignment='left', fontsize=7, \
                              fontweight='bold',transform=ccrs.PlateCarree())

                ax[i].add_feature(cfeature.GSHHSFeature(scale='high'))
                ax[i].add_feature(cfeature.LAND)
                ax[i].add_feature(cfeature.COASTLINE)

                ax[i].add_feature(states_provinces, edgecolor='gray')
                
            fig.tight_layout()
            fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace=-0.14)

            ## Add a colorbar axis at the bottom of the graph
            ymax = ax[5].get_position().ymax
            ymin = ax[5].get_position().ymin
            cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            cbar=fig.colorbar(im_diff, cax=cbar_ax,orientation='vertical')
            cbar.ax.tick_params(labelsize=8)
            cbar.ax.set_ylabel('m/s',fontsize=8)

            ymin = ax[2].get_position().ymin
            ymax = ax[2].get_position().ymax
            cbar_ax = fig.add_axes([0.99, ymin, 0.02, ymax-ymin])
            cbar=fig.colorbar(im_mean, cax=cbar_ax,orientation='vertical')
            cbar.ax.tick_params(labelsize=8)
            cbar.ax.set_ylabel('m/s',fontsize=8)

            plt.suptitle('VEL RES DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m - '+current_data.strftime("%Y%m%d")+' ('+legenda[rod]+')', y=ymax+0.03)
            fig.savefig(out_dir+'/VEL_RES_'+current_data.strftime("%Y%m%d")+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png', \
                        dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
            plt.close()
            ################## 2X3 V ###############

        current_data = current_data + datetime.timedelta(days=inc_tempo)

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

