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
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import warnings
from calendar import monthrange
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

cfg_file='dirs_plot.txt'
opt=open(cfg_file, "r")
content_list = opt.readlines()
input_dir = str(content_list[0].rstrip('\n'))
fig_output_dir = str(content_list[1].rstrip('\n'))

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

months = [1,2,3,4,5,6,7,8,9,10,11,12]
#months = [7,8,9,10,11]
#months = [1,2,3,4,5,6,7,8]
prof = np.array([[0,0], [100,100], [200,200], [400,400], [800,800]])
#prof = np.array([[0,0]])

#under='#000033'
under='#FFFFFF'
over='#330000'
bad='#FFFFFF'
cmap=matplotlib.cm.get_cmap('jet')

for k in range(0,len(prof[:,0]),1):
    print('DEPTH: '+str(prof[k,0])+"-"+str(prof[k,1])+'m')
    if '004j' in rodada[0]:
       min_lat = -33.74+2
       max_lat = -12.32-2
       min_lon = -53.33
       max_lon = -32.46-2
    
       lat_step_plot = 4
       lon_step_plot = 4
       resolution_plot = 300

    elif '008d' in rodada[0]:
       min_lat = -79.54
       max_lat =  50.27
       min_lon = -98.00
       max_lon =  45.00

       subs = np.array([[-53, -37, -34, -19], \
                        [-62, -44, -45, -30], \
                        [-42, -30, -19,  -8]])

       lat_step_plot = 3
       lon_step_plot = 3
       resolution_plot = 300
    
    elif '008i' in rodada[0]:
       min_lat = -47.93
       max_lat =  10.11
       min_lon = -70.00
       max_lon = -17.75

       subs = np.array([[-53, -37, -34, -19], \
                        [-62, -44, -45, -30], \
                        [-42, -30, -19,  -8]])

       lat_step_plot = 3
       lon_step_plot = 3
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
        current_data = data_inicial
        y = data_inicial.year
        counter = 0
    
        if rod>0:
           u_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
           u_hycom_mean = u_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
           v_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
           v_hycom_mean = v_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
           vel_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
           vel_hycom_mean = vel_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
    
        while current_data<=data_final:
            while current_data.month != months[0]:
               days = monthrange(current_data.year,current_data.month)[1]
               current_data = current_data + datetime.timedelta(days=days)
            for m in list(months):
                days = monthrange(current_data.year,current_data.month)[1]
                filename = input_dir+"/"+rodada[rod]+"/VEL_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                           str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+str(days).zfill(2)+"_"+dominio_str+"_"+ \
                           str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                if not (os.path.isfile(filename)):
                   #print(filename, ' NAO ENCONTRADO')
                   filename = input_dir+"/"+rodada[rod]+"/VEL_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                              str(m).zfill(2)+"01-"+data_final.strftime("%Y%m%d")+"_"+dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                   #print('USANDO ',filename)
                if not (os.path.isfile(filename)):
                   print(filename, ' NAO ENCONTRADO')
                   exit()
                mat_contents = scipy.io.loadmat(filename)
                if not 'lat_hycom' in locals():
                   locals()['lat_hycom'] = mat_contents['lat_hycom'].flatten()
                   locals()['lon_hycom'] = mat_contents['lon_hycom'].flatten()
    
                if y == data_inicial.year and m==(months[0]):
                   u_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   u_hycom_mean = u_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
                   v_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   v_hycom_mean = v_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
                   vel_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   vel_hycom_mean = vel_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
    
                u_hycom = mat_contents['u_hycom_mean_'+rodada[rod]]
                v_hycom = mat_contents['v_hycom_mean_'+rodada[rod]]
                vel_hycom = mat_contents['vel_hycom_mean_'+rodada[rod]]
                cond = abs(u_hycom) > 99999
                u_hycom = np.ma.masked_where(cond,u_hycom)
                cond = abs(v_hycom) > 99999
                v_hycom = np.ma.masked_where(cond,v_hycom)
                cond = abs(vel_hycom) > 99999
                vel_hycom = np.ma.masked_where(cond,vel_hycom)
                u_hycom_mean = u_hycom_mean + u_hycom
                v_hycom_mean = v_hycom_mean + v_hycom
                vel_hycom_mean = vel_hycom_mean + vel_hycom
    
                del u_hycom, v_hycom, vel_hycom, mat_contents, cond
                counter = counter+1
                current_data = current_data + datetime.timedelta(days=days)
    
        locals()['u_hycom_mean_'+rodada[rod]] = u_hycom_mean/counter
        locals()['v_hycom_mean_'+rodada[rod]] = v_hycom_mean/counter
        locals()['vel_hycom_mean_'+rodada[rod]] = vel_hycom_mean/counter
        del u_hycom_mean, v_hycom_mean, vel_hycom_mean, counter, y, m, days, current_data
    
    for region in range(0,subs.shape[0]):
        cond_lat_hycom = (lat_hycom<=subs[region,3]) & (lat_hycom>=subs[region,2])
        lat_hycom_plot = lat_hycom[cond_lat_hycom]
        cond_lon_hycom = (lon_hycom<=subs[region,1]) & (lon_hycom>=subs[region,0])
        lon_hycom_plot = lon_hycom[cond_lon_hycom]
        min_lat_plot = min(lat_hycom_plot)
        max_lat_plot = max(lat_hycom_plot)
        min_lon_plot = min(lon_hycom_plot)
        max_lon_plot = max(lon_hycom_plot)
        
        if min_lat_plot<0:
           s1 = "S"
        else:
           s1 = "N"
        if max_lat_plot<0:
           s2 = "S"
        else:
           s2 = "N"
        if min_lon_plot<0:
           s3 = "W"
        else:
           s3 = "E"
        if max_lon_plot<0:
           s4 = "W"
        else:
           s4 = "E"
        dominio_str_plot = f"{abs(min_lat_plot):.2f}"+s1+"-"+f"{abs(max_lat_plot):.2f}"+s2+"_" + f"{abs(min_lon_plot):.2f}"+s3+"-"+f"{abs(max_lon_plot):.2f}"+s4

        ####### VEL MEAN MAP ###########
        
        if len(rodada)<=3:
           fig, ax = plt.subplots(nrows=1,ncols=len(rodada),subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        elif len(rodada)==4:
           fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
        
        if len(rodada)>1:
           ax = ax.flatten()
        all_rod_save = rodada[0]
        for rod in range(1,len(rodada),1):
            all_rod_save = all_rod_save+"_"+rodada[rod]
        out_dir = fig_output_dir+"/"+all_rod_save+"/CAMPOS_VEL"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)
        
        parallels = np.arange(int(min_lat_plot),int(max_lat_plot)+0.1,lat_step_plot)
        meridians = np.arange(int(min_lon_plot),int(max_lon_plot)+0.1,lon_step_plot)
        cmap.set_under(under)
        cmap.set_over(over)
        cmap.set_bad(bad)

        for rod in range(0,len(rodada),1):
        
            ####### PLOT MEAN MAP ###########
        
            vel_hycom_plot = locals()['vel_hycom_mean_'+rodada[rod]][cond_lat_hycom,:]
            vel_hycom_plot = vel_hycom_plot[:,cond_lon_hycom]
            u_hycom_plot = locals()['u_hycom_mean_'+rodada[rod]][cond_lat_hycom,:]
            u_hycom_plot = u_hycom_plot[:,cond_lon_hycom]
            v_hycom_plot = locals()['v_hycom_mean_'+rodada[rod]][cond_lat_hycom,:]
            v_hycom_plot = v_hycom_plot[:,cond_lon_hycom]

            if len(rodada)==1:
               ax.set_extent([min_lon_plot, max_lon_plot, min_lat_plot, max_lat_plot], crs=ccrs.PlateCarree())
               teste = ax
            else:
               teste = ax[rod]
               teste.set_extent([min_lon_plot, max_lon_plot, min_lat_plot, max_lat_plot], crs=ccrs.PlateCarree())

            if (subs[region,0]==-53 and subs[region,1]==-37 and subs[region,2]==-34 and subs[region,3]==-19): #[-53, -37, -34, -19] CB
               if prof[k,1]<100:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.2,0.725,.025),cmap=cmap, \
                                       vmin=0.2,vmax=0.7,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.2,0.71,0.1)
               elif prof[k,1]<400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.15,0.625,.025),cmap=cmap, \
                                          vmin=0.15,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.61,0.1),[0],0.15)
               elif prof[k,1]==400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.10,0.525,.025), \
                                          cmap=cmap,vmin=0.10,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.1,0.51,0.1)
                  
               elif prof[k,1]>400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.10,0.41,.010), \
                                          cmap=cmap,vmin=0.10,vmax=0.40,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.1,0.40,0.05)

               if prof[k,1]<400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-1.5)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)
               elif prof[k,1]>=400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=5.0,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-1.5)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=5.00,width=0.007,headaxislength=5)

               if len(rodada)==2:
                  teste.text(min_lon_plot+1, max_lat_plot-1.2, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-2.4, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
               elif len(rodada)==4:
                  teste.text(min_lon_plot+1, max_lat_plot-1.0, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-2.4, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())

               if rod==len(rodada)-1:
                  fig.tight_layout()
                  if len(rodada)==2:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
                     y_title=0.840
                     cbar_ax = [0.990, 0.200, 0.02, 0.590]
                  elif len(rodada)==4:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.16, hspace= 0.02)
                     y_title=1.025
                     cbar_ax = [0.860, 0.01, 0.02, 0.980]

            ################ [-53, -37, -34, -19] CB ################################

            elif (subs[region,0]==-42 and subs[region,1]==-30 and subs[region,2]==-19 and subs[region,3]==-8): #[-42, -30, -19,  -8] BiCSE
               if prof[k,1]<100:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.2,0.62,.020),cmap=cmap, \
                                       vmin=0.2,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.2,0.61,0.1)
               elif prof[k,1]<400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.16,0.62,.020),cmap=cmap, \
                                          vmin=0.16,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.61,0.1),[0],0.16)
               elif prof[k,1]==400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.16,0.825,.025), \
                                          cmap=cmap,vmin=0.16,vmax=0.8,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.81,0.1),[0],0.16)

               elif prof[k,1]>400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.10,0.41,.010), \
                                          cmap=cmap,vmin=0.10,vmax=0.40,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.1,0.40,0.05)                

               if prof[k,1]<400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-1.5)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)
               elif prof[k,1]>=400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-1.5)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)

               if len(rodada)==2:
                  teste.text(min_lon_plot+1, max_lat_plot-1.2, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-2.1, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
               elif len(rodada)==4:
                  teste.text(min_lon_plot+1, max_lat_plot-1.0, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-2.2, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())

               if rod==len(rodada)-1:
                  fig.tight_layout()
                  if len(rodada)==2:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
                     y_title=0.820
                     cbar_ax = [0.990, 0.200, 0.02, 0.590]
                  elif len(rodada)==4:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.13, hspace= 0.03)
                     y_title=1.025
                     cbar_ax = [0.87, 0.005, 0.02, 0.980]


            ################## [-42, -30, -19,  -8] BiCSE ######################

            elif (subs[region,0]==-62 and subs[region,1]==-44 and subs[region,2]==-45 and subs[region,3]==-30): #[-62, -44, -45, -30] CBM
               if prof[k,1]<100:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.26,0.92,.020),cmap=cmap, \
                                       vmin=0.26,vmax=0.9,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.30,0.91,0.1),[0],0.26)
               elif prof[k,1]<400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.20,0.92,.020),cmap=cmap, \
                                          vmin=0.20,vmax=0.9,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.2,0.91,0.1)
               elif prof[k,1]==400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.15,0.825,.025), \
                                          cmap=cmap,vmin=0.15,vmax=0.8,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.81,0.1),[0],0.15)

               elif prof[k,1]>400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.20,0.71,.010), \
                                          cmap=cmap,vmin=0.20,vmax=0.70,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.2,0.75,0.05)                

               if prof[k,1]<400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:8],lat_hycom_plot[0:-1:8])
                  teste.quiver(x,y,u_hycom_plot[0:-1:8,0:-1:8],v_hycom_plot[0:-1:8,0:-1:8], \
                                   transform=ccrs.PlateCarree(),scale=10,width=0.007,headaxislength=5, alpha=0.5)
                  if len(rodada)==2:
                     x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-2.5)
                  else:
                     print('DEFINE VECTOR POSITION FOR CBM')
                     exit()
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=10.00,width=0.007,headaxislength=5)
               elif prof[k,1]>=400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:8],lat_hycom_plot[0:-1:8])
                  teste.quiver(x,y,u_hycom_plot[0:-1:8,0:-1:8],v_hycom_plot[0:-1:8,0:-1:8], \
                                   transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
                  if len(rodada)==2:
                     x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-2.5)
                  else:
                     print('DEFINE VECTOR POSITION FOR CBM')
                     exit()
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)

               if len(rodada)==2:
                  teste.text(min_lon_plot+1, max_lat_plot-1.5, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-2.2, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
               elif len(rodada)==4:
                  teste.text(min_lon_plot+1, max_lat_plot-2.8, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())
                  teste.text(min_lon_plot+1, max_lat_plot-4.2, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                                 transform=ccrs.PlateCarree())

               if rod==len(rodada)-1:
                  fig.tight_layout()
                  if len(rodada)==2:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
                     y_title=0.800
                     cbar_ax = [0.990, 0.225, 0.02, 0.525]
                  elif len(rodada)==4:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace= 0.02, hspace= 0.02)                
                     y_title=1.020
                     cbar_ax = [0.92, 0.018, 0.02, 0.972]

            ################## [-62, -44, -45, -30] CBM ######################
            else:
               if len(rodada)<=3:
                  y_title=0.825
                  cbar_ax = [0.87, 0.02, 0.02, 0.97]
               elif len(rodada)==4:
                  y_title=0.825
                  cbar_ax = [0.80, 0.015, 0.02, 0.97]

               if prof[k,1]<100:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.2,0.62,.020),cmap=cmap, \
                                       vmin=0.2,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.2,0.61,0.1)
               elif prof[k,1]<400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.16,0.62,.020),cmap=cmap, \
                                          vmin=0.16,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.61,0.1),[0],0.16)
               elif prof[k,1]==400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.16,0.52,.020), \
                                          cmap=cmap,vmin=0.16,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.insert(np.arange(0.2,0.51,0.1),[0],0.16)
             
               elif prof[k,1]>400:
                  im = teste.contourf(lon_hycom_plot, lat_hycom_plot, vel_hycom_plot, np.arange(0.10,0.36,.010), \
                                          cmap=cmap,vmin=0.10,vmax=0.35,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
                  cbar_ticks = np.arange(0.1,0.40,0.05) 

               if prof[k,1]<400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-3.3)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)
               elif prof[k,1]>=400:
                  x,y = np.meshgrid(lon_hycom_plot[0:-1:5],lat_hycom_plot[0:-1:5])
                  teste.quiver(x,y,u_hycom_plot[0:-1:5,0:-1:5],v_hycom_plot[0:-1:5,0:-1:5], \
                                   transform=ccrs.PlateCarree(),scale=5,width=0.007,headaxislength=5, alpha=0.5)
                  x,y = np.meshgrid(min_lon_plot+1,max_lat_plot-3.3)
                  teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=5.00,width=0.007,headaxislength=5)

               if rod==len(rodada)-1:
                  fig.tight_layout()
                  if len(rodada)<=2:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
                  else:
                     fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.02, hspace=-0.02)                

            teste.add_feature(cfeature.LAND)
            teste.add_feature(cfeature.COASTLINE, linewidth=0.5)
            states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines', \
                               scale='10m',facecolor='none')
            teste.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
        
            gl = teste.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.right_labels = False
            if len(rodada)<=2:
                if rod==0: 
                   gl.left_labels = True 
                elif rod==1: 
                   gl.left_labels = False
            elif len(rodada)==3:
               if rod==0:
                  gl.left_labels = False
                  gl.bottom_labels = False
               elif rod==2:
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
            gl.ylabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
            gl.xlabel_style = {'size': 8, 'color': 'black', 'weight': 'bold'}
        
            teste.tick_params(axis='both', labelsize=8)
        
            if rod==0:
                all_rod = rodada[rod]
            else:
                all_rod = all_rod+"_"+rodada[rod]

            if rod==len(rodada)-1:
               cbar_ax = fig.add_axes(cbar_ax)
               if len(months)==12:
                  plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m ('+data_inicial.strftime("%Y%m%d")+' - '+ \
                               data_final.strftime("%Y%m%d")+')', y=y_title)
               else:
                  plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m (MONTHS: '+str(months[0]).zfill(2)+'-'+ \
                               str(months[-1]).zfill(2)+'; YEARS: '+str(data_inicial.year)+'-'+str(data_final.year)+')', y=y_title)

               cbar=fig.colorbar(im, ticks = cbar_ticks, cax=cbar_ax,orientation='vertical')
               cbar.ax.set_ylabel('m/s')
        
               if len(months)==12:
                  fig.savefig(out_dir+'/VEL_MEAN_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                              '_'+dominio_str_plot+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False, \
                              bbox_inches='tight',pad_inches=0.05)
               else:
                  fig.savefig(out_dir+'/VEL_MEAN_'+all_rod+'_months_'+str(months[0]).zfill(2)+'-'+str(months[-1]).zfill(2)+'_YEARS_'+ \
                              str(data_inicial.year)+'-'+str(data_final.year)+'_'+dominio_str_plot+'_DEPTH_'+str(prof[k,0])+"-"+ \
                              str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
        
        del all_rod,fig,ax

exit()

