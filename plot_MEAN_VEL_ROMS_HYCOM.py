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

input_dir = '/prj/rodas/filipe.costa/resultados'
fig_output_dir = '/prj/rodas/filipe.costa/resultados/figs'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

months = [1,2,3,4,5,6,7,8,9,10,11,12]
months = [7,8,9,10,11,12,1,2,3,4,5,6]
prof = np.array([[0,0], [100,100], [400,400], [800,800]])

#under='#000033'
under='#FFFFFF'
over='#330000'
bad='#FFFFFF'
cmap=matplotlib.cm.get_cmap('jet')

for k in range(0,len(prof[:,0]),1):
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
        current_data = data_inicial
        y = data_inicial.year
        counter = 0
    
        if rod>0:
           u_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
           u_hycom_mean = u_hycom_mean.reshape(len(lat_roms),len(lon_roms))
           v_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
           v_hycom_mean = v_hycom_mean.reshape(len(lat_roms),len(lon_roms))
           vel_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
           vel_hycom_mean = vel_hycom_mean.reshape(len(lat_roms),len(lon_roms))
    
        while current_data<=data_final:
            while current_data.month != months[0]:
               days = monthrange(current_data.year,current_data.month)[1]
               current_data = current_data + datetime.timedelta(days=days)
            for m in list(months):
                days = monthrange(current_data.year,current_data.month)[1]
                filename = input_dir+"/"+rodada[rod]+"/VEL_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                           str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+str(days).zfill(2)+"_"+dominio_str+"_"+ \
                           str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                if not (os.path.isfile(filename)):
                   print(filename, ' NAO ENCONTRADO')
                   filename = input_dir+"/"+rodada[rod]+"/VEL_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                              str(m).zfill(2)+"01-"+data_final.strftime("%Y%m%d")+"_"+dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                   print('USANDO ',filename)
                if not (os.path.isfile(filename)):
                   print(filename, ' NAO ENCONTRADO')
                   exit()
                mat_contents = scipy.io.loadmat(filename)
                if not 'lat_roms' in locals():
                   locals()['lat_roms'] = mat_contents['lat_roms'].flatten()
                   locals()['lon_roms'] = mat_contents['lon_roms'].flatten()
    
                   u_roms_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   u_roms_mean = u_roms_mean.reshape(len(lat_roms),len(lon_roms))
                   v_roms_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   v_roms_mean = v_roms_mean.reshape(len(lat_roms),len(lon_roms))
                   vel_roms_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   vel_roms_mean = vel_roms_mean.reshape(len(lat_roms),len(lon_roms))
    
                if y == data_inicial.year and m==(months[0]):
                   u_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   u_hycom_mean = u_hycom_mean.reshape(len(lat_roms),len(lon_roms))
                   v_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   v_hycom_mean = v_hycom_mean.reshape(len(lat_roms),len(lon_roms))
                   vel_hycom_mean = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)), mask=False)
                   vel_hycom_mean = vel_hycom_mean.reshape(len(lat_roms),len(lon_roms))
    
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
    
                if rod==0:
                   u_roms = mat_contents['u_roms_mean']
                   v_roms = mat_contents['v_roms_mean']
                   vel_roms = mat_contents['vel_roms_mean']
                   cond = abs(u_roms) > 99999
                   u_roms = np.ma.masked_where(cond,u_roms)
                   cond = abs(v_roms) > 99999
                   v_roms = np.ma.masked_where(cond,v_roms)
                   cond = abs(vel_roms) > 99999
                   vel_roms = np.ma.masked_where(cond,vel_roms)
                   u_roms_mean = u_roms_mean + u_roms
                   v_roms_mean = v_roms_mean + v_roms
                   vel_roms_mean = vel_roms_mean + vel_roms
                   del u_roms, v_roms, vel_roms
    
                del u_hycom, v_hycom, vel_hycom, mat_contents, cond
                counter = counter+1
                current_data = current_data + datetime.timedelta(days=days)
    
        locals()['u_hycom_mean_'+rodada[rod]] = u_hycom_mean/counter
        locals()['v_hycom_mean_'+rodada[rod]] = v_hycom_mean/counter
        locals()['vel_hycom_mean_'+rodada[rod]] = vel_hycom_mean/counter
        if rod==0:
           u_roms_mean = u_roms_mean/counter
           v_roms_mean = v_roms_mean/counter
           vel_roms_mean = vel_roms_mean/counter
        del u_hycom_mean, v_hycom_mean, vel_hycom_mean, counter, y, m, days, current_data
    
    ####### VEL MEAN MAP ###########
    
    if len(rodada)<=2:
       fig, ax = plt.subplots(nrows=1,ncols=len(rodada)+1,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
    elif len(rodada)==3:
       fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},sharey=True,figsize=(6.4, 4.8))
    
    ax = ax.flatten()
    all_rod_save = rodada[0]
    for rod in range(1,len(rodada),1):
        all_rod_save = all_rod_save+"_"+rodada[rod]
    out_dir = fig_output_dir+"/"+all_rod_save+"/CAMPOS_VEL"
    if not (os.path.isdir(out_dir)):
       os.system("mkdir -p "+out_dir)
    
    ax[0].set_extent([min_lon, max_lon, min_lat, max_lat])
    
    parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
    meridians = np.arange(int(min_lon)+2,int(max_lon)-1,lon_step_plot)
    
    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)
    
    if prof[k,1]<100:
       im = ax[0].contourf(lon_roms, lat_roms, locals()['vel_roms_mean'], np.arange(0.2,0.82,.020),cmap=cmap,vmin=0.2,vmax=0.8, \
                           linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    elif prof[k,1]==100:
       im = ax[0].contourf(lon_roms, lat_roms, locals()['vel_roms_mean'], np.arange(0.16,0.62,.020),cmap=cmap,vmin=0.16,vmax=0.6, \
                           linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    elif prof[k,1]==400:
       im = ax[0].contourf(lon_roms, lat_roms, locals()['vel_roms_mean'], np.arange(0.16,0.52,.020),cmap=cmap,vmin=0.16,vmax=0.5, \
                           linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    elif prof[k,1]>400:
       im = ax[0].contourf(lon_roms, lat_roms, locals()['vel_roms_mean'], np.arange(0.10,0.36,.010),cmap=cmap,vmin=0.10,vmax=0.35, \
                           linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    ax[0].add_feature(cfeature.GSHHSFeature(scale='high'))
    ax[0].add_feature(cfeature.LAND)
    ax[0].add_feature(cfeature.COASTLINE)
    if prof[k,1]<400:
       x,y = np.meshgrid(lon_roms[0:-1:12],lat_roms[0:-1:12])
       ax[0].quiver(x, y, u_roms_mean[0:-1:12,0:-1:12], v_roms_mean[0:-1:12,0:-1:12], transform=ccrs.PlateCarree(),scale=12.00, \
                    width=0.007,headaxislength=5, alpha=0.5)
       x,y = np.meshgrid(min_lon+1,max_lat-3.7)
       ax[0].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=12.00,width=0.007,headaxislength=5)
    elif prof[k,1]>=400:
       x,y = np.meshgrid(lon_roms[0:-1:12],lat_roms[0:-1:12])
       ax[0].quiver(x, y, u_roms_mean[0:-1:12,0:-1:12], v_roms_mean[0:-1:12,0:-1:12], transform=ccrs.PlateCarree(),scale=8.00, \
                    width=0.007,headaxislength=5, alpha=0.5)
       x,y = np.meshgrid(min_lon+1,max_lat-3.7)
       ax[0].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)
    
    states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
    ax[0].add_feature(states_provinces, edgecolor='gray')
    
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
    ax[0].text(min_lon+1, max_lat-4.8, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
    
    all_rod = 'ROMS'
    for rod in range(0,len(rodada),1):
    
        ####### PLOT MEAN MAP ###########
    
        cmap.set_under(under)
        cmap.set_over(over)
        cmap.set_bad(bad)
    
        ax[rod+1].set_extent([min_lon, max_lon, min_lat, max_lat])
        if prof[k,1]<100:
           im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['vel_hycom_mean_'+rodada[rod]], np.arange(0.2,0.82,.020),cmap=cmap, \
                                vmin=0.2,vmax=0.8,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
           cbar_ticks = np.arange(0.2,0.81,0.1)
        elif prof[k,1]==100:
           im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['vel_hycom_mean_'+rodada[rod]], np.arange(0.16,0.62,.020),cmap=cmap, \
                                   vmin=0.16,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
           cbar_ticks = np.insert(np.arange(0.2,0.61,0.1),[0],0.16)
        elif prof[k,1]==400:
           im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['vel_hycom_mean_'+rodada[rod]], np.arange(0.16,0.52,.020), \
                                   cmap=cmap,vmin=0.16,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
           cbar_ticks = np.insert(np.arange(0.2,0.51,0.1),[0],0.16)
           
        elif prof[k,1]>400:
           im = ax[rod+1].contourf(lon_roms, lat_roms, locals()['vel_hycom_mean_'+rodada[rod]], np.arange(0.10,0.36,.010), \
                                   cmap=cmap,vmin=0.10,vmax=0.35,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
           cbar_ticks = np.arange(0.1,0.40,0.05)

        ax[rod+1].add_feature(cfeature.GSHHSFeature(scale='high'))
        ax[rod+1].add_feature(cfeature.LAND)
        ax[rod+1].add_feature(cfeature.COASTLINE)
        ax[rod+1].add_feature(states_provinces, edgecolor='gray')
        if prof[k,1]<400:
           x,y = np.meshgrid(lon_roms[0:-1:12],lat_roms[0:-1:12])
           ax[rod+1].quiver(x,y,locals()['u_hycom_mean_'+rodada[rod]][0:-1:12,0:-1:12],locals()['v_hycom_mean_'+rodada[rod]][0:-1:12,0:-1:12], \
                            transform=ccrs.PlateCarree(),scale=12,width=0.007,headaxislength=5, alpha=0.5)
           x,y = np.meshgrid(min_lon+1,max_lat-3.7)
           ax[rod+1].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=12.00,width=0.007,headaxislength=5)
        if prof[k,1]>=400:
           x,y = np.meshgrid(lon_roms[0:-1:12],lat_roms[0:-1:12])
           ax[rod+1].quiver(x,y,locals()['u_hycom_mean_'+rodada[rod]][0:-1:12,0:-1:12],locals()['v_hycom_mean_'+rodada[rod]][0:-1:12,0:-1:12], \
                            transform=ccrs.PlateCarree(),scale=8,width=0.007,headaxislength=5, alpha=0.5)
           x,y = np.meshgrid(min_lon+1,max_lat-3.7)
           ax[rod+1].quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=8.00,width=0.007,headaxislength=5)

    
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
        ax[rod+1].text(min_lon+1, max_lat-3, legenda[rod], horizontalalignment='left', fontsize=10, fontweight='bold', \
                       transform=ccrs.PlateCarree())
        ax[rod+1].text(min_lon+1, max_lat-4.8, '0.5 m/s', horizontalalignment='left', fontsize=10, fontweight='bold', \
                       transform=ccrs.PlateCarree())
    
        all_rod = all_rod+"_"+rodada[rod]
        if rod==len(rodada)-1:
           fig.tight_layout()
           
           ## Add a colorbar axis at the bottom of the graph
           if len(rodada)==1:
              fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
              cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
              plt.suptitle('VEL MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.86)
           elif len(rodada)==2:
              fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
              cbar_ax = fig.add_axes([0.99, 0.28, 0.02, 0.42])
              plt.suptitle('VEL MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=0.75)
           elif len(rodada)==3:
              xmax = ax[0].get_position().xmax
              xmin = ax[1].get_position().xmin
              if round(min_lat,2) == -33.74 and max_lat == -12.32 and min_lon == -53.33 and max_lon == -32.46:
                 fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.25, hspace=0.02)
              elif round(min_lat,2) == -31.74 and max_lat == -14.32 and min_lon == -53.33 and max_lon == -34.46:
                 fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=-0.18, hspace=0.02)
              xmax = ax[rod+1].get_position().xmax
              ymax = ax[0].get_position().ymax
              ymin = ax[rod+1].get_position().ymin
              cbar_ax = fig.add_axes([xmax+0.01, ymin, 0.02, ymax-ymin])
              if len(months)==12:
                 plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m ('+data_inicial.strftime("%Y%m%d")+' - '+ \
                              data_final.strftime("%Y%m%d")+')', y=ymax+0.035)
              else:
                 plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m (MONTHS: '+str(months[0]).zfill(2)+'-'+ \
                              str(months[-1]).zfill(2)+'; YEARS: '+str(data_inicial.year)+'-'+str(data_final.year)+')', y=ymax+0.035)
           cbar=fig.colorbar(im, ticks = cbar_ticks, cax=cbar_ax,orientation='vertical')
           cbar.ax.set_ylabel('m/s')
    
           if len(months)==12:
              fig.savefig(out_dir+'/VEL_MEAN_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                          '_'+dominio_str+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False, \
                          bbox_inches='tight',pad_inches=0.05)
           else:
              fig.savefig(out_dir+'/VEL_MEAN_'+all_rod+'_months_'+str(months[0]).zfill(2)+'-'+str(months[-1]).zfill(2)+'_YEARS_'+ \
                          str(data_inicial.year)+'-'+str(data_final.year)+'_'+dominio_str+'_DEPTH_'+str(prof[k,0])+"-"+ \
                          str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
    
    del all_rod,fig,ax

exit()

