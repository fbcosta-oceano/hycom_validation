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
months = [1,2,3,4,5,6,7,8]
prof = np.array([[0,0], [200,200], [400,400], [800,800], [2000,2000], [4000,4000]])
#prof = np.array([[0,0]])

month_label_interval = 1
line_color = ['red','blue','green','yellow','cyan','black']
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
       resolution_plot = 300

    elif '008d' in rodada[0]:
       min_lat = -79.54
       max_lat =  50.27
       min_lon = -98.00
       max_lon =  45.00

       lat_step_plot = 3
       lon_step_plot = 3
       resolution_plot = 300
    
    elif '008i' in rodada[0]:
       min_lat = -47.93
       max_lat =  10.11
       min_lon = -70.00
       max_lon = -17.75

       lat_step_plot = 8
       lon_step_plot = 8
       font_size = 12
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
           EC_hycom_mean_map = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
           EC_hycom_mean_map = EC_hycom_mean_map.reshape(len(lat_hycom),len(lon_hycom))
    
        while current_data<=data_final:
            while current_data.month != months[0]:
               days = monthrange(current_data.year,current_data.month)[1]
               current_data = current_data + datetime.timedelta(days=days)
            for m in list(months):
                days = monthrange(current_data.year,current_data.month)[1]
                filename1 = input_dir+"/"+rodada[rod]+"/VEL_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                            str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+str(days).zfill(2)+"_"+dominio_str+"_"+ \
                            str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                filename2 = input_dir+"/"+rodada[rod]+"/EC_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                            str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+str(days).zfill(2)+"_"+dominio_str+"_"+ \
                            str(prof[k,0])+"-"+str(prof[k,1])+"m"+".mat"
                if not (os.path.isfile(filename1)):
                   print(filename1, ' NAO ENCONTRADO')
                   exit()
                mat_contents1 = scipy.io.loadmat(filename1)
                mat_contents2 = scipy.io.loadmat(filename2)
                if not 'lat_hycom' in locals():
                   locals()['lat_hycom'] = mat_contents2['lat_hycom'].flatten()
                   locals()['lon_hycom'] = mat_contents2['lon_hycom'].flatten()
    
                if y == data_inicial.year and m==(months[0]):
                   u_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   u_hycom_mean = u_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
                   v_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   v_hycom_mean = v_hycom_mean.reshape(len(lat_hycom),len(lon_hycom))
                   vel_hycom_mean = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   EC_hycom_mean_map = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
                   EC_hycom_mean_map = EC_hycom_mean_map.reshape(len(lat_hycom),len(lon_hycom))
    
                u_hycom = mat_contents1['u_hycom_mean_'+rodada[rod]]
                v_hycom = mat_contents1['v_hycom_mean_'+rodada[rod]]
                EC_hycom_map = mat_contents2['EC_hycom_mean_map_'+rodada[rod]]
                EC_hycom_time = mat_contents2['EC_hycom_mean_time_'+rodada[rod]]
                cond = abs(u_hycom) > 99999
                u_hycom = np.ma.masked_where(cond,u_hycom)
                cond = abs(v_hycom) > 99999
                v_hycom = np.ma.masked_where(cond,v_hycom)
                cond = abs(EC_hycom_map) > 99999
                EC_hycom_map = np.ma.masked_where(cond,EC_hycom_map)
                cond = abs(EC_hycom_time) > 99999
                EC_hycom_time = np.ma.masked_where(cond,EC_hycom_time)
                u_hycom_mean = u_hycom_mean + u_hycom
                v_hycom_mean = v_hycom_mean + v_hycom
                EC_hycom_mean_map = EC_hycom_mean_map + EC_hycom_map
                if 'EC_hycom_time_all' in locals():
                    EC_hycom_time_all = np.ma.concatenate((EC_hycom_time_all,EC_hycom_time),axis=1)
                else:
                    EC_hycom_time_all = EC_hycom_time
    
                del u_hycom, v_hycom, EC_hycom_map, EC_hycom_time, mat_contents1, mat_contents2, cond
                counter = counter+1
                current_data = current_data + datetime.timedelta(days=days)
    
        locals()['u_hycom_mean_'+rodada[rod]] = u_hycom_mean/counter
        locals()['v_hycom_mean_'+rodada[rod]] = v_hycom_mean/counter
        locals()['EC_hycom_mean_map_'+rodada[rod]] = EC_hycom_mean_map/counter
        locals()['EC_hycom_mean_time_'+str(prof[k,0])+"-"+str(prof[k,1])+'_'+rodada[rod]] = EC_hycom_time_all
        del EC_hycom_mean_map, EC_hycom_time_all, counter, y, m, days, current_data
    
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
    out_dir = fig_output_dir+"/"+all_rod_save+"/EC"
    if not (os.path.isdir(out_dir)):
       os.system("mkdir -p "+out_dir)
    
    parallels = np.arange(int(min_lat),int(max_lat)+0.1,lat_step_plot)
    meridians = np.arange(int(min_lon),int(max_lon)+0.1,lon_step_plot)
    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    #for rod in range(0,len(rodada),1):
    #
    #    ####### PLOT MEAN MAP ###########
    #
    #    EC_hycom_plot = locals()['EC_hycom_mean_map_'+rodada[rod]]
    #    u_hycom_plot = locals()['u_hycom_mean_'+rodada[rod]]
    #    v_hycom_plot = locals()['v_hycom_mean_'+rodada[rod]]

    #    if len(rodada)==1:
    #       ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    #       teste = ax
    #       y_title=1.020
    #       #cbar_ax = [0.855, 0.014, 0.02, 0.977]
    #       #cbar_ax = [0.800, 0.012, 0.02, 0.977]
    #       #cbar_ax = [0.780, 0.010, 0.02, 0.977]
    #       cbar_ax = [0.790, 0.010, 0.02, 0.977]
    #    else:
    #       teste.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    #       teste = teste
    #       y_title=1.020
    #       cbar_ax = [0.855, 0.014, 0.02, 0.977]

    #    if prof[k,1]<100:
    #       im = teste.contourf(lon_hycom, lat_hycom, EC_hycom_plot, np.arange(0.08,0.52,.02),cmap=cmap, \
    #                            vmin=0.08,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #       cbar_ticks = np.insert(np.arange(0.1,0.51,0.1),[0],0.08)
    #    elif prof[k,1]==100:
    #       im = teste.contourf(lon_hycom, lat_hycom, EC_hycom_plot, np.arange(0.15,0.625,.025),cmap=cmap, \
    #                               vmin=0.15,vmax=0.6,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #       cbar_ticks = np.insert(np.arange(0.2,0.61,0.1),[0],0.15)
    #    elif prof[k,1]==400:
    #       im = teste.contourf(lon_hycom, lat_hycom, EC_hycom_plot, np.arange(0.10,0.525,.025), \
    #                               cmap=cmap,vmin=0.10,vmax=0.5,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #       cbar_ticks = np.arange(0.1,0.51,0.1)
    #       
    #    elif prof[k,1]>400:
    #       im = teste.contourf(lon_hycom, lat_hycom, EC_hycom_plot, np.arange(0.10,0.41,.010), \
    #                               cmap=cmap,vmin=0.10,vmax=0.40,linestyles='dashed',transform=ccrs.PlateCarree(),extend='both')
    #       cbar_ticks = np.arange(0.1,0.40,0.05)

    #    if prof[k,1]<400:
    #       x,y = np.meshgrid(lon_hycom[0:-1:15],lat_hycom[0:-1:15])
    #       teste.quiver(x,y,u_hycom_plot[0:-1:15,0:-1:15],v_hycom_plot[0:-1:15,0:-1:15], \
    #                        transform=ccrs.PlateCarree(),scale=28,width=0.006,headaxislength=3.5,headwidth=3, alpha=0.45)
    #       x,y = np.meshgrid(-68,-9)
    #       teste.quiver(x, y, 1.0, 0, transform=ccrs.PlateCarree(),scale=28.00,width=0.006,headaxislength=3.5,headwidth=2)

    #    elif prof[k,1]>=400:
    #       x,y = np.meshgrid(lon_hycom[0:-1:15],lat_hycom[0:-1:15])
    #       teste.quiver(x,y,u_hycom_plot[0:-1:15,0:-1:15],v_hycom_plot[0:-1:15,0:-1:15], \
    #                        transform=ccrs.PlateCarree(),scale=5.0,width=0.007,headaxislength=5, alpha=0.5)
    #       x,y = np.meshgrid(-68,-9)
    #       teste.quiver(x, y, 0.5, 0, transform=ccrs.PlateCarree(),scale=5.00,width=0.007,headaxislength=5)

    #    teste.text(-68, -8, legenda[rod], horizontalalignment='left', fontsize=font_size, fontweight='bold', \
    #                   transform=ccrs.PlateCarree())
    #    teste.text(-68, -11.5, '1.0 $\mathregular{m^2/s^2}$', horizontalalignment='left', fontsize=font_size, fontweight='bold', \
    #                      transform=ccrs.PlateCarree())

    #    teste.add_feature(cfeature.LAND)
    #    teste.add_feature(cfeature.COASTLINE, linewidth=0.5)
    #    states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines', \
    #                       scale='10m',facecolor='none')
    #    teste.add_feature(states_provinces, edgecolor='gray', linewidth=0.2)
    #
    #    gl = teste.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    #    gl.xlabels_top = False
    #    gl.right_labels = False
    #    if len(rodada)<=2:
    #        if rod==0: 
    #           gl.left_labels = True 
    #        elif rod==1: 
    #           gl.left_labels = False
    #    elif len(rodada)==3:
    #       if rod==0:
    #          gl.left_labels = False
    #          gl.bottom_labels = False
    #       elif rod==2:
    #          gl.left_labels = False
    #    gl.xlocator = mticker.FixedLocator(meridians)
    #    gl.ylocator = mticker.FixedLocator(parallels)
    #    gl.ylabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    #    gl.xlabel_style = {'size': font_size, 'color': 'black', 'weight': 'bold'}
    #
    #    teste.tick_params(axis='both', labelsize=10)
    #
    #    if rod==0:
    #        all_rod = rodada[rod]
    #    else:
    #        all_rod = all_rod+"_"+rodada[rod]
    #    if rod==len(rodada)-1:
    #       fig.tight_layout()
    #       
    #       ## Add a colorbar axis at the bottom of the graph
    #       if len(rodada)==1:
    #          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
    #          #cbar_ax = fig.add_axes([0.99, 0.2, 0.02, 0.6])
    #          cbar_ax = fig.add_axes(cbar_ax)
    #          #plt.suptitle('VEL MEAN ('+data_inicial.strftime("%Y%m%d")+' - '+data_final.strftime("%Y%m%d")+')', y=y_title)
    #          if len(months)==12:
    #             plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m ('+data_inicial.strftime("%Y%m%d")+' - '+ \
    #                          data_final.strftime("%Y%m%d")+')', y=y_title)
    #          else:
    #             plt.suptitle('VEL MEAN DEPTH '+str(prof[k,0])+"-"+str(prof[k,1])+'m (MONTHS: '+str(months[0]).zfill(2)+'-'+ \
    #                          str(months[-1]).zfill(2)+'; YEARS: '+str(data_inicial.year)+'-'+str(data_final.year)+')', y=y_title)
    #       elif len(rodada)==2:
    #          fig.subplots_adjust(bottom=0.0, top=0.99, left=0.0, wspace=0.02)
    #          cbar_ax = fig.add_axes([0.99, 0.20, 0.02, 0.58])
    #          plt.suptitle('VEL MEAN '+str(prof[k,0])+"-"+str(prof[k,1])+'m ('+data_inicial.strftime("%Y%m%d")+' - '+ \
    #                       data_final.strftime("%Y%m%d")+')', y=y_title)

    #       cbar=fig.colorbar(im, ticks = cbar_ticks, cax=cbar_ax,orientation='vertical')
    #       cbar.ax.set_ylabel('$\mathregular{m^2/s^2}$')
    #
    #       if len(months)==12:
    #          fig.savefig(out_dir+'/EC_MEAN_MAP_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
    #                      '_'+dominio_str+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False, \
    #                      bbox_inches='tight',pad_inches=0.05)
    #       else:
    #          fig.savefig(out_dir+'/EC_MEAN_MAP_'+all_rod+'_months_'+str(months[0]).zfill(2)+'-'+str(months[-1]).zfill(2)+'_YEARS_'+ \
    #                      str(data_inicial.year)+'-'+str(data_final.year)+'_'+dominio_str+'_DEPTH_'+str(prof[k,0])+"-"+ \
    #                      str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
    #
    #del all_rod,fig,ax

    ####### PLOT EC TIME ###########
    ####### ALL EXPERIMENTS; ONE DEPTH ###########

    if len(months) == 12:
       plt.figure(figsize=(6.4, 4.8), dpi=resolution_plot)
       plt.rcParams.update({'font.size': 10})
       for rod in range(0,len(rodada),1):
       
           if rod==0:
              all_rod = rodada[rod]
           else:
              all_rod = all_rod+"_"+rodada[rod]
           date_vec = pd.date_range(data_inicial,data_final,freq='d')
           data_plot = pd.Series(data=np.squeeze(locals()['EC_hycom_mean_time_'+str(prof[k,0])+"-"+str(prof[k,1])+'_'+rodada[rod]]), \
                                 index = pd.date_range(data_inicial,data_final,freq='d'))
           ax = data_plot.plot(label=legenda[rod], color=line_color[rod])
           if rod==len(rodada)-1:
              ax.set_yticks(np.arange(0,0.12,.02)) 
              ax.set_yticks(np.arange(0,0.11,.01),minor=True) 
              ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
              #ax.grid(False, which='minor', axis = 'y')
              ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
              ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
              ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
              ax.xaxis.set_minor_formatter(plt.NullFormatter())
       
              ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.2f'))
              plt.ylabel('$\mathregular{m^2/s^2}$')
              plt.ylim(0,0.10)
              for label in ax.get_xticklabels():
                  label.set_rotation(30)
                  label.set_horizontalalignment('right')
       
              plt.legend()
              plt.title('EC '+str(prof[k,0])+"-"+str(prof[k,1])+'m')
              plt.savefig(out_dir+'/EC_MEAN_TIME_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                          '_'+dominio_str+'_DEPTH_'+str(prof[k,0])+"-"+str(prof[k,1])+'m.png',dpi=resolution_plot,transparent=False, \
                          bbox_inches='tight',pad_inches=0.05)
           del date_vec, data_plot

if len(months) == 12:
   for rod in range(0,len(rodada),1):
       all_rod = rodada[rod]
       plt.figure(figsize=(6.4, 4.8), dpi=resolution_plot)
       plt.rcParams.update({'font.size': 10})

       for k in range(0,len(prof[:,0]),1):
           date_vec = pd.date_range(data_inicial,data_final,freq='d')
           data_plot = pd.Series(data=np.squeeze(locals()['EC_hycom_mean_time_'+str(prof[k,0])+"-"+str(prof[k,1])+'_'+rodada[rod]]), \
                                 index = pd.date_range(data_inicial,data_final,freq='d'))
           if prof[k,0] == prof[k,1]:
              ax = data_plot.plot(label=legenda[rod]+' '+str(prof[k,1])+'m', color=line_color[k])
           else:
              ax = data_plot.plot(label=legenda[rod]+' '+str(prof[k,0])+"-"+str(prof[k,1])+'m', color=line_color[k])
           if k==len(prof[:,0])-1:
              ax.set_yticks(np.arange(0,0.12,.02)) 
              ax.set_yticks(np.arange(0,0.11,.01),minor=True) 
              ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
              #ax.grid(False, which='minor', axis = 'y')
              ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
              ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
              ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
              ax.xaxis.set_minor_formatter(plt.NullFormatter())
       
              ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.2f'))
              plt.ylabel('$\mathregular{m^2/s^2}$')
              plt.ylim(0,0.10)
              for label in ax.get_xticklabels():
                  label.set_rotation(30)
                  label.set_horizontalalignment('right')
       
              plt.legend()
              plt.title('EC')
              plt.savefig(out_dir+'/EC_MEAN_TIME_'+all_rod+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                          '_'+dominio_str+'_ALL_DEPTHS.png',dpi=resolution_plot,transparent=False, \
                          bbox_inches='tight',pad_inches=0.05)
           del date_vec, data_plot

exit()

