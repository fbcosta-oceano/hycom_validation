import remo
import datetime
import time
import sys
import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import pandas as pd
import os
from mpl_toolkits.basemap import Basemap
from calendar import monthrange
jet = matplotlib.cm.get_cmap('jet')

max_prof_rmsd_time = 1000
max_prof_rmsd_profile = 1000

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

input_dir = '/home/filipe/resultados'
fig_output_dir = '/home/filipe/resultados/figs'

month_label_interval = 1

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

months = [1,2,3,4,5,6,7,8,9,10,11,12]
#months = [7,8,9,10,11,12,1,2,3,4,5,6]

if '004j' in rodada[0]:
   min_lat = -33.74+2
   max_lat = -12.32-2
   min_lon = -53.33
   max_lon = -32.46-2

   lat_step_plot = 4
   lon_step_plot = 5
   line_color = ['red','blue','green','dodgerblue','yellow']
   resolution_plot = 150

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

all_rod = rodada[0]
for rod in range(1,len(rodada),1):
    all_rod = all_rod+"_"+rodada[rod]
out_dir = fig_output_dir+"/"+all_rod+"/TS_ROMS_HYCOM"
if not (os.path.isdir(out_dir)):
   os.system("mkdir -p "+out_dir)

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    counter = 0
    while current_data<=data_final:
        while current_data.month != months[0]:
           days = monthrange(current_data.year,current_data.month)[1]
           current_data = current_data + datetime.timedelta(days=days)
        for m in list(months):
            days = monthrange(current_data.year,current_data.month)[1]
            filename = input_dir+"/"+rodada[rod]+"/RMSD_ROMS_HYCOM_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                       str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+str(days).zfill(2)+"_"+dominio_str+".mat"
            if not (os.path.isfile(filename)):
               print(filename, ' NAO ENCONTRADO')
               exit()

            mat_contents = scipy.io.loadmat(filename)
            if not 'levels' in locals():
               levels = np.squeeze(mat_contents['levels'])
            if not 'mean_profile_temp_'+rodada[rod] in locals():
               locals()['mean_profile_temp_'+rodada[rod]] = mat_contents['mean_profile_temp_'+rodada[rod]]
               locals()['mean_profile_salt_'+rodada[rod]] = mat_contents['mean_profile_salt_'+rodada[rod]]
               locals()['mean_time_temp_'+rodada[rod]] = mat_contents['mean_time_temp_'+rodada[rod]]
               locals()['mean_time_salt_'+rodada[rod]] = mat_contents['mean_time_salt_'+rodada[rod]]
               locals()['rmsd_profile_temp_'+rodada[rod]] = mat_contents['rmsd_profile_temp_'+rodada[rod]]
               locals()['rmsd_profile_salt_'+rodada[rod]] = mat_contents['rmsd_profile_salt_'+rodada[rod]]
               locals()['rmsd_time_temp_'+rodada[rod]] = mat_contents['rmsd_time_temp_'+rodada[rod]]
               locals()['rmsd_time_salt_'+rodada[rod]] = mat_contents['rmsd_time_salt_'+rodada[rod]]
               if rod==0:
                  mean_profile_roms_temp = mat_contents['mean_profile_roms_temp']
                  mean_profile_roms_salt = mat_contents['mean_profile_roms_salt']
                  mean_time_roms_temp = mat_contents['mean_time_roms_temp']
                  mean_time_roms_salt = mat_contents['mean_time_roms_salt']
            else:
               locals()['mean_profile_temp_'+rodada[rod]] = locals()['mean_profile_temp_'+rodada[rod]] + mat_contents['mean_profile_temp_'+rodada[rod]]
               locals()['mean_profile_salt_'+rodada[rod]] = locals()['mean_profile_salt_'+rodada[rod]] + mat_contents['mean_profile_salt_'+rodada[rod]]
               locals()['mean_time_temp_'+rodada[rod]] = np.ma.concatenate((locals()['mean_time_temp_'+rodada[rod]],mat_contents['mean_time_temp_'+rodada[rod]]),axis=0)
               locals()['mean_time_salt_'+rodada[rod]] = np.ma.concatenate((locals()['mean_time_salt_'+rodada[rod]],mat_contents['mean_time_salt_'+rodada[rod]]),axis=0)
               locals()['rmsd_profile_temp_'+rodada[rod]] = locals()['rmsd_profile_temp_'+rodada[rod]] + mat_contents['rmsd_profile_temp_'+rodada[rod]]
               locals()['rmsd_profile_salt_'+rodada[rod]] = locals()['rmsd_profile_salt_'+rodada[rod]] + mat_contents['rmsd_profile_salt_'+rodada[rod]]
               locals()['rmsd_time_temp_'+rodada[rod]] = np.ma.concatenate((locals()['rmsd_time_temp_'+rodada[rod]],mat_contents['rmsd_time_temp_'+rodada[rod]]),axis=0)
               locals()['rmsd_time_salt_'+rodada[rod]] = np.ma.concatenate((locals()['rmsd_time_salt_'+rodada[rod]],mat_contents['rmsd_time_salt_'+rodada[rod]]),axis=0)
               if rod==0:
                  mean_profile_roms_temp = mean_profile_roms_temp + mat_contents['mean_profile_roms_temp']
                  mean_profile_roms_salt = mean_profile_roms_salt + mat_contents['mean_profile_roms_salt']
                  mean_time_roms_temp = np.ma.concatenate((mean_time_roms_temp,mat_contents['mean_time_roms_temp']),axis=0)
                  mean_time_roms_salt = np.ma.concatenate((mean_time_roms_salt,mat_contents['mean_time_roms_salt']),axis=0)

            counter = counter+1
            current_data = current_data + datetime.timedelta(days=days)
            if current_data>data_final:
               break

    if rod==0:
       cond = abs(mean_profile_roms_temp) > 99999
       mean_profile_roms_temp = np.ma.masked_where(cond,mean_profile_roms_temp)/counter
       cond = abs(mean_profile_roms_salt) > 99999
       mean_profile_roms_salt = np.ma.masked_where(cond,mean_profile_roms_salt)/counter
       cond = abs(mean_time_roms_temp) > 99999
       mean_time_roms_temp = np.ma.masked_where(cond,mean_time_roms_temp)
       cond = abs(mean_time_roms_salt) > 99999
       mean_time_roms_temp = np.ma.masked_where(cond,mean_time_roms_salt)


    cond = abs(locals()['mean_profile_temp_'+rodada[rod]]) > 99999
    locals()['mean_profile_temp_'+rodada[rod]] = np.ma.masked_where(cond,locals()['mean_profile_temp_'+rodada[rod]])/counter
    cond = abs(locals()['mean_profile_salt_'+rodada[rod]]) > 99999
    locals()['mean_profile_salt_'+rodada[rod]] = np.ma.masked_where(cond,locals()['mean_profile_salt_'+rodada[rod]])/counter
    cond = abs(locals()['mean_time_temp_'+rodada[rod]]) > 99999
    locals()['mean_time_temp_'+rodada[rod]] = np.ma.masked_where(cond,locals()['mean_time_temp_'+rodada[rod]])
    cond = abs(locals()['mean_time_salt_'+rodada[rod]]) > 99999
    locals()['mean_time_salt_'+rodada[rod]] = np.ma.masked_where(cond,locals()['mean_time_salt_'+rodada[rod]])

    cond = abs(locals()['rmsd_profile_temp_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_temp_'+rodada[rod]] = np.sqrt(np.ma.masked_where(cond,locals()['rmsd_profile_temp_'+rodada[rod]])/counter)
    cond = abs(locals()['rmsd_profile_salt_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_salt_'+rodada[rod]] = np.sqrt(np.ma.masked_where(cond,locals()['rmsd_profile_salt_'+rodada[rod]])/counter)
    cond = abs(locals()['rmsd_time_temp_'+rodada[rod]]) > 99999
    locals()['rmsd_time_temp_'+rodada[rod]] = np.sqrt(np.ma.masked_where(cond,locals()['rmsd_time_temp_'+rodada[rod]]))
    cond = abs(locals()['rmsd_time_salt_'+rodada[rod]]) > 99999
    locals()['rmsd_time_salt_'+rodada[rod]] = np.sqrt(np.ma.masked_where(cond,locals()['rmsd_time_salt_'+rodada[rod]]))



plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
for rod in range(0,len(rodada),1):
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    cond = np.array(levels)<=max_prof_rmsd_time
    temp = locals()['rmsd_time_temp_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq='d'))
    ax = data_plot.plot(label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0.8,2.4,.2))
       ax.set_yticks(np.arange(0.8,2.3,.1),minor=True)
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       #ax.grid(False, which='minor', axis = 'y')
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.1f'))
       plt.ylabel('\N{DEGREE SIGN}C')
       plt.ylim(0.8,2.2)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('TEMPERATURE RMSD (0-'+str(max_prof_rmsd_time)+'m; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_TEMP_TIME_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()
    del date_vec, data_plot

plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
for rod in range(0,len(rodada),1):
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    cond = np.array(levels)<=max_prof_rmsd_time
    temp = locals()['rmsd_time_salt_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq='d'))
    ax = data_plot.plot(label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0.18,0.38,.02))
       ax.set_yticks(np.arange(0.18,0.37,.01),minor=True)
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       #ax.grid(False, which='minor', axis = 'y')
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.2f'))
       plt.ylabel('PSU')
       plt.ylim(0.18,0.36)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('SALINITY RMSD (0-'+str(max_prof_rmsd_time)+'m; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_SALT_TIME_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()
    del date_vec, data_plot

plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['rmsd_profile_temp_'+rodada[rod]])
    temp = temp[cond]
    ax = plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('\N{DEGREE SIGN}C')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.title('TEMPERATURE RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_TEMP_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()

plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['rmsd_profile_salt_'+rodada[rod]])
    temp = temp[cond]
    ax = plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('PSU')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.title('SALINITY RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_SALT_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()

plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
cond = np.array(levels)<=max_prof_rmsd_profile
temp = np.squeeze(mean_profile_roms_temp)
temp = temp[cond]
ax = plt.plot(temp, levels[cond], label='NATURE RUN', color='black')
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['mean_profile_temp_'+rodada[rod]])
    temp = temp[cond]
    ax = plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('\N{DEGREE SIGN}C')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.title('TEMPERATURE ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_MEAN_TEMP_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()

plt.figure(figsize=(16, 8), dpi=150)
plt.rcParams.update({'font.size': 20})
cond = np.array(levels)<=max_prof_rmsd_profile
temp = np.squeeze(mean_profile_roms_salt)
temp = temp[cond]
ax = plt.plot(temp, levels[cond], label='NATURE RUN', color='black')
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['mean_profile_salt_'+rodada[rod]])
    temp = temp[cond]
    ax = plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('PSU')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.title('SALINITY ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_MEAN_SALT_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+"_0-"+str(max_prof_rmsd_time)+"m.png",dpi=150)
       plt.close()

