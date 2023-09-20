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
jet = matplotlib.cm.get_cmap('jet')

max_prof_rmsd_time = 1600
max_prof_rmsd_profile = 1600

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

input_dir = '/home/filipe/resultados'
fig_output_dir = '/home/filipe/resultados/figs'

month_label_interval = 1

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

if '004j' in rodada[0]:
   min_lat = -33.74
   max_lat = -12.32
   min_lon = -53.33
   max_lon = -32.46

   lat_step_plot = 4
   lon_step_plot = 5
   line_color = ['red','blue']
   resolution_plot = 150

all_rod = rodada[0]
for rod in range(1,len(rodada),1):
    all_rod = all_rod+"_"+rodada[rod]
out_dir = fig_output_dir+"/"+all_rod+"/TS_ROMS_HYCOM"
if not (os.path.isdir(out_dir)):
   os.system("mkdir -p "+out_dir)

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
filename = "RMSD_ROMS_HYCOM_ANABACK_ARGOPOINTS_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str+".mat"

for rod in range(0,len(rodada),1):
    mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[rod]+"/"+filename)
    #locals()['rmsd_profile_temp_'+rodada[rod]] = mat_contents['rmsd_profile_temp_'+rodada[rod]]
    #locals()['rmsd_profile_saln_'+rodada[rod]] = mat_contents['rmsd_profile_saln_'+rodada[rod]]
    #locals()['rmsd_time_temp_'+rodada[rod]] = mat_contents['rmsd_time_temp_'+rodada[rod]]
    #locals()['rmsd_time_saln_'+rodada[rod]] = mat_contents['rmsd_time_saln_'+rodada[rod]]
    q_buoys = mat_contents['q_buoys_tot']
    levels = np.squeeze(mat_contents['levels'])

    cond = abs(mat_contents['rmsd_profile_temp_back_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_temp_back_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_temp_back_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_profile_salt_back_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_salt_back_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_salt_back_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_time_temp_back_'+rodada[rod]]) > 99999
    locals()['rmsd_time_temp_back_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_temp_back_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_time_salt_back_'+rodada[rod]]) > 99999
    locals()['rmsd_time_salt_back_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_salt_back_'+rodada[rod]])

    cond = abs(mat_contents['rmsd_profile_temp_ana_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_temp_ana_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_temp_ana_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_profile_salt_ana_'+rodada[rod]]) > 99999
    locals()['rmsd_profile_salt_ana_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_salt_ana_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_time_temp_ana_'+rodada[rod]]) > 99999
    locals()['rmsd_time_temp_ana_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_temp_ana_'+rodada[rod]])
    cond = abs(mat_contents['rmsd_time_salt_ana_'+rodada[rod]]) > 99999
    locals()['rmsd_time_salt_ana_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_salt_ana_'+rodada[rod]])

plt.figure(figsize=(16, 8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    cond = np.array(levels)<=max_prof_rmsd_time
    temp = locals()['rmsd_time_temp_back_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq=str(int(inc_tempo))+'d'))
    ax = data_plot.plot(label=legenda[rod]+'_BACK', color=line_color[rod], linestyle='dashed')

    temp = locals()['rmsd_time_temp_ana_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq=str(int(inc_tempo))+'d'))
    ax = data_plot.plot(label=legenda[rod]+'_ANA', color=line_color[rod], linestyle='solid')
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0.0,2.5,.2))
       ax.set_yticks(np.arange(0.0,2.4,.1),minor=True)
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.1f'))
       plt.ylabel('\N{DEGREE SIGN}C')
       plt.ylim(0.0,2.3)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('TEMPERATURE RMSD (0-'+str(max_prof_rmsd_time)+'m; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_ANABACK_TEMP_TIME_"+all_rod+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
    del date_vec, data_plot

plt.figure(figsize=(16, 8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['rmsd_profile_temp_back_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+'_BACK', color=line_color[rod],linestyle='dashed')

    temp = np.squeeze(locals()['rmsd_profile_temp_ana_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+'_ANA', color=line_color[rod],linestyle='solid')

    temp = np.squeeze(locals()['rmsd_profile_temp_ana_'+rodada[rod]]-locals()['rmsd_profile_temp_back_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+': ANA-BACK', color=line_color[rod],linestyle='dashdot')
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('\N{DEGREE SIGN}C')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.axvline(x = 0, color = 'k', linestyle='dashed')
       plt.title('TEMPERATURE RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_ANABACK_TEMP_PROFILE_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()

plt.figure(figsize=(16, 8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):
    date_vec = pd.date_range(data_inicial,data_final,freq='d')
    cond = np.array(levels)<=max_prof_rmsd_time
    temp = locals()['rmsd_time_salt_back_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq=str(int(inc_tempo))+'d'))
    ax = data_plot.plot(label=legenda[rod]+'_BACK', color=line_color[rod], linestyle='dashed')

    temp = locals()['rmsd_time_salt_ana_'+rodada[rod]]
    temp = np.mean(temp[:,cond],axis=1)
    data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq=str(int(inc_tempo))+'d'))
    ax = data_plot.plot(label=legenda[rod]+'_ANA', color=line_color[rod], linestyle='solid')
    if rod==len(rodada)-1:
       ax.set_yticks(np.arange(0.0,0.5,.1))
       ax.set_yticks(np.arange(0.0,0.45,.05),minor=True)
       ax.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=15, interval=month_label_interval))
       ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=np.arange(5,31,5)))
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%Y'))
       ax.xaxis.set_minor_formatter(plt.NullFormatter())

       ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%1.2f'))
       plt.ylabel('PSU')
       plt.ylim(0.0,0.4)
       for label in ax.get_xticklabels():
           label.set_rotation(30)
           label.set_horizontalalignment('right')

       plt.legend()
       plt.title('SALINITY RMSD (0-'+str(max_prof_rmsd_time)+'m; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_ANABACK_SALT_TIME_"+all_rod+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
    del date_vec, data_plot

plt.figure(figsize=(16, 8), dpi=resolution_plot)
plt.rcParams.update({'font.size': 10})
for rod in range(0,len(rodada),1):
    cond = np.array(levels)<=max_prof_rmsd_profile
    temp = np.squeeze(locals()['rmsd_profile_salt_back_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+'_BACK', color=line_color[rod],linestyle='dashed')

    temp = np.squeeze(locals()['rmsd_profile_salt_ana_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+'_ANA', color=line_color[rod],linestyle='solid')

    temp = np.squeeze(locals()['rmsd_profile_salt_ana_'+rodada[rod]]-locals()['rmsd_profile_salt_back_'+rodada[rod]])
    temp = temp[cond]
    plt.plot(temp, levels[cond], label=legenda[rod]+': ANA-BACK', color=line_color[rod],linestyle='dashdot')
    if rod==len(rodada)-1:
       plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
       #plt.yticks(np.arange(0,1600,100),minor=True)
       plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
       ax = plt.gca()
       ax.invert_yaxis()
       plt.xlabel('PSU')
       plt.ylabel('DEPTH')
       plt.legend()
       plt.axvline(x = 0, color = 'k', linestyle='dashed')
       plt.title('SALINITY RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
       plt.savefig(out_dir+"/ROMS_HYCOM_RMSD_ANABACK_SALT_PROFILE_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                  "_"+dominio_str+".png",dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
       plt.close()
