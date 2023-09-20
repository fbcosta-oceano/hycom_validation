import remo
import datetime
import time
import os
import sys
import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
jet = matplotlib.cm.get_cmap('jet')

max_prof_rmsd_time = 1500
max_prof_rmsd_profile = 1500
whole_domain = True

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

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

if '004j' in rodada[0]:
   min_lat = -33.74
   max_lat = -12.32
   min_lon = -53.33
   max_lon = -32.46

   line_color = ['red','blue','green','yellow']
   resolution_plot = 150

elif '008d' in rodada[0]:
   subs = np.array([[-98,  45, -79.54, 50.27], \
                    [-68, -18, -45,  10], \
                    [-60, -20, -36,   7], \
                    [-53, -33, -33, -13], \
                    [-81, -56,  26,  41], \
                    [-51, -31,  26,  46], \
                    [-26,  -9,  26,  46], \
                    [-98, -81,  18,  30], \
                    [-61, -36,  -5,  18], \
                    [-31, -16,  -5,  18], \
                    [-14,  10, -12,   6], \
                    [-51, -18, -30, -10], \
                    [-15,  14, -30, -15], \
                    [-67, -40, -49, -34], \
                    [-36,  -5, -49, -34], \
                    [  0,  25, -49, -34], \
                    [ 28,  43, -49, -19]]) 

   line_color = ['red','blue','green','yellow']
   resolution_plot = 150

elif '008i' in rodada[0]:
   subs = np.array([[-70, -17.75, -47.93, 10.11], \
                    [-60, -20, -36,   7], \
                    [-53, -33, -33, -13], \
                    [-51, -18, -30, -10], \
                    [-67, -40, -47, -34], \
                    [-36,  -5, -47, -34]])

   line_color = ['red','blue','green','yellow']
   resolution_plot = 150

for region in range(0,subs.shape[0]):
    if(not(whole_domain) and (region==0)):
        continue
    if subs[region,2]<0:
       s1 = "S"
    else:
       s1 = "N"
    if subs[region,3]<0:
       s2 = "S"
    else:
       s2 = "N"
    if subs[region,0]<0:
       s3 = "W"
    else:
       s3 = "E"
    if subs[region,1]<0:
       s4 = "W"
    else:
       s4 = "E"
    dominio_str = f"{abs(subs[region,2]):.2f}"+s1+"-"+f"{abs(subs[region,3]):.2f}"+s2+"_" + f"{abs(subs[region,0]):.2f}"+s3+"-"+f"{abs(subs[region,1]):.2f}"+s4
    filename = "RMSD_ARGO_HYCOM_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str+".mat"
    
    for rod in range(0,len(rodada),1):
        mat_contents = scipy.io.loadmat(input_dir+"/"+rodada[rod]+"/"+filename)
        locals()['rmsd_profile_temp_'+rodada[rod]] = mat_contents['rmsd_profile_temp_'+rodada[rod]]
        locals()['rmsd_profile_saln_'+rodada[rod]] = mat_contents['rmsd_profile_saln_'+rodada[rod]]
        locals()['rmsd_time_temp_'+rodada[rod]] = mat_contents['rmsd_time_temp_'+rodada[rod]]
        locals()['rmsd_time_saln_'+rodada[rod]] = mat_contents['rmsd_time_saln_'+rodada[rod]]
        q_buoys = mat_contents['q_buoys']
        levels = np.squeeze(mat_contents['levels'])
    
        cond = abs(mat_contents['rmsd_profile_temp_'+rodada[rod]]) > 99999
        locals()['rmsd_profile_temp_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_temp_'+rodada[rod]])
        cond = abs(mat_contents['rmsd_profile_saln_'+rodada[rod]]) > 99999
        locals()['rmsd_profile_saln_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_profile_saln_'+rodada[rod]])
        cond = abs(mat_contents['rmsd_time_temp_'+rodada[rod]]) > 99999
        locals()['rmsd_time_temp_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_temp_'+rodada[rod]])
        cond = abs(mat_contents['rmsd_time_saln_'+rodada[rod]]) > 99999
        locals()['rmsd_time_saln_'+rodada[rod]] = np.ma.masked_where(cond,mat_contents['rmsd_time_saln_'+rodada[rod]])
    
        if rod==0:
           all_rod = rodada[rod]
        else:
           all_rod = all_rod+"_"+rodada[rod]
    
    plt.figure(figsize=(16, 8), dpi=150)
    plt.rcParams.update({'font.size': 20})
    for rod in range(0,len(rodada),1):
        out_dir = fig_output_dir+"/"+all_rod+"/RMSD_ARGO"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)
    
        date_vec = pd.date_range(data_inicial,data_final,freq='d')
        cond = np.array(levels)<=max_prof_rmsd_time
        temp = locals()['rmsd_time_temp_'+rodada[rod]]
        temp = np.mean(temp[:,cond],axis=1)
        data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq='d'))
        data_plot.plot(label=legenda[rod], color=line_color[rod])
        if rod==len(rodada)-1:
           plt.ylabel('\N{DEGREE SIGN}C')
           plt.legend()
           plt.title('TEMPERATURE RMSD ('+dominio_str+')')
           plt.savefig(out_dir+"/ARGO_HYCOM_RMSD_TEMP_TIME_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                      "_"+dominio_str+".png",dpi=resolution_plot)
           plt.close()
        del date_vec, data_plot

    plt.figure(figsize=(16, 8), dpi=150)
    plt.rcParams.update({'font.size': 20})
    for rod in range(0,len(rodada),1):
        out_dir = fig_output_dir+"/"+all_rod+"/RMSD_ARGO"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)

        date_vec = pd.date_range(data_inicial,data_final,freq='d')
        cond = np.array(levels)<=max_prof_rmsd_time
        temp = locals()['rmsd_time_saln_'+rodada[rod]]
        temp = np.mean(temp[:,cond],axis=1)
        data_plot = pd.Series(data=temp, index = pd.date_range(data_inicial,data_final,freq='d'))
        data_plot.plot(label=legenda[rod], color=line_color[rod])
        if rod==len(rodada)-1:
           plt.ylabel('PSU')
           plt.legend()
           plt.title('SALINITY RMSD ('+dominio_str+')')
           plt.savefig(out_dir+"/ARGO_HYCOM_RMSD_SALT_TIME_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                      "_"+dominio_str+".png",dpi=resolution_plot)
           plt.close()
        del date_vec, data_plot        
    
    plt.figure(figsize=(16, 8), dpi=150)
    plt.rcParams.update({'font.size': 20})
    for rod in range(0,len(rodada),1):
        cond = np.array(levels)<=max_prof_rmsd_profile
        temp = np.squeeze(locals()['rmsd_profile_temp_'+rodada[rod]])
        temp = temp[cond]
        plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
        if rod==len(rodada)-1:
           plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
           plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
           ax = plt.gca()
           ax.invert_yaxis()
           plt.xlabel('\N{DEGREE SIGN}C')
           plt.ylabel('DEPTH')
           plt.legend()
           plt.title('TEMPERATURE RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
           plt.savefig(out_dir+"/ARGO_HYCOM_RMSD_TEMP_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                      "_"+dominio_str+".png",dpi=resolution_plot)
           plt.close()
    
    plt.figure(figsize=(16, 8), dpi=150)
    plt.rcParams.update({'font.size': 20})
    for rod in range(0,len(rodada),1):
        cond = np.array(levels)<=max_prof_rmsd_profile
        temp = np.squeeze(locals()['rmsd_profile_saln_'+rodada[rod]])
        temp = temp[cond]
        plt.plot(temp, levels[cond], label=legenda[rod], color=line_color[rod])
        if rod==len(rodada)-1:
           plt.yticks(np.arange(0,max_prof_rmsd_profile,200))
           plt.grid(which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
           ax = plt.gca()
           ax.invert_yaxis()
           plt.xlabel('PSU')
           plt.ylabel('DEPTH')
           plt.legend()
           plt.title('SALINITY RMSD ('+data_inicial.strftime("%d/%m/%Y")+"-"+data_final.strftime("%d/%m/%Y")+'; '+dominio_str+')')
           plt.savefig(out_dir+"/ARGO_HYCOM_RMSD_SALT_PROFILE_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+ \
                      "_"+dominio_str+".png",dpi=resolution_plot)
           plt.close()

