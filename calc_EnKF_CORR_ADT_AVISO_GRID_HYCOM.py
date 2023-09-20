from netCDF4 import Dataset
from scipy import interpolate
from scipy.stats import pearsonr
import scipy.io
import remo
import datetime
import sys
import ast
import os
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
jet = matplotlib.cm.get_cmap('jet')

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = ast.literal_eval(content_list[0])
rodada = ast.literal_eval(content_list[1])
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))
expt = np.array(expt,dtype=object)
rodada = np.array(rodada,dtype=object)

output_dir = '/home/filipe.costa/resultados'
dir_hycom = '/home/filipe.costa/previsao/hycom_2_2_18/proc'
dir_obs = '/home/filipe.costa/dados_obs/get_AVISO/adt'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

if '004j' in rodada[0][0]:
   IDM = 502
   JDM = 564
   dir_ATL = dir_hycom+'/ATLj0.04/'
elif '008d' in rodada[0][0]:
   IDM = 1717
   JDM = 2345
   dir_ATL = dir_hycom+'/ATLd0.08/'
elif '008i' in rodada[0][0]:
   IDM = 628
   JDM = 780
   dir_ATL = dir_hycom+'/ATLi0.08/'
else:
   print('CASE '+rodada[0][0]+' NOT DEFINED')
   exit()

IJDM=IDM*JDM
npad=4096-(IJDM%4096)
counter = -1

offset = np.ma.array(np.zeros((data_final - data_inicial).days+1), mask=True)

for ii in range(0,len(rodada),1):
    current_data = data_inicial
    while(current_data <= data_final):
        print("DIA: "+current_data.strftime("%d-%m-%Y"))
        for rod in range(0,len(rodada[ii]),1):
            dir_expt = dir_ATL+expt[ii][rod]
            if (current_data==data_inicial) and (rod==0):
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

            arq_obs = dir_obs+'/'+current_data.strftime("%Y")+'/dt_global_allsat_phy_l4_'+current_data.strftime("%Y%m%d")+'.nc'
            if not (os.path.isfile(arq_obs)):
               print('OBS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_obs)
               current_data = current_data + datetime.timedelta(days=inc_tempo)
               exit()

            arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
                            +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
            if not (os.path.isfile(arq_dia_hycom)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom)
               current_data = current_data + datetime.timedelta(days=inc_tempo)
               exit()

            if (rod==0):
               counter = counter+1

               nc = Dataset(arq_obs, 'r')
               if (current_data==data_inicial):
                  lat_aviso = nc.variables['latitude'][:]
                  lon_aviso = nc.variables['longitude'][:]
                  cond_lat = (lat_aviso<=max_lat_hycom) & (lat_aviso>=min_lat_hycom)
                  cond_lon = (lon_aviso<=max_lon_hycom) & (lon_aviso>=min_lon_hycom)
                  lat_aviso = lat_aviso[cond_lat]
                  lon_aviso = lon_aviso[cond_lon]
               ssh = nc.variables['adt'][:]
               ssh = np.moveaxis(ssh,0,-1)
               cond = abs(ssh) > 99999
               ssh = np.ma.masked_where(cond,ssh)
               ssh = ssh[cond_lat,:,0]
               ssh = ssh[:,cond_lon]
               if 'ssh_aviso' in locals():
                  ssh_aviso = np.ma.concatenate((ssh_aviso,ssh[:,:,np.newaxis]),axis=2)
               else:
                  ssh_aviso = ssh[:,:,np.newaxis]

            f = open(arq_dia_hycom,'rb')
            f.seek(1*4*(IJDM+npad))
            field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
            field = np.reshape(field,(JDM,IDM))
            field = field/9.806
            I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
            field = I(lon_aviso,lat_aviso)
            cond = abs(field) > 999999
            field = np.ma.masked_where(cond,field)
            cond = abs(ssh_aviso[:,:,counter]) > 999999
            field = np.ma.masked_where(cond,field)
            if 'field_sum' in locals():
               field_sum = field_sum+field
            else:
               field_sum = field
            f.close()
            del field, cond
        field_sum = field_sum/len(rodada[ii])
        if counter==0:
           ssh_hycom = field_sum[:,:,np.newaxis]
        else:
           ssh_hycom = np.ma.concatenate((ssh_hycom,field_sum[:,:,np.newaxis]),axis=2)
        del field_sum

        offset[counter] = remo.nanmean(ssh_hycom[:,:,counter].compressed()) - remo.nanmean(ssh_aviso[:,:,counter].compressed())

        current_data = current_data + datetime.timedelta(days=inc_tempo)

    ssh_hycom_mean = np.ma.array(np.zeros(ssh_hycom[:,:,0].shape), mask=True)
    ssh_hycom_std = np.ma.array(np.zeros(ssh_hycom[:,:,0].shape), mask=True)
    ssh_aviso_mean = np.ma.array(np.zeros(ssh_aviso[:,:,0].shape), mask=True)
    ssh_aviso_std = np.ma.array(np.zeros(ssh_aviso[:,:,0].shape), mask=True)
    
    ssh_corr_map = np.ma.array(np.zeros(ssh_aviso[:,:,0].shape), mask=True)
    ssh_rmsd_map = np.ma.array(np.zeros(ssh_aviso[:,:,0].shape), mask=True)
    ssh_corr_time = np.ma.array(np.zeros(len(ssh_aviso[0,0,:])), mask=True)
    ssh_rmsd_time = np.ma.array(np.zeros(len(ssh_aviso[0,0,:])), mask=True)
    
    for lat in range(0,ssh_hycom[:,0,0].size,1):
        for lon in range(0,ssh_hycom[0,:,0].size,1):
            if not np.ma.is_masked(ssh_hycom[lat,lon]):
               ssh_hycom_mean[lat,lon] = remo.nanmean(ssh_hycom[lat,lon,:])
               ssh_hycom_std[lat,lon] = np.std(ssh_hycom[lat,lon,:])
            if not np.ma.is_masked(ssh_aviso[lat,lon]):
               ssh_aviso_mean[lat,lon] = remo.nanmean(ssh_aviso[lat,lon,:])
               ssh_aviso_std[lat,lon] = np.std(ssh_aviso[lat,lon,:])
            if not np.ma.is_masked(ssh_hycom[lat,lon]) and not np.ma.is_masked(ssh_aviso[lat,lon]):
               corr = pearsonr(ssh_hycom[lat,lon,:], ssh_aviso[lat,lon,:])
               ssh_corr_map[lat,lon] = corr[0]
               ssh_rmsd_map[lat,lon] = remo.rmsd(ssh_aviso[lat,lon,:],ssh_hycom[lat,lon,:]-offset[counter])
    
    current_data = data_inicial
    counter = -1
    while(current_data <= data_final):
        arq_obs = dir_obs+'/'+current_data.strftime("%Y")+'/dt_global_allsat_phy_l4_'+current_data.strftime("%Y%m%d")+'.nc'
        if not (os.path.isfile(arq_obs)):
           print('OBS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           exit()
    
        #arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
        #                +str((current_data).timetuple().tm_year).zfill(4)+'_' \
        #                +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        #if not (os.path.isfile(arq_dia_hycom)):
        #   print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
        #   current_data = current_data + datetime.timedelta(days=inc_tempo)
        #   continue
    
        counter = counter+1
    
        ssh_hycom_temp = np.ma.masked_where(ssh_aviso[:,:,-1].mask,ssh_hycom[:,:,counter])
        ssh_aviso_temp = np.ma.masked_where(ssh_hycom[:,:,-1].mask,ssh_aviso[:,:,counter])
        ssh_rmsd_time[counter] = np.mean((ssh_hycom_temp-offset[counter]-ssh_aviso_temp)**2,axis=(0,1))
        corr = pearsonr(ssh_hycom_temp[:,:].compressed(), ssh_aviso_temp[:,:].compressed())
        ssh_corr_time[counter] = corr[0]
    
        del ssh_hycom_temp, ssh_aviso_temp, corr
        current_data = current_data + datetime.timedelta(days=inc_tempo)
    
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
    dominio_str = f"{abs(min_lat_hycom):.2f}"+s1+"-"+f"{abs(max_lat_hycom):.2f}"+s2+"_"+ \
                  f"{abs(min_lon_hycom):.2f}"+s3+"-"+f"{abs(max_lon_hycom):.2f}"+s4
    
    all_rod = rodada[ii][0]
    for rod in range(1,len(rodada[ii]),1):
        all_rod = all_rod+"_"+rodada[ii][rod][5:]
    
    filename = "CORR_ADT_AVISO_GRID_HYCOM_"+all_rod+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+all_rod)):
       os.system("mkdir -p "+output_dir+"/"+rodada[ii][0])
    scipy.io.savemat(output_dir+"/"+rodada[ii][0]+"/"+filename+".mat", mdict={'lat_aviso': lat_aviso, 'lon_aviso': lon_aviso, \
                     'ssh_hycom_mean_'+rodada[ii][0]: ssh_hycom_mean, 'ssh_aviso_mean': ssh_aviso_mean, \
                     'ssh_hycom_std_'+rodada[ii][0]: ssh_hycom_std, 'ssh_aviso_std': ssh_aviso_std, \
                     'ssh_corr_time_'+rodada[ii][0]: ssh_corr_time, 'ssh_corr_map_'+rodada[ii][0]: ssh_corr_map, \
                     'ssh_rmsd_time_'+rodada[ii][0]: ssh_rmsd_time, 'ssh_rmsd_map_'+rodada[ii][0]: ssh_rmsd_map})
    
    del ssh_aviso, ssh_hycom_mean, ssh_aviso_mean, ssh_hycom_std, ssh_aviso_std, ssh_corr_time, ssh_corr_map, ssh_rmsd_time
    del ssh_rmsd_map, ssh_hycom, dominio_str, filename, current_data, counter

