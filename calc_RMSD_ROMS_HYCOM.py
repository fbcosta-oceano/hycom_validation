from netCDF4 import Dataset
from scipy import interpolate
import scipy.io
import remo
import locale
import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange
import sys
import os
import time 
import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
jet = matplotlib.cm.get_cmap('jet')

lat_step_plot = 4
lon_step_plot = 5
resolution_plot = 150

pd.options.display.float_format = '{:.2f}'.format
np.set_printoptions(threshold=sys.maxsize)

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))
assim_step = float(content_list[6].rstrip('\n')) #A CADA QUANTOS DIAS ASSIMILOU NO EXPT
data_first_assim = str(content_list[7].rstrip('\n')) #DATA PRIMEIRA ASSIM

output_dir = '/home/filipe/resultados'
dir_hycom = '/disco1/remo/data/hycom_ufba'
dir_obs_roms = '/disco1/remo/leo/DADOS_ROMS/ROMS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
data_first_assim = datetime.datetime(*time.strptime(data_first_assim,"%d/%m/%Y")[0:3])
while(data_first_assim < data_inicial):
    data_first_assim = data_first_assim + datetime.timedelta(days=assim_step)

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       KDM = 32
       dir_ATL = dir_hycom+'/ATLj0.04/'
       record_thkn = 10
       record_temp = 11
       record_saln = 12
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1
    data_last_assim = data_first_assim

    while(current_data <= data_final):
        for assim in range(0,int(assim_step),1):
            if (current_data + datetime.timedelta(days=assim)>data_final):
               break
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))

            arq_roms = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/ROMS_'+current_data.strftime("%Y%m%d")+'.nc'
            if not (os.path.isfile(arq_roms)):
               print('ROMS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_roms)
               exit()

            if (current_data + datetime.timedelta(days=assim)==data_last_assim):
               arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data - datetime.timedelta(days=assim_step)).strftime("%Y%m%d")+'00/archv.' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            else:
               arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data).strftime("%Y%m%d")+'00/archv.' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'

            if not (os.path.isfile(arq_dia_hycom)):
               arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1+assim)).strftime("%Y%m%d")+'00/archv.' \
                               +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                               +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            #arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
            #                +str((current_data).timetuple().tm_year).zfill(4)+'_' \
            #                +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            if not (os.path.isfile(arq_dia_hycom)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom)
               exit()

            print(arq_dia_hycom)
            counter = counter+1
            if (current_data==data_inicial) and (assim==0):
               nc = Dataset(arq_dia_hycom, 'r')
               lon_hycom = nc.variables['Longitude'][:]
               lat_hycom = nc.variables['Latitude'][:]
               depth_hycom = nc.variables['Depth'][:]
               nc.close()
               del nc

               max_lon_hycom = max(lon_hycom)-2
               min_lon_hycom = min(lon_hycom)
               max_lat_hycom = max(lat_hycom)-2
               min_lat_hycom = min(lat_hycom)+2
            
               nc = Dataset(arq_roms, 'r')
               lat_roms = nc.variables['lat_rho'][:] #564
               lon_roms = nc.variables['lon_rho'][:] #502
               cond_lat_roms = (lat_roms<=max_lat_hycom) & (lat_roms>=min_lat_hycom)
               cond_lon_roms = (lon_roms<=max_lon_hycom) & (lon_roms>=min_lon_hycom)
               lat_roms = lat_roms[cond_lat_roms]
               lon_roms = lon_roms[cond_lon_roms]
               z_roms = nc.variables['depth'][:] #102
               nc.close()
               del nc

               mean_profile_roms_temp = np.ma.array(np.zeros(len(z_roms)), mask=False)
               mean_profile_roms_temp = mean_profile_roms_temp.reshape(len(z_roms))
               mean_profile_roms_salt = np.ma.array(np.zeros(len(z_roms)), mask=False)
               mean_profile_roms_salt = mean_profile_roms_salt.reshape(len(z_roms))

               mean_profile_temp = np.ma.array(np.zeros(len(z_roms)), mask=False)
               mean_profile_temp = mean_profile_temp.reshape(len(z_roms))
               mean_profile_salt = np.ma.array(np.zeros(len(z_roms)), mask=False)
               mean_profile_salt = mean_profile_salt.reshape(len(z_roms))

               rmsd_profile_temp = np.ma.array(np.zeros(len(z_roms)), mask=False)
               rmsd_profile_temp = rmsd_profile_temp.reshape(len(z_roms))
               rmsd_profile_salt = np.ma.array(np.zeros(len(z_roms)), mask=False)
               rmsd_profile_salt = rmsd_profile_salt.reshape(len(z_roms))

               mean_time_roms_temp = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               mean_time_roms_temp = mean_time_roms_temp.reshape((data_final - data_inicial).days+1,len(z_roms))
               mean_time_roms_salt = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               mean_time_roms_salt = mean_time_roms_salt.reshape((data_final - data_inicial).days+1,len(z_roms))

               mean_time_temp = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               mean_time_temp = mean_time_temp.reshape((data_final - data_inicial).days+1,len(z_roms))
               mean_time_salt = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               mean_time_salt = mean_time_salt.reshape((data_final - data_inicial).days+1,len(z_roms))

               rmsd_time_temp = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               rmsd_time_temp = rmsd_time_temp.reshape((data_final - data_inicial).days+1,len(z_roms))
               rmsd_time_salt = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(z_roms)), mask=False)
               rmsd_time_salt = rmsd_time_salt.reshape((data_final - data_inicial).days+1,len(z_roms))

            nc = Dataset(arq_roms, 'r')
            t_roms = nc.variables['temp'][:] #(564, 502, 102)
            t_roms = t_roms[cond_lat_roms,:,:]
            t_roms = t_roms[:,cond_lon_roms,:]
            cond = abs(t_roms) > 999999999
            t_roms = np.ma.masked_where(cond,t_roms)
            s_roms = nc.variables['salt'][:] #(564, 502, 102)
            s_roms = s_roms[cond_lat_roms,:,:]
            s_roms = s_roms[:,cond_lon_roms,:]
            cond = abs(s_roms) > 999999999
            s_roms = np.ma.masked_where(cond,s_roms)

            nc = Dataset(arq_dia_hycom, 'r')
            t_hycom = nc.variables['temperature'][0,:,:,:] #(33, 564, 502)
            t_hycom = np.moveaxis(t_hycom,0,-1)
            s_hycom = nc.variables['salinity'][0,:,:,:]
            s_hycom = np.moveaxis(s_hycom,0,-1)

            nc.close()
            del nc, cond

            t_interp = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*t_hycom.shape[-1]), mask=True)
            t_interp = t_interp.reshape(len(lat_roms),len(lon_roms),t_hycom.shape[-1])

            s_interp = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*s_hycom.shape[-1]), mask=True)
            s_interp = s_interp.reshape(len(lat_roms),len(lon_roms),s_hycom.shape[-1])
            for lev in range(0,len(depth_hycom),1):
                I = interpolate.interp2d(lon_hycom,lat_hycom,t_hycom[:,:,lev],kind='linear')
                t_interp[:,:,lev] = I(lon_roms,lat_roms)
                I = interpolate.interp2d(lon_hycom,lat_hycom,s_hycom[:,:,lev],kind='linear')
                s_interp[:,:,lev] = I(lon_roms,lat_roms)

            I = interpolate.interp1d(depth_hycom,t_interp,axis=-1)
            t_interp = I(z_roms.data)
            I = interpolate.interp1d(depth_hycom,s_interp,axis=-1)
            s_interp = I(z_roms.data)
            cond = abs(t_interp) > 999999999
            t_interp = np.ma.masked_where(cond,t_interp)
            cond = abs(s_interp) > 999999999
            s_interp = np.ma.masked_where(cond,s_interp)

            del t_hycom, s_hycom, I, cond

            #mean_profile_roms_temp = mean_profile_roms_temp + np.mean(t_roms,axis=(0,1))
            #mean_profile_roms_salt = mean_profile_roms_salt + np.mean(s_roms,axis=(0,1))
            #mean_profile_temp = mean_profile_temp + np.mean(t_interp,axis=(0,1))
            #mean_profile_salt = mean_profile_salt + np.mean(s_interp,axis=(0,1))

            mean_time_roms_temp[counter,:] = np.mean(t_roms,axis=(0,1))
            mean_time_roms_salt[counter,:] = np.mean(s_roms,axis=(0,1))
            mean_time_temp[counter,:] = np.mean(t_interp,axis=(0,1))
            mean_time_salt[counter,:] = np.mean(s_interp,axis=(0,1))

            #rmsd_profile_temp = rmsd_profile_temp + np.mean((t_interp-t_roms)**2,axis=(0,1))
            #rmsd_profile_salt = rmsd_profile_salt + np.mean((s_interp-s_roms)**2,axis=(0,1))
            #rmsd_time_temp[counter,:] = np.sqrt(np.mean((t_interp-t_roms)**2,axis=(0,1)))
            #rmsd_time_salt[counter,:] = np.sqrt(np.mean((s_interp-s_roms)**2,axis=(0,1)))
            rmsd_time_temp[counter,:] = np.mean((t_interp-t_roms)**2,axis=(0,1))
            rmsd_time_salt[counter,:] = np.mean((s_interp-s_roms)**2,axis=(0,1))

        #print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" --OK")
        current_data = current_data + datetime.timedelta(days=assim_step)

        del t_roms, s_roms, t_interp, s_interp

    #rmsd_profile_temp = np.sqrt(np.divide(rmsd_profile_temp,(data_final - data_inicial).days+1))
    #rmsd_profile_salt = np.sqrt(np.divide(rmsd_profile_salt,(data_final - data_inicial).days+1))

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

    counter_inic = 0
    current_data = data_inicial
    while(current_data <= data_final):
        days = monthrange(current_data.year,current_data.month)[1]
        counter_fim = counter_inic+days
        counter_fim = min(counter_fim,len(rmsd_time_temp[:,0]))

        rmsd_profile_temp = np.mean(rmsd_time_temp[counter_inic:counter_fim,:],axis=0)
        rmsd_profile_salt = np.mean(rmsd_time_salt[counter_inic:counter_fim,:],axis=0)

        mean_profile_temp = np.mean(mean_time_temp[counter_inic:counter_fim,:],axis=0)
        mean_profile_salt = np.mean(mean_time_salt[counter_inic:counter_fim,:],axis=0)
        mean_profile_roms_temp = np.mean(mean_time_roms_temp[counter_inic:counter_fim,:],axis=0)
        mean_profile_roms_salt = np.mean(mean_time_roms_salt[counter_inic:counter_fim,:],axis=0)

        if not (os.path.isdir(output_dir+"/"+rodada[rod])):
           os.system("mkdir -p "+output_dir+"/"+rodada[rod])
        filename = "RMSD_ROMS_HYCOM_"+rodada[rod]+"_"+current_data.strftime("%Y%m%d")+"-"+ \
                   (current_data + datetime.timedelta(days=counter_fim-counter_inic-1)).strftime("%Y%m%d")+"_"+ \
                   dominio_str
        scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'rmsd_profile_temp_'+rodada[rod]: rmsd_profile_temp, \
                         'rmsd_profile_salt_'+rodada[rod]: rmsd_profile_salt, \
                         'rmsd_time_temp_'+rodada[rod]: rmsd_time_temp[counter_inic:counter_fim,:], \
                         'rmsd_time_salt_'+rodada[rod]: rmsd_time_salt[counter_inic:counter_fim,:], \
                         'mean_profile_roms_temp': mean_profile_roms_temp, \
                         'mean_profile_roms_salt': mean_profile_roms_salt, \
                         'mean_time_roms_temp': mean_time_roms_temp[counter_inic:counter_fim,:], \
                         'mean_time_roms_salt': mean_time_roms_salt[counter_inic:counter_fim,:], \
                         'mean_profile_temp_'+rodada[rod]: mean_profile_temp, \
                         'mean_profile_salt_'+rodada[rod]: mean_profile_salt, \
                         'mean_time_temp_'+rodada[rod]: mean_time_temp[counter_inic:counter_fim,:], \
                         'mean_time_salt_'+rodada[rod]: mean_time_salt[counter_inic:counter_fim,:], 'levels': z_roms})

        counter_inic = counter_fim
        current_data = current_data + datetime.timedelta(days=days)

    #filename = "RMSD_ROMS_HYCOM_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    #if not (os.path.isdir(output_dir+"/"+rodada[rod])):
    #   os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    #scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'rmsd_profile_temp_'+rodada[rod]: rmsd_profile_temp, \
    #                 'rmsd_profile_salt_'+rodada[rod]: rmsd_profile_salt, 'rmsd_time_temp_'+rodada[rod]: rmsd_time_temp, \
    #                 'rmsd_time_salt_'+rodada[rod]: rmsd_time_salt, 'mean_profile_roms_temp': mean_profile_roms_temp, \
    #                 'mean_profile_roms_salt': mean_profile_roms_salt, 'mean_time_roms_temp': mean_time_roms_temp, \
    #                 'mean_time_roms_salt': mean_time_roms_salt, 'mean_profile_temp_'+rodada[rod]: mean_profile_temp, \
    #                 'mean_profile_salt_'+rodada[rod]: mean_profile_salt, 'mean_time_temp_'+rodada[rod]: mean_time_temp, \
    #                 'mean_time_salt_'+rodada[rod]: mean_time_salt, 'levels': z_roms})
    del rmsd_profile_salt, rmsd_profile_temp, rmsd_time_temp, rmsd_time_salt, mean_profile_roms_temp, mean_profile_roms_salt, \
        mean_time_roms_temp, mean_time_roms_salt, mean_profile_temp, mean_profile_salt, mean_time_temp, mean_time_salt

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

exit()

