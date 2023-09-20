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

#output_dir = '/prj/rodas/filipe.costa/resultados'
#dir_hycom = '/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc'
output_dir = '/home/filipe/resultados'
dir_hycom = '/disco1/remo/data/hycom_ufba'
dir_obs = 'assim/slasst_gridded_map_error_or_argoZ_stats/Check'

woa_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, \
              85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, \
              400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, \
              1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, \
              1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, \
              2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, \
              3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, \
              4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]
woa_levels = np.array(woa_levels)

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       dir_ATL = dir_hycom+'/ATLj0.04/'
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    counter = -1

    while(current_data <= data_final):
        print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))
        counter = counter+1

        dir_obs_full = dir_ATL+expt[rod]+'/'+dir_obs+'/'+current_data.strftime("%Y%m%d")+'00'
        arq_obs = dir_obs_full+'/ROMS_ARGO_CHECK.ascii'
        if not (os.path.isfile(arq_obs)):
           print('ARGO FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_obs)
           rmsd_time_temp_back[counter,:].mask = True
           rmsd_time_salt_back[counter,:].mask = True
           rmsd_time_temp_ana[counter,:].mask = True
           rmsd_time_salt_ana[counter,:].mask = True
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue
           exit()

        for i in range(0,int(inc_tempo),1):
            arq_dia_back_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1-i)).strftime("%Y%m%d")+'00/archv.' \
                                 +str((current_data + datetime.timedelta(days=-i)).timetuple().tm_year).zfill(4)+'_' \
                                 +str((current_data + datetime.timedelta(days=-i)).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            #print(arq_dia_back_hycom)
            if not (os.path.isfile(arq_dia_back_hycom)):
               print('HYCOM BACK FILE FOR DAY: '+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_back_hycom)
               exit()

        arq_dia_ana_hycom = dir_expt+'/output/Ncdf/'+current_data.strftime("%Y%m%d")+'00/archv.' \
                             +str(current_data.timetuple().tm_year).zfill(4)+'_' \
                             +str(current_data.timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
        #print(arq_dia_ana_hycom)
        if not (os.path.isfile(arq_dia_ana_hycom)):
           print('HYCOM ANA FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_dia_ana_hycom)
           exit()

        opt=open(arq_obs, "r")
        content_list = opt.readlines()
        obs = np.ma.array(np.zeros((len(content_list),10)), mask=True)
        for i in range(0,len(content_list),1):
            cont = content_list[i].rstrip('\n').split(' ')
            cont = list(filter(None,cont))
            cont.pop(8)  #ELIMINANDO OS STRINGS
            cont.pop(8)  #ELIMINANDO OS STRINGS
            cont = np.array(cont)
            obs[i] = cont[:]

        if (current_data==data_inicial):
           nc = Dataset(arq_dia_back_hycom, 'r')
           lon_hycom = nc.variables['Longitude'][:]
           lat_hycom = nc.variables['Latitude'][:]
           depth_hycom = nc.variables['Depth'][:]
           nc.close()
           del nc

           #max_lon_hycom = max(lon_hycom)-2
           #min_lon_hycom = min(lon_hycom)
           #max_lat_hycom = max(lat_hycom)-2
           #min_lat_hycom = min(lat_hycom)+2
           max_lon_hycom = max(lon_hycom)
           min_lon_hycom = min(lon_hycom)
           max_lat_hycom = max(lat_hycom)
           min_lat_hycom = min(lat_hycom)
           niv_max_hycom = max(depth_hycom)

           mean_profile_obs_temp = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_obs_temp = mean_profile_obs_temp.reshape(len(woa_levels))
           mean_profile_obs_salt = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_obs_salt = mean_profile_obs_salt.reshape(len(woa_levels))

           mean_profile_temp_back = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_temp_back = mean_profile_temp_back.reshape(len(woa_levels))
           mean_profile_salt_back = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_salt_back = mean_profile_salt_back.reshape(len(woa_levels))
           mean_profile_temp_ana = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_temp_ana = mean_profile_temp_ana.reshape(len(woa_levels))
           mean_profile_salt_ana = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           mean_profile_salt_ana = mean_profile_salt_ana.reshape(len(woa_levels))

           rmsd_profile_temp_back = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           rmsd_profile_temp_back = rmsd_profile_temp_back.reshape(len(woa_levels))
           rmsd_profile_salt_back = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           rmsd_profile_salt_back = rmsd_profile_salt_back.reshape(len(woa_levels))
           rmsd_profile_temp_ana = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           rmsd_profile_temp_ana = rmsd_profile_temp_ana.reshape(len(woa_levels))
           rmsd_profile_salt_ana = np.ma.array(np.zeros(len(woa_levels)), mask=False)
           rmsd_profile_salt_ana = rmsd_profile_salt_ana.reshape(len(woa_levels))

           mean_time_obs_temp = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_obs_temp = mean_time_obs_temp.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           mean_time_obs_salt = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_obs_salt = mean_time_obs_salt.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))

           mean_time_temp_back = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_temp_back = mean_time_temp_back.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           mean_time_salt_back = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_salt_back = mean_time_salt_back.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           mean_time_temp_ana = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_temp_ana = mean_time_temp_ana.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           mean_time_salt_ana = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           mean_time_salt_ana = mean_time_salt_ana.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))

           rmsd_time_temp_back = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           rmsd_time_temp_back = rmsd_time_temp_back.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           rmsd_time_salt_back = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           rmsd_time_salt_back = rmsd_time_salt_back.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           rmsd_time_temp_ana = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           rmsd_time_temp_ana = rmsd_time_temp_ana.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))
           rmsd_time_salt_ana = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           rmsd_time_salt_ana = rmsd_time_salt_ana.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))

           q_buoys_tot = np.ma.array(np.zeros(int(((data_final - data_inicial).days)/inc_tempo+1)*len(woa_levels)), mask=False)
           q_buoys_tot = q_buoys_tot.reshape(int(((data_final - data_inicial).days)/inc_tempo+1),len(woa_levels))

        cond = (obs[:,3] < max_lon_hycom) & (obs[:,3] > min_lon_hycom) & \
               (obs[:,4] < max_lat_hycom) & (obs[:,4] > min_lat_hycom)
        obs = obs[cond,:]
        del cond

        if not obs.any():
           rmsd_time_temp_back[counter,:].mask = True
           rmsd_time_salt_back[counter,:].mask = True
           rmsd_time_temp_ana[counter,:].mask = True
           rmsd_time_salt_ana[counter,:].mask = True
           print("NO OBS DIA: "+current_data.strftime("%d-%m-%Y"))
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        nc = Dataset(arq_dia_ana_hycom, 'r')
        t_ana_hycom = nc.variables['temperature'][0,:,:,:] #(33, 564, 502)
        t_ana_hycom = np.moveaxis(t_ana_hycom,0,-1)
        s_ana_hycom = nc.variables['salinity'][0,:,:,:]
        s_ana_hycom = np.moveaxis(s_ana_hycom,0,-1)

        nc.close()
        del nc

        dif_dia = np.diff(obs[:,5])
        cond = (dif_dia!=0)
        ind_dia = [i for i, val in enumerate(cond) if val]
        ind_dia.insert(0,-1)

        q_buoys_dia = np.ma.array(np.zeros((len(woa_levels))), mask=False)
        for j in range(0,len(ind_dia),1):
            data_obs = datetime.datetime(*time.strptime(str(int(obs[ind_dia[j]+1,5])).zfill(2) \
                       +str(int(obs[ind_dia[j]+1,6])).zfill(2)+str(int(obs[ind_dia[j]+1,7])),"%d%m%Y")[0:3])
            arq_back_hycom = dir_expt+'/output/Ncdf/'+(data_obs + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                             +str(data_obs.timetuple().tm_year).zfill(4)+'_' \
                             +str(data_obs.timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            nc = Dataset(arq_back_hycom, 'r')
            t_back_hycom = nc.variables['temperature'][0,:,:,:] #(33, 564, 502)
            t_back_hycom = np.moveaxis(t_back_hycom,0,-1)
            s_back_hycom = nc.variables['salinity'][0,:,:,:]
            s_back_hycom = np.moveaxis(s_back_hycom,0,-1)

            nc.close()
            del nc

            cond = (obs[:,5] == obs[ind_dia[j]+1,5]) & (obs[:,6] == obs[ind_dia[j]+1,6]) & (obs[:,7] == obs[ind_dia[j]+1,7])
            obs_dia = obs[cond,:]
            dif_lon = np.diff(obs_dia[:,3])
            dif_lat = np.diff(obs_dia[:,4])
            cond = (dif_lon!=0) | (dif_lat!=0)
            ind = [i for i, val in enumerate(cond) if val]
            ind.insert(0,-1)
            for i in range(0,len(ind),1):
                temp_back_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                salt_back_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                temp_ana_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                salt_ana_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                if i==len(ind)-1:
                   niv_max_obs = max(obs_dia[ind[i]+1:,0])
                   niv_min = min(obs_dia[ind[i]+1:,0])
                   cond = (woa_levels<=niv_max_obs) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom)
                   I = interpolate.interp1d(obs_dia[ind[i]+1:,0],obs_dia[ind[i]+1:,1])
                   obs_interp_temp = I(woa_levels[cond])
                   I = interpolate.interp1d(obs_dia[ind[i]+1:,0],obs_dia[ind[i]+1:,2])
                   obs_interp_salt = I(woa_levels[cond])
                else:
                   niv_max_obs = max(obs_dia[ind[i]+1:ind[i+1]+1,0])
                   niv_min = min(obs_dia[ind[i]+1:ind[i+1]+1,0])
                   cond = (woa_levels<=niv_max_obs) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom)
                   I = interpolate.interp1d(obs_dia[ind[i]+1:ind[i+1]+1,0],obs_dia[ind[i]+1:ind[i+1]+1,1])
                   obs_interp_temp = I(woa_levels[cond])
                   I = interpolate.interp1d(obs_dia[ind[i]+1:ind[i+1]+1,0],obs_dia[ind[i]+1:ind[i+1]+1,2])
                   obs_interp_salt = I(woa_levels[cond])

                ind_levs = [k for k, val in enumerate(cond) if val]
                ind_levs = np.array(ind_levs)
                for k in range(0,len(depth_hycom),1):
                    I = interpolate.interp2d(lon_hycom,lat_hycom,t_back_hycom[:,:,k],kind='linear')
                    temp_back_interp[k] = I(obs_dia[ind[i]+1,3] ,obs_dia[ind[i]+1,4])
                    I = interpolate.interp2d(lon_hycom,lat_hycom,s_back_hycom[:,:,k],kind='linear')
                    salt_back_interp[k] = I(obs_dia[ind[i]+1,3] ,obs_dia[ind[i]+1,4])
                    I = interpolate.interp2d(lon_hycom,lat_hycom,t_ana_hycom[:,:,k],kind='linear')
                    temp_ana_interp[k] = I(obs_dia[ind[i]+1,3] ,obs_dia[ind[i]+1,4])
                    I = interpolate.interp2d(lon_hycom,lat_hycom,s_ana_hycom[:,:,k],kind='linear')
                    salt_ana_interp[k] = I(obs_dia[ind[i]+1,3] ,obs_dia[ind[i]+1,4])
                I = interpolate.interp1d(depth_hycom,temp_back_interp)
                temp_back_interp = I(woa_levels[cond])
                I = interpolate.interp1d(depth_hycom,salt_back_interp)
                salt_back_interp = I(woa_levels[cond])
                I = interpolate.interp1d(depth_hycom,temp_ana_interp)
                temp_ana_interp = I(woa_levels[cond])
                I = interpolate.interp1d(depth_hycom,salt_ana_interp)
                salt_ana_interp = I(woa_levels[cond])
                cond = (abs(temp_back_interp) < 999999999) & (abs(salt_back_interp) < 999999999)

                q_buoys_tot[counter,ind_levs[cond]] = q_buoys_tot[counter,ind_levs[cond]]+1
                q_buoys_dia[ind_levs[cond]] = q_buoys_dia[ind_levs[cond]]+1

                mean_profile_obs_temp[ind_levs[cond]] = mean_profile_obs_temp[ind_levs[cond]] + obs_interp_temp[cond]
                mean_profile_obs_salt[ind_levs[cond]] = mean_profile_obs_salt[ind_levs[cond]] + obs_interp_salt[cond]
                mean_time_obs_temp[counter,ind_levs[cond]] = mean_time_obs_temp[counter,ind_levs[cond]] + obs_interp_temp[cond]
                mean_time_obs_salt[counter,ind_levs[cond]] = mean_time_obs_salt[counter,ind_levs[cond]] + obs_interp_salt[cond]

                mean_profile_temp_back[ind_levs[cond]] = mean_profile_temp_back[ind_levs[cond]] + temp_back_interp[cond]
                mean_profile_salt_back[ind_levs[cond]] = mean_profile_salt_back[ind_levs[cond]] + salt_back_interp[cond]
                mean_profile_temp_ana[ind_levs[cond]] = mean_profile_temp_ana[ind_levs[cond]] + temp_ana_interp[cond]
                mean_profile_salt_ana[ind_levs[cond]] = mean_profile_salt_ana[ind_levs[cond]] + salt_ana_interp[cond]

                mean_time_temp_back[counter,ind_levs[cond]] = mean_time_temp_back[counter,ind_levs[cond]] + temp_back_interp[cond]
                mean_time_salt_back[counter,ind_levs[cond]] = mean_time_salt_back[counter,ind_levs[cond]] + salt_back_interp[cond]
                mean_time_temp_ana[counter,ind_levs[cond]] = mean_time_temp_ana[counter,ind_levs[cond]] + temp_ana_interp[cond]
                mean_time_salt_ana[counter,ind_levs[cond]] = mean_time_salt_ana[counter,ind_levs[cond]] + salt_ana_interp[cond]

                rmsd_profile_temp_back[ind_levs[cond]] = rmsd_profile_temp_back[ind_levs[cond]] + \
                                                         (temp_back_interp[cond]-obs_interp_temp[cond])**2
                rmsd_profile_salt_back[ind_levs[cond]] = rmsd_profile_salt_back[ind_levs[cond]] + \
                                                         (salt_back_interp[cond]-obs_interp_salt[cond])**2
                rmsd_profile_temp_ana[ind_levs[cond]] = rmsd_profile_temp_ana[ind_levs[cond]] + \
                                                         (temp_ana_interp[cond]-obs_interp_temp[cond])**2
                rmsd_profile_salt_ana[ind_levs[cond]] = rmsd_profile_salt_ana[ind_levs[cond]] + \
                                                         (salt_ana_interp[cond]-obs_interp_salt[cond])**2

                rmsd_time_temp_back[counter,ind_levs[cond]] = rmsd_time_temp_back[counter,ind_levs[cond]] + \
                                                              (temp_back_interp[cond]-obs_interp_temp[cond])**2
                rmsd_time_salt_back[counter,ind_levs[cond]] = rmsd_time_salt_back[counter,ind_levs[cond]] + \
                                                              (salt_back_interp[cond]-obs_interp_salt[cond])**2
                rmsd_time_temp_ana[counter,ind_levs[cond]] = rmsd_time_temp_ana[counter,ind_levs[cond]] + \
                                                              (temp_ana_interp[cond]-obs_interp_temp[cond])**2
                rmsd_time_salt_ana[counter,ind_levs[cond]] = rmsd_time_salt_ana[counter,ind_levs[cond]] + \
                                                              (salt_ana_interp[cond]-obs_interp_salt[cond])**2

        cond = q_buoys_dia!=0
        mean_time_obs_temp[counter,cond] = np.sqrt(np.divide(mean_time_obs_temp[counter,cond],q_buoys_dia[cond]))
        mean_time_obs_salt[counter,cond] = np.sqrt(np.divide(mean_time_obs_salt[counter,cond],q_buoys_dia[cond]))
        mean_time_temp_back[counter,cond] = np.sqrt(np.divide(mean_time_temp_back[counter,cond],q_buoys_dia[cond]))
        mean_time_salt_back[counter,cond] = np.sqrt(np.divide(mean_time_salt_back[counter,cond],q_buoys_dia[cond]))
        mean_time_temp_ana[counter,cond] = np.sqrt(np.divide(mean_time_temp_ana[counter,cond],q_buoys_dia[cond]))
        mean_time_salt_ana[counter,cond] = np.sqrt(np.divide(mean_time_salt_ana[counter,cond],q_buoys_dia[cond]))

        rmsd_time_temp_back[counter,cond] = np.sqrt(np.divide(rmsd_time_temp_back[counter,cond],q_buoys_dia[cond]))
        rmsd_time_salt_back[counter,cond] = np.sqrt(np.divide(rmsd_time_salt_back[counter,cond],q_buoys_dia[cond]))
        rmsd_time_temp_ana[counter,cond] = np.sqrt(np.divide(rmsd_time_temp_ana[counter,cond],q_buoys_dia[cond]))
        rmsd_time_salt_ana[counter,cond] = np.sqrt(np.divide(rmsd_time_salt_ana[counter,cond],q_buoys_dia[cond]))

        cond = q_buoys_dia==0
        mean_time_obs_temp[counter,:] = np.ma.masked_where(cond,mean_time_obs_temp[counter,:])
        mean_time_obs_salt[counter,:] = np.ma.masked_where(cond,mean_time_obs_salt[counter,:])
        mean_time_temp_back[counter,:] = np.ma.masked_where(cond,mean_time_temp_back[counter,:])
        mean_time_salt_back[counter,:] = np.ma.masked_where(cond,mean_time_salt_back[counter,:])
        mean_time_temp_ana[counter,:] = np.ma.masked_where(cond,mean_time_temp_ana[counter,:])
        mean_time_salt_ana[counter,:] = np.ma.masked_where(cond,mean_time_salt_ana[counter,:])

        rmsd_time_temp_back[counter,:] = np.ma.masked_where(cond,rmsd_time_temp_back[counter,:])
        rmsd_time_salt_back[counter,:] = np.ma.masked_where(cond,rmsd_time_salt_back[counter,:])
        rmsd_time_temp_ana[counter,:] = np.ma.masked_where(cond,rmsd_time_temp_ana[counter,:])
        rmsd_time_salt_ana[counter,:] = np.ma.masked_where(cond,rmsd_time_salt_ana[counter,:])

        print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" --OK")
        current_data = current_data + datetime.timedelta(days=inc_tempo)

    tot = np.sum(q_buoys_tot,axis=0)
    cond = tot==0
    mean_profile_obs_temp = np.ma.masked_where(cond,mean_profile_obs_temp)
    mean_profile_obs_salt = np.ma.masked_where(cond,mean_profile_obs_salt)

    rmsd_profile_temp_back = np.ma.masked_where(cond,rmsd_profile_temp_back)
    rmsd_profile_salt_back = np.ma.masked_where(cond,rmsd_profile_salt_back)
    rmsd_profile_temp_ana = np.ma.masked_where(cond,rmsd_profile_temp_ana)
    rmsd_profile_salt_ana = np.ma.masked_where(cond,rmsd_profile_salt_ana)

    cond = tot!=0
    mean_profile_obs_temp[cond] = np.sqrt(np.divide(mean_profile_obs_temp[cond],tot[cond]))
    mean_profile_obs_salt[cond] = np.sqrt(np.divide(mean_profile_obs_salt[cond],tot[cond]))

    rmsd_profile_temp_back[cond] = np.sqrt(np.divide(rmsd_profile_temp_back[cond],tot[cond]))
    rmsd_profile_salt_back[cond] = np.sqrt(np.divide(rmsd_profile_salt_back[cond],tot[cond]))
    rmsd_profile_temp_ana[cond] = np.sqrt(np.divide(rmsd_profile_temp_ana[cond],tot[cond]))
    rmsd_profile_salt_ana[cond] = np.sqrt(np.divide(rmsd_profile_salt_ana[cond],tot[cond]))

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

    filename = "RMSD_ROMS_HYCOM_ANABACK_ARGOPOINTS_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+rodada[rod])):
       os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'rmsd_profile_temp_back_'+rodada[rod]: rmsd_profile_temp_back, \
                    'rmsd_profile_temp_ana_'+rodada[rod]: rmsd_profile_temp_ana, 'rmsd_profile_salt_back_'+rodada[rod]: rmsd_profile_salt_back, \
                    'rmsd_profile_salt_ana_'+rodada[rod]: rmsd_profile_salt_ana, 'rmsd_time_temp_back_'+rodada[rod]: rmsd_time_temp_back, \
                    'rmsd_time_temp_ana_'+rodada[rod]: rmsd_time_temp_ana, 'rmsd_time_salt_back_'+rodada[rod]: rmsd_time_salt_back, \
                    'rmsd_time_salt_ana_'+rodada[rod]: rmsd_time_salt_ana, \
                    'mean_profile_obs_temp': mean_profile_obs_temp, 'mean_profile_obs_salt': mean_profile_obs_salt, \
                    'mean_time_obs_temp': mean_time_obs_temp, 'mean_time_obs_salt': mean_time_obs_salt, \
                    'mean_profile_temp_back_'+rodada[rod]: mean_profile_temp_back, 'mean_profile_temp_ana_'+rodada[rod]: mean_profile_temp_ana,\
                    'mean_profile_salt_back_'+rodada[rod]: mean_profile_salt_back, 'mean_profile_salt_ana_'+rodada[rod]: mean_profile_salt_ana, \
                    'mean_time_temp_back_'+rodada[rod]: mean_time_temp_back, 'mean_time_temp_ana_'+rodada[rod]: mean_time_temp_ana,\
                    'mean_time_salt_back_'+rodada[rod]: mean_time_salt_back, 'mean_time_salt_ana_'+rodada[rod]: mean_time_salt_ana, \
                    'levels': woa_levels, 'q_buoys_tot': q_buoys_tot})
    del rmsd_profile_salt_back, rmsd_profile_temp_back, rmsd_time_temp_back, rmsd_time_salt_back, mean_profile_obs_temp, mean_profile_obs_salt, \
        mean_time_obs_temp, mean_time_obs_salt, mean_profile_temp_back, mean_profile_salt_back, mean_time_temp_back, mean_time_salt_back, \
        rmsd_profile_salt_ana, rmsd_profile_temp_ana, rmsd_time_temp_ana, rmsd_time_salt_ana, mean_profile_temp_ana, mean_profile_salt_ana, \
        mean_time_temp_ana, mean_time_salt_ana

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

exit()

