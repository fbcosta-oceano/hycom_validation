from netCDF4 import Dataset
from scipy import interpolate
from scipy.stats import pearsonr
import scipy.io
import remo
import datetime
import sys
import os
import time
import numpy as np

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))  #INCREMENTO DA DATA PARA CALCULO ABAIXO
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
#if (data_first_assim>data_inicial):
#    data_first_assim = data_first_assim - datetime.timedelta(days=assim_step)

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       dir_ATL = dir_hycom+'/ATLj0.04/'
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       exit()

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1
    data_last_assim = data_first_assim

    while(current_data <= data_final):
        for assim in range(0,int(assim_step),1):
            if (current_data + datetime.timedelta(days=assim)>data_final):
               break
            print(rodada[rod]+" DIA: "+(current_data + datetime.timedelta(days=assim)).strftime("%d-%m-%Y"))
            if (current_data==data_inicial) and (assim==0):
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

               max_lon_hycom = max(lon_hycom)-2
               min_lon_hycom = min(lon_hycom)
               max_lat_hycom = max(lat_hycom)-2
               min_lat_hycom = min(lat_hycom)+2

            if (current_data + datetime.timedelta(days=assim)==data_last_assim):
               arq_dia_hycom = dir_expt+'/output/ab/'+(current_data - datetime.timedelta(days=assim_step)).strftime("%Y%m%d")+'/archv.' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00.a'
            else:
               arq_dia_hycom = dir_expt+'/output/ab/'+(current_data).strftime("%Y%m%d")+'/archv.' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00.a'
            if not (os.path.isfile(arq_dia_hycom)):
               arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1+assim)).strftime("%Y%m%d")+'/archv.' \
                               +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_year).zfill(4)+'_' \
                               +str((current_data + datetime.timedelta(days=assim)).timetuple().tm_yday).zfill(3)+'_00.a'
               if not (os.path.isfile(arq_dia_hycom)):
                  print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
                  print(arq_dia_hycom)
                  exit()

            print(arq_dia_hycom)
            counter = counter+1

            if (current_data==data_inicial) and (assim==0):
               #arq_roms1 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_1de2.nc4'
               #arq_roms2 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_2de2.nc4'
               arq_roms1 = dir_obs_roms+'/ROMS_2011/roms_avg_1de2.nc4'
               arq_roms2 = dir_obs_roms+'/ROMS_2011/roms_avg_2de2.nc4'
               arq_roms3 = dir_obs_roms+'/ROMS_2012/roms_avg_1de2.nc4'
               arq_roms4 = dir_obs_roms+'/ROMS_2012/roms_avg_2de2.nc4'
               if not (os.path.isfile(arq_roms1)) or not (os.path.isfile(arq_roms2)):
                  print('ROMS FILE FOR YEAR: '+current_data.strftime("%Y")+' DOES NOT EXIST')
                  print(arq_roms1)
                  print(arq_roms2)
                  print(arq_roms3)
                  print(arq_roms4)
                  exit()
               nc1 = Dataset(arq_roms1, 'r')
               nc2 = Dataset(arq_roms2, 'r')
               nc3 = Dataset(arq_roms3, 'r')
               nc4 = Dataset(arq_roms4, 'r')
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
               time_roms3 = nc3.variables['ocean_time'][:]
               time_roms4 = nc4.variables['ocean_time'][:]
               time_roms = np.concatenate((time_roms1,time_roms2,time_roms3,time_roms4),axis=0)
               time_roms = time_roms/(60*60*24) #AJUSTANDO PARA DIAS
               time_roms = time_roms-.5 #SUBTRAINDO 12 h PARA FICAR EM 00Z
               ref = datetime.datetime(*time.strptime('01/01/1970',"%d/%m/%Y")[0:3])
               date_roms = [ref + datetime.timedelta(days=time_roms[i]) for i in range(0,len(time_roms),1)]
               ind_inic_date_roms = np.where(date_roms == np.datetime64(datetime.datetime(data_inicial.year, data_inicial.month, data_inicial.day, 0, 0)))
               ind_fim_date_roms = np.where(date_roms == np.datetime64(datetime.datetime(data_final.year, data_final.month, data_final.day, 0, 0)))

               ssh1 = nc1.variables['zeta'][:] #(238, 564, 502)
               ssh2 = nc2.variables['zeta'][:]
               ssh3 = nc3.variables['zeta'][:]
               ssh4 = nc4.variables['zeta'][:]
               ssh = np.ma.concatenate((ssh1,ssh2,ssh3,ssh4),axis=0)
               ssh = np.moveaxis(ssh,0,-1)
               ssh = ssh[cond_lat_roms,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1]
               ssh_roms = ssh[:,cond_lon_roms,:]
               del ssh1, ssh2, ssh3, ssh4, ssh, cond_lat_roms, cond_lon_roms, time_roms1, time_roms2, time_roms3, time_roms4
               nc1.close()
               nc2.close()
               nc3.close()
               nc4.close()
               del nc1, nc2, nc3, nc4

               offset = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*((data_final - data_inicial).days+1)), mask=True)
               offset = np.reshape(offset,ssh_roms.shape)

            f = open(arq_dia_hycom,'rb')
            f.seek(1*4*(IJDM+npad))
            field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
            field = np.reshape(field,(JDM,IDM))
            field = field/9.806
            I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
            field = I(lon_roms,lat_roms)
            cond = abs(field) > 999999
            field = np.ma.masked_where(cond,field)
            if 'ssh_hycom' in locals():
               ssh_hycom = np.ma.concatenate((ssh_hycom,field[:,:,np.newaxis]),axis=2)
            else:
               ssh_hycom = field[:,:,np.newaxis]
            f.close()
            offset[:,:,counter] = remo.nanmean(ssh_hycom[:,:,counter].compressed()) - remo.nanmean(ssh_roms[:,:,counter].compressed())

        data_last_assim = data_last_assim + datetime.timedelta(days=assim_step)
        current_data = current_data + datetime.timedelta(days=assim_step)

    ssh_hycom_mean = np.ma.array(np.zeros(ssh_hycom[:,:,0].shape), mask=True)
    ssh_hycom_std = np.ma.array(np.zeros(ssh_hycom[:,:,0].shape), mask=True)
    ssh_roms_mean = np.ma.array(np.zeros(ssh_roms[:,:,0].shape), mask=True)
    ssh_roms_std = np.ma.array(np.zeros(ssh_roms[:,:,0].shape), mask=True)

    ssh_hycom_mean = np.mean(ssh_hycom,axis=2)
    ssh_hycom_std = np.std(ssh_hycom,axis=2)
    ssh_roms_mean = np.mean(ssh_roms,axis=2)
    ssh_roms_std = np.std(ssh_roms,axis=2)

    ssh_corr_map = np.ma.array(np.zeros(ssh_roms[:,:,0].shape), mask=True)
    ssh_rmsd_map = np.ma.array(np.zeros(ssh_roms[:,:,0].shape), mask=True)
    ssh_corr_time = np.ma.array(np.zeros(len(ssh_roms[0,0,:])), mask=True)
    ssh_rmsd_time = np.ma.array(np.zeros(len(ssh_roms[0,0,:])), mask=True)

    ssh_rmsd_map = np.sqrt(np.mean((ssh_hycom-offset-ssh_roms)**2,axis=2))
    ssh_rmsd_time = np.sqrt(np.mean((ssh_hycom-offset-ssh_roms)**2,axis=(0,1)))

    for lat in range(0,ssh_hycom[:,0,0].size,1):
        for lon in range(0,ssh_hycom[0,:,0].size,1):
            if not np.ma.is_masked(ssh_hycom[lat,lon]) and not np.ma.is_masked(ssh_roms[lat,lon]):
               corr = pearsonr(ssh_hycom[lat,lon,:], ssh_roms[lat,lon,:])
               ssh_corr_map[lat,lon] = corr[0]
               #ssh_rmsd_map[lat,lon] = remo.rmsd(ssh_roms[lat,lon,:],ssh_hycom[lat,lon,:]-offset)

    current_data = data_inicial
    counter = -1
    while(current_data <= data_final):
        #arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
        #                +str((current_data).timetuple().tm_year).zfill(4)+'_' \
        #                +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        #if not (os.path.isfile(arq_dia_hycom)):
        #   print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
        #   print(arq_dia_hycom)
        #   exit()

        counter = counter+1

        ssh_hycom_temp = np.ma.masked_where(ssh_roms[:,:,-1].mask,ssh_hycom[:,:,counter])
        ssh_roms_temp = np.ma.masked_where(ssh_hycom[:,:,-1].mask,ssh_roms[:,:,counter])
        #ssh_rmsd_time[counter] = remo.rmsd(ssh_hycom_temp[:,:].compressed()-offset[counter],ssh_roms_temp[:,:].compressed())
        corr = pearsonr(ssh_hycom_temp[:,:].compressed(), ssh_roms_temp[:,:].compressed())
        ssh_corr_time[counter] = corr[0]

        del ssh_hycom_temp, ssh_roms_temp, corr
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
    filename = "CORR_ADT_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+rodada[rod])):
       os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_roms': lat_roms, 'lon_roms': lon_roms, \
                     'ssh_hycom_mean_'+rodada[rod]: ssh_hycom_mean, 'ssh_roms_mean': ssh_roms_mean, \
                     'ssh_hycom_std_'+rodada[rod]: ssh_hycom_std, 'ssh_roms_std': ssh_roms_std, \
                     'ssh_corr_time_'+rodada[rod]: ssh_corr_time, 'ssh_corr_map_'+rodada[rod]: ssh_corr_map, \
                     'ssh_rmsd_time_'+rodada[rod]: ssh_rmsd_time, 'ssh_rmsd_map_'+rodada[rod]: ssh_rmsd_map})

    del ssh_roms, ssh_hycom_mean, ssh_roms_mean, ssh_hycom_std, ssh_roms_std, ssh_corr_time, ssh_corr_map, ssh_rmsd_time
    del ssh_rmsd_map, ssh_hycom, dominio_str, filename, current_data, counter

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

