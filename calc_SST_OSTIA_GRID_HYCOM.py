from netCDF4 import Dataset
from scipy import interpolate
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
inc_tempo = float(content_list[5].rstrip('\n'))
inc_assim = float(content_list[6].rstrip('\n'))

cfg_file='dirs_calc.txt'
opt=open(cfg_file, "r")
content_list = opt.readlines()
output_dir = str(content_list[0].rstrip('\n'))
dir_hycom = str(content_list[1].rstrip('\n'))
dir_obs = str(content_list[3].rstrip('\n'))

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       record = 11
       dir_ATL = dir_hycom+'/ATLj0.04/'
    elif '008d' in rodada[rod]:
       IDM = 1717
       JDM = 2345
       record = 14
       dir_ATL = dir_hycom+'/ATLd0.08/'
    elif '008i' in rodada[rod]:
       IDM = 628
       JDM = 780
       record = 11
       dir_ATL = dir_hycom+'/ATLi0.08/'
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1

    while(current_data <= data_final):
        print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))
        if (current_data==data_inicial):
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

           #max_lon_hycom = max(lon_hycom)-2
           #min_lon_hycom = min(lon_hycom)
           #max_lat_hycom = max(lat_hycom)-2
           #min_lat_hycom = min(lat_hycom)+2

           max_lon_hycom = max(lon_hycom)
           min_lon_hycom = min(lon_hycom)
           max_lat_hycom = max(lat_hycom)
           min_lat_hycom = min(lat_hycom)

        arq_obs = dir_obs+'/'+current_data.strftime("%Y")+'/'+current_data.strftime("%Y%m%d")+ \
                  '120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB_REP-v02.0-fv02.0.nc'
        if not (os.path.isfile(arq_obs)):
           print('OBS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
                        +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                        +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        if not (os.path.isfile(arq_dia_hycom)):
           for i in range(2,int(inc_assim)+1,1):
               arq_dia_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'/archv.' \
                               +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                               +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
               if (os.path.isfile(arq_dia_hycom)):
                  break        
        if not (os.path.isfile(arq_dia_hycom)):
           print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_dia_hycom)
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        counter = counter+1

        nc = Dataset(arq_obs, 'r')
        if (current_data==data_inicial):
           lat_ostia = nc.variables['lat'][:]
           lon_ostia = nc.variables['lon'][:]
           cond_lat = (lat_ostia<=max_lat_hycom) & (lat_ostia>=min_lat_hycom)
           cond_lon = (lon_ostia<=max_lon_hycom) & (lon_ostia>=min_lon_hycom)
           lat_ostia = lat_ostia[cond_lat]
           lon_ostia = lon_ostia[cond_lon]
        sst = nc.variables['analysed_sst'][0,cond_lat,cond_lon]
        sst = sst-273.15
        cond = abs(sst) > 99999
        sst = np.ma.masked_where(cond,sst)
        if 'sst_ostia' in locals():
           sst_ostia = np.ma.concatenate((sst_ostia,sst[:,:,np.newaxis]),axis=2)
        else:
           sst_ostia = sst[:,:,np.newaxis]

        f = open(arq_dia_hycom,'rb')
        f.seek(record*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_ostia,lat_ostia)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        if 'sst_hycom' in locals():
           sst_hycom = np.ma.concatenate((sst_hycom,field[:,:,np.newaxis]),axis=2)
        else:
           sst_hycom = field[:,:,np.newaxis]
        f.close()

        current_data = current_data + datetime.timedelta(days=inc_tempo)

    sst_hycom_mean = np.ma.array(np.zeros(sst_hycom[:,:,0].shape), mask=True)
    sst_hycom_std = np.ma.array(np.zeros(sst_hycom[:,:,0].shape), mask=True)
    sst_ostia_mean = np.ma.array(np.zeros(sst_ostia[:,:,0].shape), mask=True)
    sst_ostia_std = np.ma.array(np.zeros(sst_ostia[:,:,0].shape), mask=True)

    sst_rmsd_map = np.ma.array(np.zeros(sst_ostia[:,:,0].shape), mask=True)
    sst_rmsd_time = np.ma.array(np.zeros(len(sst_ostia[0,0,:])), mask=True)

    for lat in range(0,sst_hycom[:,0,0].size,1):
        for lon in range(0,sst_hycom[0,:,0].size,1):
            if not np.ma.is_masked(sst_hycom[lat,lon]):
               sst_hycom_mean[lat,lon] = remo.nanmean(sst_hycom[lat,lon,:])
               sst_hycom_std[lat,lon] = np.std(sst_hycom[lat,lon,:])
            if not np.ma.is_masked(sst_ostia[lat,lon]):
               sst_ostia_mean[lat,lon] = remo.nanmean(sst_ostia[lat,lon,:])
               sst_ostia_std[lat,lon] = np.std(sst_ostia[lat,lon,:])
            if not np.ma.is_masked(sst_hycom[lat,lon]) and not np.ma.is_masked(sst_ostia[lat,lon]):
               sst_rmsd_map[lat,lon] = remo.rmsd(sst_ostia[lat,lon,:],sst_hycom[lat,lon,:])

    current_data = data_inicial
    counter = -1
    while(current_data <= data_final):
        counter = counter+1

        sst_hycom_temp = np.ma.masked_where(sst_ostia[:,:,-1].mask,sst_hycom[:,:,counter])
        sst_ostia_temp = np.ma.masked_where(sst_hycom[:,:,-1].mask,sst_ostia[:,:,counter])
        sst_rmsd_time[counter] = remo.rmsd(sst_hycom_temp[:,:].compressed(),sst_ostia_temp[:,:].compressed())

        del sst_hycom_temp, sst_ostia_temp
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
    filename = "SST_OSTIA_GRID_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+rodada[rod])):
       os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_ostia': lat_ostia, 'lon_ostia': lon_ostia, \
                     'sst_hycom_mean_'+rodada[rod]: sst_hycom_mean, 'sst_ostia_mean': sst_ostia_mean, \
                     'sst_hycom_std_'+rodada[rod]: sst_hycom_std, 'sst_ostia_std': sst_ostia_std, \
                     'sst_rmsd_time_'+rodada[rod]: sst_rmsd_time, 'sst_rmsd_map_'+rodada[rod]: sst_rmsd_map})

    del sst_ostia, sst_hycom_mean, sst_ostia_mean, sst_hycom_std, sst_ostia_std, sst_rmsd_time
    del sst_rmsd_map, sst_hycom, dominio_str, filename, current_data, counter

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

