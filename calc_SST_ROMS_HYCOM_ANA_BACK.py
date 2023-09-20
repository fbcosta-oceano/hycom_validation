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

output_dir = '/prj/rodas/filipe.costa/resultados'
dir_hycom = '/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc'
dir_obs_roms = '/prj/rodas/leonardo.pires/ROMS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       record = 11
       dir_ATL = dir_hycom+'/ATLj0.04/'
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1

    while(current_data <= data_final):
        print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))
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

           max_lon_hycom = max(lon_hycom)-2
           min_lon_hycom = min(lon_hycom)
           max_lat_hycom = max(lat_hycom)-2
           min_lat_hycom = min(lat_hycom)+2

        arq_back_hycom = dir_expt+'/output/ab/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'/archv.' \
                        +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                        +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        arq_ana_hycom = dir_expt+'/output/ab/'+current_data.strftime("%Y%m%d")+'/archv.' \
                        +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                        +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        if not (os.path.isfile(arq_back_hycom)) or not (os.path.isfile(arq_ana_hycom)):
           print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_back_hycom)
           print(arq_ana_hycom)
           exit()

        counter = counter+1

        if (current_data==data_inicial):
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

           sst1 = nc1.variables['temp'][:,31,:,:] #(238, 32, 564, 502)
           sst2 = nc2.variables['temp'][:,31,:,:]
           sst3 = nc3.variables['temp'][:,31,:,:] #(238, 32, 564, 502)
           sst4 = nc4.variables['temp'][:,31,:,:]
           sst = np.ma.concatenate((sst1,sst2,sst3,sst4),axis=0)
           sst = np.moveaxis(sst,0,-1)
           sst = sst[cond_lat_roms,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1:int(inc_tempo)]
           sst_roms = sst[:,cond_lon_roms,:]
           nc1.close()
           nc2.close()
           nc3.close()
           nc4.close()
           del sst1, sst2, sst3, sst4, sst, nc1, nc2, nc3, nc4

        f = open(arq_back_hycom,'rb')
        f.seek(record*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_roms,lat_roms)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        if 'sst_hycom_back' in locals():
           sst_hycom_back = np.ma.concatenate((sst_hycom_back,field[:,:,np.newaxis]),axis=2)
        else:
           sst_hycom_back = field[:,:,np.newaxis]
        f.close()

        f = open(arq_back_hycom,'rb')
        f.seek(record*4*(IJDM+npad))
        field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
        field = np.reshape(field,(JDM,IDM))
        I = interpolate.interp2d(lon_hycom,lat_hycom,field,kind='linear')
        field = I(lon_roms,lat_roms)
        cond = abs(field) > 999999
        field = np.ma.masked_where(cond,field)
        if 'sst_hycom_ana' in locals():
           sst_hycom_ana = np.ma.concatenate((sst_hycom_ana,field[:,:,np.newaxis]),axis=2)
        else:
           sst_hycom_ana = field[:,:,np.newaxis]
        f.close()

        current_data = current_data + datetime.timedelta(days=inc_tempo)

    sst_hycom_mean_back = np.ma.array(np.zeros(sst_hycom_back[:,:,0].shape), mask=True)
    sst_hycom_std_back = np.ma.array(np.zeros(sst_hycom_back[:,:,0].shape), mask=True)
    sst_hycom_mean_ana = np.ma.array(np.zeros(sst_hycom_ana[:,:,0].shape), mask=True)
    sst_hycom_std_ana = np.ma.array(np.zeros(sst_hycom_ana[:,:,0].shape), mask=True)
    sst_roms_mean = np.ma.array(np.zeros(sst_roms[:,:,0].shape), mask=True)
    sst_roms_std = np.ma.array(np.zeros(sst_roms[:,:,0].shape), mask=True)

    sst_rmsd_map_back = np.ma.array(np.zeros(sst_roms[:,:,0].shape), mask=True)
    sst_rmsd_time_back = np.ma.array(np.zeros(len(sst_roms[0,0,:])), mask=True)
    sst_rmsd_map_ana = np.ma.array(np.zeros(sst_roms[:,:,0].shape), mask=True)
    sst_rmsd_time_ana = np.ma.array(np.zeros(len(sst_roms[0,0,:])), mask=True)

    sst_hycom_mean_back = np.mean(sst_hycom_back,axis=2)
    sst_hycom_mean_time_back = np.mean(sst_hycom_back,axis=(0,1))
    sst_hycom_inov_time = np.mean((sst_roms-sst_hycom_back),axis=(0,1))
    sst_hycom_std_back = np.std(sst_hycom_back,axis=2)
    sst_hycom_mean_ana = np.mean(sst_hycom_ana,axis=2)
    sst_hycom_mean_time_ana = np.mean(sst_hycom_ana,axis=(0,1))
    sst_hycom_std_ana = np.std(sst_hycom_ana,axis=2)
    sst_roms_mean = np.mean(sst_roms,axis=2)
    sst_roms_time_mean = np.mean(sst_roms,axis=(0,1))
    sst_roms_std = np.std(sst_roms,axis=2)
    sst_rmsd_map_back = np.sqrt(np.mean((sst_hycom_back-sst_roms)**2,axis=2))
    sst_rmsd_time_back =  np.sqrt(np.mean((sst_hycom_back-sst_roms)**2,axis=(0,1)))
    sst_rmsd_map_ana = np.sqrt(np.mean((sst_hycom_ana-sst_roms)**2,axis=2))
    sst_rmsd_time_ana =  np.sqrt(np.mean((sst_hycom_ana-sst_roms)**2,axis=(0,1)))

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
    filename = "SST_ROMS_GRID_HYCOM_BACK_ANA_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+ \
               str(int(inc_tempo))+"d_step_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+rodada[rod])):
       os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_roms': lat_roms, 'lon_roms': lon_roms, \
                     'sst_hycom_mean_back_'+rodada[rod]: sst_hycom_mean_back, 'sst_hycom_mean_ana_'+rodada[rod]: sst_hycom_mean_ana, \
                     'sst_roms_mean': sst_roms_mean, 'sst_hycom_std_back_'+rodada[rod]: sst_hycom_std_back, \
                     'sst_hycom_std_ana_'+rodada[rod]: sst_hycom_std_ana, 'sst_roms_std': sst_roms_std, \
                     'sst_rmsd_time_back_'+rodada[rod]: sst_rmsd_time_back, 'sst_rmsd_time_ana_'+rodada[rod]: sst_rmsd_time_ana, \
                     'sst_inov_time_'+rodada[rod]: sst_hycom_inov_time, \
                     'sst_rmsd_map_back_'+rodada[rod]: sst_rmsd_map_back, 'sst_rmsd_map_ana_'+rodada[rod]: sst_rmsd_map_ana})

    del sst_roms, sst_hycom_mean_back, sst_hycom_mean_ana, sst_roms_mean, sst_hycom_std_back, sst_hycom_std_ana, sst_roms_std
    del sst_rmsd_time_back, sst_rmsd_time_ana, sst_rmsd_map_back, sst_rmsd_map_ana, sst_hycom_back, sst_hycom_ana, dominio_str
    del filename, current_data, counter

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

