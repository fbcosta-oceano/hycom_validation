from netCDF4 import Dataset
from scipy import interpolate
import scipy.io
import remo
import datetime
import sys
import os
import time
import numpy as np
from dateutil.relativedelta import relativedelta
from calendar import monthrange
import warnings
warnings.filterwarnings('ignore')

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
dir_obs_roms = '/prj/rodas/filipe.costa/ROMS/mma'

prof = np.array([[0,0], [0,50], [50,50], [100,100], [200,200], [300,300], [400,400], [800,800], [2000,2000]])
prof = np.array([[0,0], [0,50], [50,50], [100,100], [200,200], [300,300], [400,400], [800,800]])
prof = np.array([[0,0], [100,100], [400,400], [800,800]])
k_inic = 3

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    for k in range(k_inic,len(prof[:,0]),1):
        current_data = data_inicial
        if '004j' in rodada[rod]:
           dir_ATL = dir_hycom+'/ATLj0.04/'
        else:
           print('CASE '+rodada[rod]+' NOT DEFINED')
           continue

        dir_expt = dir_ATL+expt[rod]
        counter = -1

        while(current_data <= data_final):
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" PROF: "+str(prof[k,0])+"-"+str(prof[k,1]))
            arq_roms = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y%m%d")+'.nc'
            if not (os.path.isfile(arq_roms)):
               print('ROMS FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_roms)
               exit()

            arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                            +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'
            if not (os.path.isfile(arq_dia_hycom)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom)
               exit()

            counter = counter+1
            if (current_data==data_inicial):
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
            if (current_data==data_inicial):
               lat_roms = nc.variables['lat_rho'][:] #564
               lon_roms = nc.variables['lon_rho'][:] #502
               cond_lat_roms = (lat_roms<=max_lat_hycom) & (lat_roms>=min_lat_hycom)
               cond_lon_roms = (lon_roms<=max_lon_hycom) & (lon_roms>=min_lon_hycom)
               lat_roms = lat_roms[cond_lat_roms]
               lon_roms = lon_roms[cond_lon_roms]
               z_roms = nc.variables['depth'][:] #102
               ind1 = [x for x in enumerate(np.flip(z_roms)) if x[1] <= prof[k,0]][0]
               ind2 = [x for x in enumerate(z_roms) if x[1] >= prof[k,1]][0]
               ind1 = abs(ind1[0]+1-len(z_roms))
               indz_roms = np.arange(ind1,ind2[0]+1,1)
            u = nc.variables['u'][:] #(564, 502, 102)
            u = u[cond_lat_roms,:,:]
            u = u[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               u = np.squeeze(u[:,:,indz_roms])
            else:
               u = np.mean(u[:,:,indz_roms],axis=-1)
            if 'u_roms' in locals():
               u_roms = np.ma.concatenate((u_roms,u[:,:,np.newaxis]),axis=2)
            else:
               u_roms = u[:,:,np.newaxis]

            v = nc.variables['v'][:] #(564, 502, 102)
            v = v[cond_lat_roms,:,:]
            v = v[:,cond_lon_roms,:]
            if len(indz_roms)==1:
               v = np.squeeze(v[:,:,indz_roms])
            else:
               v = np.mean(v[:,:,indz_roms],axis=-1)
            if 'v_roms' in locals():
               v_roms = np.ma.concatenate((v_roms,v[:,:,np.newaxis]),axis=2)
            else:
               v_roms = v[:,:,np.newaxis]
            nc.close()
            del u, v, nc

            nc = Dataset(arq_dia_hycom, 'r')
            u = nc.variables['u'][0,:,:,:] #(33, 564, 502)
            u = np.moveaxis(u,0,-1)
            v = nc.variables['v'][0,:,:,:]
            v = np.moveaxis(v,0,-1)

            nc.close()
            del nc

            I = interpolate.interp1d(depth_hycom,u,axis=-1)
            u_interp = I(z_roms[indz_roms].data)
            I = interpolate.interp1d(depth_hycom,v,axis=-1)
            v_interp = I(z_roms[indz_roms].data)
            del u, v, I
            u = np.mean(u_interp,axis=-1)
            I = interpolate.interp2d(lon_hycom,lat_hycom,u,kind='linear')
            u = I(lon_roms,lat_roms)
            cond = abs(u) > 999999999
            u = np.ma.masked_where(cond,u)
            del I, cond

            v = np.mean(v_interp,axis=-1)
            I = interpolate.interp2d(lon_hycom,lat_hycom,v,kind='linear')
            v = I(lon_roms,lat_roms)
            cond = abs(v) > 999999999
            v = np.ma.masked_where(cond,v)
            del I, cond

            if 'u_hycom' in locals():
               u_hycom = np.ma.concatenate((u_hycom,u[:,:,np.newaxis]),axis=2)
               v_hycom = np.ma.concatenate((v_hycom,v[:,:,np.newaxis]),axis=2)
            else:
               u_hycom = u[:,:,np.newaxis]
               v_hycom = v[:,:,np.newaxis]

            del u, v, u_interp, v_interp
            current_data = current_data + datetime.timedelta(days=inc_tempo)

        vel_hycom = np.ma.array(np.zeros(v_hycom.shape), mask=True)
        vel_roms = np.ma.array(np.zeros(v_roms.shape), mask=True)
        vel_hycom = np.sqrt(u_hycom**2+v_hycom**2)
        vel_roms = np.sqrt(u_roms**2+v_roms**2)

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
            counter_fim = min(counter_fim,len(u_roms[0,0,:]))

            u_hycom_mean = np.ma.array(np.zeros(u_hycom[:,:,0].shape), mask=True)
            u_roms_mean = np.ma.array(np.zeros(u_roms[:,:,0].shape), mask=True)
            v_hycom_mean = np.ma.array(np.zeros(v_hycom[:,:,0].shape), mask=True)
            v_roms_mean = np.ma.array(np.zeros(v_roms[:,:,0].shape), mask=True)
            vel_hycom_mean = np.ma.array(np.zeros(v_hycom[:,:,0].shape), mask=True)
            vel_roms_mean = np.ma.array(np.zeros(v_roms[:,:,0].shape), mask=True)

            u_hycom_mean = np.mean(u_hycom[:,:,counter_inic:counter_fim],axis=2)
            u_roms_mean = np.mean(u_roms[:,:,counter_inic:counter_fim],axis=2)
            v_hycom_mean = np.mean(v_hycom[:,:,counter_inic:counter_fim],axis=2)
            v_roms_mean = np.mean(v_roms[:,:,counter_inic:counter_fim],axis=2)
            vel_hycom_mean = np.mean(vel_hycom[:,:,counter_inic:counter_fim],axis=2)
            vel_roms_mean = np.mean(vel_roms[:,:,counter_inic:counter_fim],axis=2)

            if not (os.path.isdir(output_dir+"/"+rodada[rod])):
               os.system("mkdir -p "+output_dir+"/"+rodada[rod])
            filename = "VEL_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+current_data.strftime("%Y%m%d")+"-"+ \
                       (current_data + datetime.timedelta(days=counter_fim-counter_inic-1)).strftime("%Y%m%d")+"_"+ \
                       dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"
            scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_roms': lat_roms, 'lon_roms': lon_roms, \
                             'u_hycom_mean_'+rodada[rod]: u_hycom_mean, 'v_hycom_mean_'+rodada[rod]: v_hycom_mean, \
                             'u_roms_mean': u_roms_mean, 'v_roms_mean': v_roms_mean, \
                             'vel_hycom_mean_'+rodada[rod]: vel_hycom_mean, 'vel_roms_mean': vel_roms_mean})

            del u_hycom_mean,u_roms_mean,v_hycom_mean,v_roms_mean,vel_hycom_mean,vel_roms_mean

            counter_inic = counter_fim
            current_data = current_data + datetime.timedelta(days=days)

        u_hycom_mean = np.ma.array(np.zeros(u_hycom[:,:,0].shape), mask=True)
        u_roms_mean = np.ma.array(np.zeros(u_roms[:,:,0].shape), mask=True)
        v_hycom_mean = np.ma.array(np.zeros(v_hycom[:,:,0].shape), mask=True)
        v_roms_mean = np.ma.array(np.zeros(v_roms[:,:,0].shape), mask=True)
        vel_hycom_mean = np.ma.array(np.zeros(v_hycom[:,:,0].shape), mask=True)
        vel_roms_mean = np.ma.array(np.zeros(v_roms[:,:,0].shape), mask=True)

        u_hycom_mean = np.mean(u_hycom,axis=2)
        v_hycom_mean = np.mean(v_hycom,axis=2)
        u_roms_mean = np.mean(u_roms,axis=2)
        v_roms_mean = np.mean(v_roms,axis=2)
        vel_hycom_mean = np.mean(vel_hycom,axis=2)
        vel_roms_mean = np.mean(vel_roms,axis=2)

        if not (os.path.isdir(output_dir+"/"+rodada[rod])):
           os.system("mkdir -p "+output_dir+"/"+rodada[rod])
        filename = "VEL_ROMS_GRID_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+ \
                   data_final.strftime("%Y%m%d")+"_"+dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"
        scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_roms': lat_roms, 'lon_roms': lon_roms, \
                         'u_hycom_mean_'+rodada[rod]: u_hycom_mean, 'v_hycom_mean_'+rodada[rod]: v_hycom_mean, \
                         'u_roms_mean': u_roms_mean, 'v_roms_mean': v_roms_mean, \
                         'vel_hycom_mean_'+rodada[rod]: vel_hycom_mean, 'vel_roms_mean': vel_roms_mean})

        del u_hycom_mean,v_hycom_mean,vel_hycom_mean,vel_hycom,u_hycom,v_hycom,filename
        del u_roms_mean, v_roms_mean, vel_roms_mean, lat_roms, lon_roms,u_roms,v_roms,vel_roms,cond_lat_roms,cond_lon_roms

        print()
        print("EXPT: "+rodada[rod]+" PROF:"+str(prof[k,0])+"-"+str(prof[k,1])+" ---- OK")
        print()

exit()
