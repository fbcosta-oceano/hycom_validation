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
inc_assim = float(content_list[6].rstrip('\n'))

output_dir = '/home/filipe.costa/resultados'
dir_hycom = '/home/filipe.costa/previsao/hycom_2_2_18/proc'

prof = np.array([[0,0], [0,50], [50,50], [100,100], [200,200], [300,300], [400,400], [800,800], [2000,2000], [3000,3000], [4000,4000]])
#prof = np.array([[0,0], [0,50], [50,50], [100,100], [200,200], [300,300], [400,400], [800,800]])
#prof = np.array([[0,0], [100,100], [400,400], [800,800]])
k_inic = 0

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    for k in range(k_inic,len(prof[:,0]),1):
        current_data = data_inicial
        if '004j' in rodada[rod]:
           dir_ATL = dir_hycom+'/ATLj0.04/'
        elif '008d' in rodada[rod]:
           IDM = 1717
           JDM = 2345
           dir_ATL = dir_hycom+'/ATLd0.08/'
        elif '008i' in rodada[rod]:
           IDM = 628
           JDM = 780
           dir_ATL = dir_hycom+'/ATLi0.08/'
        else:
           print('CASE '+rodada[rod]+' NOT DEFINED')
           continue

        dir_expt = dir_ATL+expt[rod]
        counter = -1

        while(current_data <= data_final):
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" PROF: "+str(prof[k,0])+"-"+str(prof[k,1]))

            arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                            +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'
            if not (os.path.isfile(arq_dia_hycom)):
               for i in range(2,int(inc_assim)+1,1):
                   arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'00/archv.' \
                                   +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                   +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'
                   if (os.path.isfile(arq_dia_hycom)):
                      break                            
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

               max_lon_hycom = max(lon_hycom)
               min_lon_hycom = min(lon_hycom)
               max_lat_hycom = max(lat_hycom)
               min_lat_hycom = min(lat_hycom)

            nc = Dataset(arq_dia_hycom, 'r')
            u = nc.variables['u'][0,:,:,:] #(33, 564, 502)
            u = np.moveaxis(u,0,-1)
            v = nc.variables['v'][0,:,:,:]
            v = np.moveaxis(v,0,-1)
            depth_hycom_temp = depth_hycom

            if not((depth_hycom_temp == prof[k,0]).any()):
                I = interpolate.interp1d(depth_hycom_temp,v,axis=-1)
                v_interp = I(prof[k,0])
                I = interpolate.interp1d(depth_hycom_temp,u,axis=-1)
                u_interp = I(prof[k,0])
                cond = abs(v_interp) > 999999999
                v_interp = np.ma.masked_where(cond,v_interp)
                u_interp = np.ma.masked_where(cond,v_interp)

                ind = np.squeeze(np.array(np.where(depth_hycom_temp<prof[k,0])),axis=0)
                depth_hycom_temp = np.insert(depth_hycom_temp,ind[-1]+1,prof[k,0])
                v = np.insert(v,ind[-1]+1,v_interp,axis=-1)
                u = np.insert(u,ind[-1]+1,u_interp,axis=-1)
                del I, cond, ind, v_interp, u_interp

            if not((depth_hycom_temp == prof[k,1]).any()):
                I = interpolate.interp1d(depth_hycom_temp,v,axis=-1)
                v_interp = I(prof[k,1])
                I = interpolate.interp1d(depth_hycom_temp,u,axis=-1)
                u_interp = I(prof[k,1])
                cond = abs(v_interp) > 999999999
                v_interp = np.ma.masked_where(cond,v_interp)
                u_interp = np.ma.masked_where(cond,v_interp)

                ind = np.squeeze(np.array(np.where(depth_hycom_temp<prof[k,1])),axis=0)
                depth_hycom_temp = np.insert(depth_hycom_temp,ind[-1]+1,prof[k,1])
                v = np.insert(v,ind[-1]+1,v_interp,axis=-1)
                u = np.insert(u,ind[-1]+1,u_interp,axis=-1)
                del I, cond, ind, v_interp, u_interp

            cond = (depth_hycom_temp <= prof[k,1]) & (depth_hycom_temp >= prof[k,0])
            v = np.mean(v[:,:,cond],axis=-1)
            u = np.mean(u[:,:,cond],axis=-1)
            if 'v_hycom' in locals():
                v_hycom = np.ma.concatenate((v_hycom,v[:,:,np.newaxis]),axis=2)
                u_hycom = np.ma.concatenate((u_hycom,v[:,:,np.newaxis]),axis=2)
            else:
                v_hycom = v[:,:,np.newaxis]
                u_hycom = u[:,:,np.newaxis]

            nc.close()
            del nc,v,u,cond,depth_hycom_temp
            current_data = current_data + datetime.timedelta(days=inc_tempo)

        EC_hycom = np.ma.array(np.zeros(v_hycom.shape), mask=True)
        EC_hycom = (u_hycom**2+v_hycom**2)/2
        del u_hycom, v_hycom

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
            counter_fim = min(counter_fim,len(EC_hycom[0,0,:]))

            EC_hycom_mean_map = np.ma.array(np.zeros(EC_hycom[:,:,0].shape), mask=True)
            EC_hycom_mean_time = np.ma.array(np.zeros(EC_hycom[0,0,counter_inic:counter_fim].shape), mask=True)

            EC_hycom_mean_map = np.mean(EC_hycom[:,:,counter_inic:counter_fim],axis=2)
            EC_hycom_mean_time = np.mean(EC_hycom[:,:,counter_inic:counter_fim],axis=(0,1))

            if not (os.path.isdir(output_dir+"/"+rodada[rod])):
               os.system("mkdir -p "+output_dir+"/"+rodada[rod])
            filename = "EC_HYCOM_"+rodada[rod]+"_"+current_data.strftime("%Y%m%d")+"-"+ \
                       (current_data + datetime.timedelta(days=counter_fim-counter_inic-1)).strftime("%Y%m%d")+"_"+ \
                       dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"
            scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_hycom': lat_hycom, 'lon_hycom': lon_hycom, \
                             'EC_hycom_mean_map_'+rodada[rod]: EC_hycom_mean_map, 'EC_hycom_mean_time_'+rodada[rod]: EC_hycom_mean_time})

            del EC_hycom_mean_map, EC_hycom_mean_time

            counter_inic = counter_fim
            current_data = current_data + datetime.timedelta(days=days)

        EC_hycom_mean_map = np.ma.array(np.zeros(EC_hycom[:,:,0].shape), mask=True)
        EC_hycom_mean_time = np.ma.array(np.zeros(EC_hycom[0,0,:].shape), mask=True)

        EC_hycom_mean_map = np.mean(EC_hycom,axis=2)
        EC_hycom_mean_time = np.mean(EC_hycom,axis=(0,1))

        if not (os.path.isdir(output_dir+"/"+rodada[rod])):
           os.system("mkdir -p "+output_dir+"/"+rodada[rod])
        filename = "EC_HYCOM_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+ \
                   data_final.strftime("%Y%m%d")+"_"+dominio_str+"_"+str(prof[k,0])+"-"+str(prof[k,1])+"m"
        scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'lat_hycom': lat_hycom, 'lon_hycom': lon_hycom, \
                         'EC_hycom_mean_map_'+rodada[rod]: EC_hycom_mean_map, 'EC_hycom_mean_time_'+rodada[rod]: EC_hycom_mean_time})

        del EC_hycom_mean_map,EC_hycom_mean_time,filename
        del lat_hycom, lon_hycom

        print()
        print("EXPT: "+rodada[rod]+" PROF:"+str(prof[k,0])+"-"+str(prof[k,1])+" ---- OK")
        print()

exit()
