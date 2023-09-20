from netCDF4 import Dataset
from scipy import interpolate
import scipy.io
import remo
import locale
import datetime
from dateutil.relativedelta import relativedelta
import sys
import os
import time 
import math
import pandas as pd
import numpy as np

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
inc_assim = float(content_list[6].rstrip('\n'))

woa_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, \
              85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, \
              400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, \
              1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, \
              1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, \
              2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, \
              3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, \
              4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]
#print(type(woa_levels))
woa_levels = np.array(woa_levels)
#print(type(woa_levels))
output_dir = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/resultados'
dir_hycom = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/previsao/hycom_2_2_18/proc'
dir_obs = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/dados_obs/ARGO/controle_qualidade_argo'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

#locale.setlocale(locale.LC_ALL, 'pt_BR.utf8')
q_mes = (data_final.year - data_inicial.year) * 12 + data_final.month - data_inicial.month + 1
current_data = data_inicial
for m in range(0,q_mes,1):
    arq_obs = dir_obs+'/'+current_data.strftime("%Y")+'/'+current_data.strftime("%B")+'/theta_sal_argo_' \
              +current_data.strftime("%Y")+current_data.strftime("%m")+'_formatted.ascii'

    opt=open(arq_obs, "r")
    content_list = opt.readlines()
    temp = np.ma.array(np.zeros((len(content_list),8)), mask=True)
    for i in range(0,len(content_list),1):
        cont = content_list[i].rstrip('\n').split(' ')
        cont = list(filter(None,cont))
        cont = np.array(cont)
        temp[i] = cont[:]
    if m==0:
       argo = temp
    else:
       argo = np.ma.concatenate((argo,temp),axis=0)
    
    del opt, content_list, cont, temp

    current_data = current_data + relativedelta(months=+1)

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       dir_ATL = dir_hycom+'/ATLj0.04/'

    elif '008i' in rodada[rod]:
       dir_ATL = dir_hycom+'/ATLi0.08/'
       whole_domain = True
       ### [min_lon, max_lon, min_lat, max_lat]
       subs = np.array([[-60, -20, -36,   7], \
                        [-53, -33, -33, -13], \
                        [-51, -18, -30, -10], \
                        [-67, -40, -47, -34], \
                        [-36,  -5, -47, -34]])

    elif '008d' in rodada[rod]:
       dir_ATL = dir_hycom+'/ATLd0.08/'
       whole_domain = True
       ### [min_lon, max_lon, min_lat, max_lat]  ### [ALTi0.08], [METAREA V], [ATLj0.04], [SUB A], [SUB B].... [SUB M]
       subs = np.array([[-68, -18, -45,  10], \
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

    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    for region in range(0,subs.shape[0]+1):
        if not whole_domain and region==subs.shape[0]:
           continue 
        counter = -1

        q_buoys = np.ma.array(np.zeros((len(woa_levels))), mask=False)
        rmsd_profile_temp = np.ma.array(np.zeros((len(woa_levels))), mask=False)
        rmsd_profile_saln = np.ma.array(np.zeros((len(woa_levels))), mask=False)
        rmsd_time_temp = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(woa_levels)), mask=False)
        rmsd_time_saln = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(woa_levels)), mask=False)
        rmsd_time_temp = np.reshape(rmsd_time_temp,((data_final - data_inicial).days+1,len(woa_levels)))
        rmsd_time_saln = np.reshape(rmsd_time_saln,((data_final - data_inicial).days+1,len(woa_levels)))

        current_data = data_inicial
        while(current_data <= data_final):
            counter = counter+1
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y"))

            arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                            +str(current_data.timetuple().tm_year).zfill(4)+'_' \
                            +str(current_data.timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            if not (os.path.isfile(arq_dia_hycom)):
               for i in range(2,int(inc_assim)+1,1):
                   arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'00/archv.' \
                                   +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                   +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
                   if (os.path.isfile(arq_dia_hycom)):
                      break            
            if not (os.path.isfile(arq_dia_hycom)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               rmsd_time_temp[counter,:].mask=True
               rmsd_time_saln[counter,:].mask=True
               current_data = current_data + datetime.timedelta(days=inc_tempo)
               continue

            if (current_data==data_inicial):
               nc = Dataset(arq_dia_hycom, 'r')
               lon_hycom = nc.variables['Longitude'][:]
               lat_hycom = nc.variables['Latitude'][:]
               depth_hycom = nc.variables['Depth'][:]
               if region != subs.shape[0]:
                  cond_lat_hycom = (lat_hycom<=subs[region,3]) & (lat_hycom>=subs[region,2])
                  lat_hycom = lat_hycom[cond_lat_hycom]
                  cond_lon_hycom = (lon_hycom<=subs[region,1]) & (lon_hycom>=subs[region,0])
                  lon_hycom = lon_hycom[cond_lon_hycom]

               max_lon_hycom = max(lon_hycom)
               min_lon_hycom = min(lon_hycom)
               max_lat_hycom = max(lat_hycom)
               min_lat_hycom = min(lat_hycom)
               nc.close()
               del nc

               cond = (argo[:,3] < max_lon_hycom) & (argo[:,3] > min_lon_hycom) & \
                      (argo[:,4] < max_lat_hycom) & (argo[:,4] > min_lat_hycom)
               argo_work = argo[cond,:]
               del cond

            #del thkn, temp, saln

            cond = (argo_work[:,5] == float(current_data.strftime("%d"))) & (argo_work[:,6] == float(current_data.strftime("%m"))) & \
                   (argo_work[:,7] == float(current_data.strftime("%Y")))
            argo_dia = argo_work[cond,:]
            if argo_dia.any():
               print('ARGO ENCONTRADO')

               nc = Dataset(arq_dia_hycom, 'r')
               temp = nc.variables['temperature'][0,:,:,:] #(KDM, JDM, IDM)
               temp = np.moveaxis(temp,0,-1)
               saln = nc.variables['salinity'][0,:,:,:] #(KDM, JDM, IDM)
               saln = np.moveaxis(saln,0,-1)
               if region != subs.shape[0]:
                   temp = temp[cond_lat_hycom,:,:]
                   temp = temp[:,cond_lon_hycom,:]
                   saln = saln[cond_lat_hycom,:,:]
                   saln = saln[:,cond_lon_hycom,:]
               nc.close()
               del nc

               dif_lon = np.diff(argo_dia[:,3])
               dif_lat = np.diff(argo_dia[:,4])
               cond = (dif_lon!=0) | (dif_lat!=0)
               ind = [i for i, val in enumerate(cond) if val]
               ind.insert(0,-1)
               q_buoys_dia = np.ma.array(np.zeros((len(woa_levels))), mask=False)

               for i in range(0,len(ind),1):
                   temp_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                   saln_interp = np.ma.array(np.zeros(len(depth_hycom)), mask=True)
                   temp_interp_woa = np.ma.array(np.zeros(len(woa_levels)), mask=True)
                   saln_interp_woa = np.ma.array(np.zeros(len(woa_levels)), mask=True)
                   for k in range(0,len(depth_hycom),1):
                       I = interpolate.interp2d(lon_hycom,lat_hycom,temp[:,:,k],kind='linear')
                       temp_interp[k] = I(argo_dia[ind[i]+1,3] ,argo_dia[ind[i]+1,4])

                       I = interpolate.interp2d(lon_hycom,lat_hycom,saln[:,:,k],kind='linear')
                       saln_interp[k] = I(argo_dia[ind[i]+1,3] ,argo_dia[ind[i]+1,4])

                   niv_max_hycom = max(depth_hycom)

                   I = interpolate.interp1d(depth_hycom,temp_interp[:])
                   temp_interp_woa[:] = I(woa_levels)
                   I = interpolate.interp1d(depth_hycom,saln_interp[:])
                   saln_interp_woa[:] = I(woa_levels)
                   if i==len(ind)-1:
                      niv_max_argo = max(argo_dia[ind[i]+1:,0])
                      niv_min = min(argo_dia[ind[i]+1:,0])
                      cond = (woa_levels<=niv_max_argo) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom) & \
                             (abs(temp_interp_woa) < 99) & (abs(saln_interp_woa) < 99)
                      I = interpolate.interp1d(argo_dia[ind[i]+1:,0],argo_dia[ind[i]+1:,1])
                      argo_interp_temp = I(woa_levels[cond])
                      I = interpolate.interp1d(argo_dia[ind[i]+1:,0],argo_dia[ind[i]+1:,2])
                      argo_interp_saln = I(woa_levels[cond])
                   else:
                      niv_max_argo = max(argo_dia[ind[i]+1:ind[i+1]+1,0])
                      niv_min = min(argo_dia[ind[i]+1:ind[i+1]+1,0])
                      cond = (woa_levels<=niv_max_argo) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom) & \
                             (abs(temp_interp_woa) < 99) & (abs(saln_interp_woa) < 99)
                      I = interpolate.interp1d(argo_dia[ind[i]+1:ind[i+1]+1,0],argo_dia[ind[i]+1:ind[i+1]+1,1])
                      argo_interp_temp = I(woa_levels[cond])
                      I = interpolate.interp1d(argo_dia[ind[i]+1:ind[i+1]+1,0],argo_dia[ind[i]+1:ind[i+1]+1,2])
                      argo_interp_saln = I(woa_levels[cond])
                   ind_levs = [i for i, val in enumerate(cond) if val]

                   rmsd_profile_temp[ind_levs] = rmsd_profile_temp[ind_levs] + (temp_interp_woa[ind_levs]-argo_interp_temp)**2
                   rmsd_profile_saln[ind_levs] = rmsd_profile_saln[ind_levs] + (saln_interp_woa[ind_levs]-argo_interp_saln)**2
                   rmsd_time_temp[counter,ind_levs] = rmsd_time_temp[counter,ind_levs] + (temp_interp_woa[ind_levs]-argo_interp_temp)**2
                   rmsd_time_saln[counter,ind_levs] = rmsd_time_saln[counter,ind_levs] + (saln_interp_woa[ind_levs]-argo_interp_saln)**2
                   if temp_interp_woa.any():
                      q_buoys[ind_levs] = q_buoys[ind_levs]+1
                      q_buoys_dia[ind_levs] = q_buoys_dia[ind_levs]+1
                   del temp_interp_woa, saln_interp_woa
               cond = q_buoys_dia!=0
               rmsd_time_temp[counter,cond] = np.sqrt(np.divide(rmsd_time_temp[counter,cond],q_buoys_dia[cond]))
               rmsd_time_saln[counter,cond] = np.sqrt(np.divide(rmsd_time_saln[counter,cond],q_buoys_dia[cond]))
               cond = q_buoys_dia==0
               rmsd_time_temp[counter,:] = np.ma.masked_where(cond,rmsd_time_temp[counter,:])
               rmsd_time_saln[counter,:] = np.ma.masked_where(cond,rmsd_time_saln[counter,:])

               del temp, saln, q_buoys_dia, cond, ind_levs, temp_interp, saln_interp
               del argo_interp_temp, argo_interp_saln

            else:
               rmsd_time_temp[counter,:].mask=True
               rmsd_time_saln[counter,:].mask=True
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" --OK")
            current_data = current_data + datetime.timedelta(days=inc_tempo)
        cond = q_buoys!=0
        rmsd_profile_temp[cond] = np.sqrt(np.divide(rmsd_profile_temp[cond],q_buoys[cond]))
        rmsd_profile_saln[cond] = np.sqrt(np.divide(rmsd_profile_saln[cond],q_buoys[cond]))
        cond = q_buoys==0
        rmsd_profile_temp = np.ma.masked_where(cond,rmsd_profile_temp)
        rmsd_profile_saln = np.ma.masked_where(cond,rmsd_profile_saln)

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
        if whole_domain and region==subs.shape[0]:   
           dominio_str = f"{abs(min_lat_hycom):.2f}"+s1+"-"+f"{abs(max_lat_hycom):.2f}"+s2+"_"+ \
                         f"{abs(min_lon_hycom):.2f}"+s3+"-"+f"{abs(max_lon_hycom):.2f}"+s4
        else:
           dominio_str = f"{abs(subs[region,2]):.2f}"+s1+"-"+f"{abs(subs[region,3]):.2f}"+s2+"_"+ \
                         f"{abs(subs[region,0]):.2f}"+s3+"-"+f"{abs(subs[region,1]):.2f}"+s4
        filename = "RMSD_ARGO_HYCOM_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
        if not (os.path.isdir(output_dir+"/"+rodada[rod])):
           os.system("mkdir -p "+output_dir+"/"+rodada[rod])
        scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'rmsd_profile_temp_'+rodada[rod]: rmsd_profile_temp, \
                         'rmsd_profile_saln_'+rodada[rod]: rmsd_profile_saln, 'rmsd_time_temp_'+rodada[rod]: rmsd_time_temp, \
                         'rmsd_time_saln_'+rodada[rod]: rmsd_time_saln, 'q_buoys': q_buoys, 'levels': woa_levels})
        del rmsd_profile_saln, rmsd_profile_temp, rmsd_time_temp, rmsd_time_saln, q_buoys

        print()
        print("EXPT: "+rodada[rod]+"  REGION: "+str(region)+" ---- OK")
        print()

exit()

