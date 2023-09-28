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

woa_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, \
              85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, \
              400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, \
              1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, \
              1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, \
              2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, \
              3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, \
              4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]
woa_levels = np.array(woa_levels)

cfg_file='dirs_calc.txt'
opt=open(cfg_file, "r")
content_list = opt.readlines()
output_dir = str(content_list[0].rstrip('\n'))
dir_hycom = str(content_list[1].rstrip('\n'))
dir_obs = str(content_list[4].rstrip('\n'))

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
    #print(argo.shape)
    #print(argo[-1,:])
    
    del opt, content_list, cont, temp

    current_data = current_data + relativedelta(months=+1)

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
    elif '008d' in rodada[rod]:
       IDM = 1717
       JDM = 2345
       KDM = 32
       dir_ATL = dir_hycom+'/ATLd0.08/'
       record_thkn = 13
       record_temp = 14
       record_saln = 15
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1

    q_buoys = np.ma.array(np.zeros((len(woa_levels))), mask=False)
    rmsd_profile_temp = np.ma.array(np.zeros((len(woa_levels))), mask=False)
    rmsd_profile_saln = np.ma.array(np.zeros((len(woa_levels))), mask=False)
    rmsd_time_temp = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(woa_levels)), mask=False)
    rmsd_time_saln = np.ma.array(np.zeros(((data_final - data_inicial).days+1)*len(woa_levels)), mask=False)
    rmsd_time_temp = np.reshape(rmsd_time_temp,((data_final - data_inicial).days+1,len(woa_levels)))
    rmsd_time_saln = np.reshape(rmsd_time_saln,((data_final - data_inicial).days+1,len(woa_levels)))

    while(current_data <= data_final):
        counter = counter+1
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

           max_lon_hycom = max(lon_hycom)
           min_lon_hycom = min(lon_hycom)
           max_lat_hycom = max(lat_hycom)
           min_lat_hycom = min(lat_hycom)

        arq_dia_hycom = dir_expt+'/output/ab/'+current_data.strftime("%Y%m%d")+'/archv.' \
                        +str((current_data + datetime.timedelta(days=inc_tempo)).timetuple().tm_year).zfill(4)+'_' \
                        +str((current_data + datetime.timedelta(days=inc_tempo)).timetuple().tm_yday).zfill(3)+'_00.a'
        if not (os.path.isfile(arq_dia_hycom)):
           print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           rmsd_time_temp[counter,:].mask=True
           rmsd_time_saln[counter,:].mask=True
           current_data = current_data + datetime.timedelta(days=inc_tempo)
           continue

        if (current_data==data_inicial):
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

           f = open(arq_dia_hycom,'rb')
           for k in range(0,KDM,1):
               f.seek((record_thkn+5*k)*4*(IJDM+npad))
               field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
               field = np.reshape(field,(JDM,IDM))/9806
               cond = abs(field) > 999999999
               field = np.ma.masked_where(cond,field)
               if 'thkn' in locals():
                   thkn = np.ma.concatenate((thkn,field[:,:,np.newaxis]),axis=2)
               else:
                   thkn = field[:,:,np.newaxis]
               del field, cond

               f.seek((record_temp+5*k)*4*(IJDM+npad))
               field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
               field = np.reshape(field,(JDM,IDM))
               cond = abs(field) > 999999999
               field = np.ma.masked_where(cond,field)
               if 'temp' in locals():
                   temp = np.ma.concatenate((temp,field[:,:,np.newaxis]),axis=2)
               else:
                   temp = field[:,:,np.newaxis]
               del field, cond

               f.seek((record_saln+5*k)*4*(IJDM+npad))
               field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
               field = np.reshape(field,(JDM,IDM))
               cond = abs(field) > 999999999
               field = np.ma.masked_where(cond,field)
               if 'saln' in locals():
                   saln = np.ma.concatenate((saln,field[:,:,np.newaxis]),axis=2)
               else:
                   saln = field[:,:,np.newaxis]
               del field, cond

           dif_lon = np.diff(argo_dia[:,3])
           dif_lat = np.diff(argo_dia[:,4])
           cond = (dif_lon!=0) | (dif_lat!=0)
           ind = [i for i, val in enumerate(cond) if val]
           ind.insert(0,-1)
           q_buoys_dia = np.ma.array(np.zeros((len(woa_levels))), mask=False)
           temp_interp_layer = np.ma.array(np.zeros((len(ind),KDM)), mask=True)
           saln_interp_layer = np.ma.array(np.zeros((len(ind),KDM)), mask=True)
           thkn_interp_layer = np.ma.array(np.zeros((len(ind),KDM)), mask=True)
           pres_interp_layer = np.ma.array(np.zeros((len(ind),KDM)), mask=True)

           temp_interp_pres = np.ma.array(np.zeros((len(ind),len(woa_levels))), mask=True)
           saln_interp_pres = np.ma.array(np.zeros((len(ind),len(woa_levels))), mask=True)
           for i in range(0,len(ind),1):
               for k in range(0,KDM,1):
                   I = interpolate.interp2d(lon_hycom,lat_hycom,temp[:,:,k],kind='linear')
                   temp_interp_layer[i,k] = I(argo_dia[ind[i]+1,3] ,argo_dia[ind[i]+1,4])

                   I = interpolate.interp2d(lon_hycom,lat_hycom,saln[:,:,k],kind='linear')
                   saln_interp_layer[i,k] = I(argo_dia[ind[i]+1,3] ,argo_dia[ind[i]+1,4])

                   I = interpolate.interp2d(lon_hycom,lat_hycom,thkn[:,:,k],kind='linear')
                   thkn_interp_layer[i,k] = I(argo_dia[ind[i]+1,3] ,argo_dia[ind[i]+1,4])
                   if (k==0):
                      #pres_interp_layer[i,k] = thkn_interp_layer[i,k]/2
                      pres_interp_layer[i,k] = 0
                   else:
                      pres_interp_layer[i,k] = sum(thkn_interp_layer[i,0:k])+thkn_interp_layer[i,k]/2

               niv_max_hycom = max(pres_interp_layer[i,:])
               if i==len(ind)-1:
                  niv_max_argo = max(argo_dia[ind[i]+1:,0])
                  niv_min = min(argo_dia[ind[i]+1:,0])
                  cond = (woa_levels<=niv_max_argo) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom)
                  I = interpolate.interp1d(argo_dia[ind[i]+1:,0],argo_dia[ind[i]+1:,1])
                  argo_interp_temp = I(woa_levels[cond])
                  I = interpolate.interp1d(argo_dia[ind[i]+1:,0],argo_dia[ind[i]+1:,2])
                  argo_interp_saln = I(woa_levels[cond])
               else:
                  niv_max_argo = max(argo_dia[ind[i]+1:ind[i+1]+1,0])
                  niv_min = min(argo_dia[ind[i]+1:ind[i+1]+1,0])
                  cond = (woa_levels<=niv_max_argo) & (woa_levels>=niv_min) & (woa_levels<=niv_max_hycom)
                  I = interpolate.interp1d(argo_dia[ind[i]+1:ind[i+1]+1,0],argo_dia[ind[i]+1:ind[i+1]+1,1])
                  argo_interp_temp = I(woa_levels[cond])
                  I = interpolate.interp1d(argo_dia[ind[i]+1:ind[i+1]+1,0],argo_dia[ind[i]+1:ind[i+1]+1,2])
                  argo_interp_saln = I(woa_levels[cond])
               I = interpolate.interp1d(pres_interp_layer[i,:],temp_interp_layer[i,:])
               temp_interp_pres = I(woa_levels[cond])
               I = interpolate.interp1d(pres_interp_layer[i,:],saln_interp_layer[i,:])
               saln_interp_pres = I(woa_levels[cond])
               ind_levs = [i for i, val in enumerate(cond) if val]
               cond = abs(temp_interp_pres) > 999999999
               temp_interp_pres = np.ma.masked_where(cond,temp_interp_pres)
               cond = abs(saln_interp_pres) > 999999999
               saln_interp_pres = np.ma.masked_where(cond,saln_interp_pres)
               rmsd_profile_temp[ind_levs] = rmsd_profile_temp[ind_levs] + (temp_interp_pres-argo_interp_temp)**2
               rmsd_profile_saln[ind_levs] = rmsd_profile_saln[ind_levs] + (saln_interp_pres-argo_interp_saln)**2
               rmsd_time_temp[counter,ind_levs] = rmsd_time_temp[counter,ind_levs] + (temp_interp_pres-argo_interp_temp)**2
               rmsd_time_saln[counter,ind_levs] = rmsd_time_saln[counter,ind_levs] + (saln_interp_pres-argo_interp_saln)**2
               if temp_interp_pres.any():
                  q_buoys[ind_levs] = q_buoys[ind_levs]+1
                  q_buoys_dia[ind_levs] = q_buoys_dia[ind_levs]+1
           cond = q_buoys_dia!=0
           rmsd_time_temp[counter,cond] = np.sqrt(np.divide(rmsd_time_temp[counter,cond],q_buoys_dia[cond]))
           rmsd_time_saln[counter,cond] = np.sqrt(np.divide(rmsd_time_saln[counter,cond],q_buoys_dia[cond]))
           cond = q_buoys_dia==0
           rmsd_time_temp[counter,:] = np.ma.masked_where(cond,rmsd_time_temp[counter,:])
           rmsd_time_saln[counter,:] = np.ma.masked_where(cond,rmsd_time_saln[counter,:])

           del temp, saln, temp_interp_pres, saln_interp_pres, q_buoys_dia, cond, ind_levs, temp_interp_layer, saln_interp_layer
           del thkn_interp_layer, pres_interp_layer, argo_interp_temp, argo_interp_saln

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
    dominio_str = f"{abs(min_lat_hycom):.2f}"+s1+"-"+f"{abs(max_lat_hycom):.2f}"+s2+"_"+ \
                  f"{abs(min_lon_hycom):.2f}"+s3+"-"+f"{abs(max_lon_hycom):.2f}"+s4
    filename = "RMSD_ARGO_HYCOM_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+dominio_str
    if not (os.path.isdir(output_dir+"/"+rodada[rod])):
       os.system("mkdir -p "+output_dir+"/"+rodada[rod])
    scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'rmsd_profile_temp_'+rodada[rod]: rmsd_profile_temp, \
                     'rmsd_profile_saln_'+rodada[rod]: rmsd_profile_saln, 'rmsd_time_temp_'+rodada[rod]: rmsd_time_temp, \
                     'rmsd_time_saln_'+rodada[rod]: rmsd_time_saln, 'q_buoys': q_buoys, 'levels': woa_levels})
    print(rmsd_profile_temp)
    del rmsd_profile_saln, rmsd_profile_temp, rmsd_time_temp, rmsd_time_saln, q_buoys

    print()
    print("EXPT: "+rodada[rod]+" ---- OK")
    print()

exit()

