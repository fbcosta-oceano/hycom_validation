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

data_inicial = "01/01/2011"
data_final = "10/01/2011"
	
woa_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, \
	      85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, \
	      400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, \
	      1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, \
	      1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, \
	      2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, \
	      3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, \
	      4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]
woa_levels = np.array(woa_levels)
dir_obs_roms = '/prj/rodas/leonardo.pires/ROMS'
dir_obs_roms_out = '/prj/rodas/filipe.costa/ROMS/mma'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

current_data = data_inicial
counter = -1


while(current_data <= data_final):
    print(current_data.strftime("%d-%m-%Y"))
    if (current_data==data_inicial):
    
       arq_roms1 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_1de2.nc4'
       arq_roms2 = dir_obs_roms+'/ROMS_'+current_data.strftime("%Y")+'/roms_avg_2de2.nc4'
       if not (os.path.isfile(arq_roms1)) or not (os.path.isfile(arq_roms2)):
          print('ROMS FILE FOR YEAR: '+current_data.strftime("%Y")+' DOES NOT EXIST')
          print(arq_roms1)
          print(arq_roms2)
          exit()
       nc1 = Dataset(arq_roms1, 'r')
       nc2 = Dataset(arq_roms2, 'r')
       lat_u_roms = nc1.variables['lat_u'][:,0] #564
       lon_u_roms = nc1.variables['lon_u'][0,:] #502
       lat_v_roms = nc1.variables['lat_v'][:,0] #564
       lon_v_roms = nc1.variables['lon_v'][0,:] #502
       lat_roms = nc1.variables['lat_rho'][:,0] #564
       lon_roms = nc1.variables['lon_rho'][0,:] #502
       time_roms1 = nc1.variables['ocean_time'][:]
       time_roms2 = nc2.variables['ocean_time'][:]
       time_roms = np.concatenate((time_roms1,time_roms2),axis=0)
       time_roms = time_roms/(60*60*24) #AJUSTANDO PARA DIAS
       time_roms = time_roms-.5 #SUBTRAINDO 12 h PARA FICAR EM 00Z
       ref = datetime.datetime(*time.strptime('01/01/1970',"%d/%m/%Y")[0:3])
       date_roms = [ref + datetime.timedelta(days=time_roms[i]) for i in range(0,len(time_roms),1)]
       ind_inic_date_roms = np.where(date_roms == np.datetime64(datetime.datetime(data_inicial.year, data_inicial.month, data_inicial.day, 0, 0)))
       ind_fim_date_roms = np.where(date_roms == np.datetime64(datetime.datetime(data_final.year, data_final.month, data_final.day, 0, 0)))
    
       print('CARREGANDO T')
       temp_roms = nc1.variables['temp'][0,:,:,:] #(238, 32, 564, 502)
       temp1 = nc1.variables['temp'][:] #(238, 32, 564, 502)
       temp2 = nc2.variables['temp'][:]
       temp = np.ma.concatenate((temp1,temp2),axis=0)
       temp = np.flip(temp,axis=1)
       temp = np.moveaxis(temp,0,-1)
       temp = np.moveaxis(temp,0,-1)
       temp_roms = temp[:,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1,:]
       del temp1, temp2, temp
    
       print('CARREGANDO S')
       salt1 = nc1.variables['salt'][:] #(238, 32, 564, 502)
       salt2 = nc2.variables['salt'][:]
       salt = np.ma.concatenate((salt1,salt2),axis=0)
       salt = np.moveaxis(salt,0,-1)
       salt = np.moveaxis(salt,0,-1) #(564, 502, 365, 32)
       salt_roms = salt[:,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1,:]
       del salt1, salt2, salt
    
       print('CARREGANDO U')
       u1 = nc1.variables['u'][:] #(238, 32, 564, 502)
       u2 = nc2.variables['u'][:]
       u = np.ma.concatenate((u1,u2),axis=0)
       u = np.moveaxis(u,0,-1)
       u = np.moveaxis(u,0,-1) #(564, 502, 365, 32)
       u = u[:,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1,:]
       field = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(u[0,0,:,0])*len(u[0,0,0,:])), mask=True)
       field = field.reshape(len(lat_roms),len(lon_roms),len(u[0,0,:,0]),len(u[0,0,0,:]))
       for i in range(0,len(u[0,0,:,0]),1):
           for j in range(0,len(u[0,0,0,:]),1):
               I = interpolate.interp2d(lon_u_roms,lat_u_roms,u[:,:,i,j],kind='linear')
               field[:,:,i,j] = I(lon_roms,lat_roms)
               del I
       cond = abs(field) > 999999
       u_roms = np.ma.masked_where(cond,field)
       del u1, u2, u, lon_u_roms, lat_u_roms, cond

       print('CARREGANDO V')
       v1 = nc1.variables['v'][:] #(238, 32, 564, 502)
       v2 = nc2.variables['v'][:]
       v = np.ma.concatenate((v1,v2),axis=0)
       v = np.moveaxis(v,0,-1)
       v = np.moveaxis(v,0,-1) #(564, 502, 365, 32)
       v = v[:,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1,:]
       field = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(v[0,0,:,0])*len(v[0,0,0,:])), mask=True)
       field = field.reshape(len(lat_roms),len(lon_roms),len(v[0,0,:,0]),len(v[0,0,0,:]))
       for i in range(0,len(v[0,0,:,0]),1):
           for j in range(0,len(v[0,0,0,:]),1):
               I = interpolate.interp2d(lon_v_roms,lat_v_roms,v[:,:,i,j],kind='linear')
               field[:,:,i,j] = I(lon_roms,lat_roms)
               del I
       cond = abs(field) > 999999
       v_roms = np.ma.masked_where(cond,field)
       del v1, v2, v, lon_v_roms, lat_v_roms, cond

       C  = nc1.variables['Cs_w'][:].data #(33)
       hc = nc1.variables['hc'][:].data
       h  = nc1.variables['h'][:].data #(564, 502)
       sigma_w  = nc1.variables['s_w'][:].data #(33)
       zeta1  = nc1.variables['zeta'][:].data #(238, 564, 502)
       zeta2  = nc2.variables['zeta'][:].data
       zeta = np.ma.concatenate((zeta1,zeta2),axis=0)
    
       S = hc*sigma_w[:,None,None] + np.dot(C[:,None],(h-hc)[:,None,:])
       z = S + zeta[:,None,:,:] *   (1 + (S/h))
       #z.shape (365, 33, 564, 502)
       cond = abs(z)>=999999
       z = np.ma.masked_where(cond,z)
       z = np.flip(z,axis=1)*-1
       z_roms = np.moveaxis(z,0,-1)
       z_roms = np.moveaxis(z_roms,0,-1) #(564, 502, 365, 33)
       z_roms = z_roms[:,:,int(ind_inic_date_roms[0]):int(ind_fim_date_roms[0])+1,1:33]
       z_roms[:,:,:,0] = 0
       del C, hc, sigma_w, zeta1, zeta2, S, z, cond
       nc1.close()
       nc2.close()
       del time_roms1, time_roms2, nc1, nc2

    counter = counter+1
    temp = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(woa_levels)), mask=True)
    temp = temp.reshape(len(lat_roms),len(lon_roms),len(woa_levels))
    salt = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(woa_levels)), mask=True)
    salt = salt.reshape(len(lat_roms),len(lon_roms),len(woa_levels))
    u = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(woa_levels)), mask=True)
    u = u.reshape(len(lat_roms),len(lon_roms),len(woa_levels))
    v = np.ma.array(np.zeros(len(lat_roms)*len(lon_roms)*len(woa_levels)), mask=True)
    v = v.reshape(len(lat_roms),len(lon_roms),len(woa_levels))
    print('INTERPOLANDO PARA NIVEIS WOA')
    for lat in range(0,len(lat_roms),1):
        for lon in range(0,len(lon_roms),1):
            cond = woa_levels<=max(z_roms[lat,lon,counter,:])

            I = interpolate.interp1d(z_roms[lat,lon,counter,:],temp_roms[lat,lon,counter,:])
            temp[lat,lon,cond] = I(woa_levels[cond])
            del I

            I = interpolate.interp1d(z_roms[lat,lon,counter,:],salt_roms[lat,lon,counter,:])
            salt[lat,lon,cond] = I(woa_levels[cond])
            del I

            I = interpolate.interp1d(z_roms[lat,lon,counter,:],u_roms[lat,lon,counter,:])
            u[lat,lon,cond] = I(woa_levels[cond])
            del I

            I = interpolate.interp1d(z_roms[lat,lon,counter,:],v_roms[lat,lon,counter,:])
            v[lat,lon,cond] = I(woa_levels[cond])
            del I

    roms = Dataset(dir_obs_roms_out+'/ROMS_'+current_data.strftime("%Y%m%d")+'.nc','w',format='NETCDF3_CLASSIC')
    roms.set_fill_on()
    xi_rho = roms.createDimension('xi_rho',len(lon_roms))
    eta_rho = roms.createDimension('eta_rho',len(lat_roms))
    z_rho = roms.createDimension('z_rho',len(woa_levels))
    
    longitude = roms.createVariable('lon_rho','f',('xi_rho',))
    latitude = roms.createVariable('lat_rho','f',('eta_rho',))
    z = roms.createVariable('depth','f',('z_rho',))
    t = roms.createVariable('temp','f',('eta_rho','xi_rho','z_rho',),fill_value=9.96920996838687e+36)
    s = roms.createVariable('salt','f',('eta_rho','xi_rho','z_rho',),fill_value=9.96920996838687e+36)
    u1 = roms.createVariable('u','f',('eta_rho','xi_rho','z_rho',),fill_value=9.96920996838687e+36)
    v1 = roms.createVariable('v','f',('eta_rho','xi_rho','z_rho',),fill_value=9.96920996838687e+36)
    ssh = roms.createVariable('zeta','f',('eta_rho','xi_rho',),fill_value=9.96920996838687e+36)
    
    latitude[:] = lat_roms[:]
    longitude[:] = lon_roms[:]
    z[:] = woa_levels[:]
    t[:] = temp[:]
    s[:] = salt[:]
    u1[:] = u[:]
    v1[:] = v[:]
    ssh[:] = zeta[counter,:,:]
    
    roms.close
    del roms, xi_rho, eta_rho, z_rho, longitude, latitude, z, t, s, u, u1, v, v1, ssh

exit()

