import warnings
import datetime
import sys
import os
import time
#from netCDF4 import nc.Dataset
import netCDF4 as nc
import matplotlib
import geopy.distance
from scipy import interpolate
import scipy.io
import numpy as np
warnings.filterwarnings('ignore')
jet = matplotlib.cm.get_cmap('jet')

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))
inc_assim = float(content_list[6].rstrip('\n'))

output_dir = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/resultados'
dir_hycom = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/previsao/hycom_2_2_18/proc'
sobs = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/dados_obs/ARGO/controle_qualidade_argo'

lats = np.array([-9, -10, -11, -13, -14, -15, -16, -17, -18, -19, -21, -22, -23,
                 -24, -25, -26, -27, -28, -29, -30, -32, -44])
#lats = np.array([-24, -25, -26, -27, -28, -29, -30, -32, -44])
#lats = np.array([-9, -10, -22])
#print(lats)


data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
days = (data_final-data_inicial).days +1
months = (data_final.year - data_inicial.year)*12 + data_final.month - data_inicial.month +1

for ll in range(0,len(lats),1):
    if lats[ll]==-9:
        p1 = [-08.80, -35.20] #[lat lon]
        p2 = [-09.50, -32.50] #[lat lon]
    elif lats[ll]==-10:
        p1 = [-09.60, -35.50] #[lat lon]
        p2 = [-11.10, -33.00] #[lat lon]
    elif lats[ll]==-11:
        p1 = [-10.50, -36.50] #[lat lon]
        p2 = [-12.00, -34.30] #[lat lon]
    elif lats[ll]==-13:
        p1 = [-12.60, -38.20] #[lat lon]
        p2 = [-13.60, -36.50] #[lat lon]
    elif lats[ll]==-14:
        p1 = [-13.70, -39.10] #[lat lon]
        p2 = [-14.90, -37.00] #[lat lon]
    elif lats[ll]==-15:
        p1 = [-15.00, -39.10] #[lat lon]
        p2 = [-15.00, -37.00] #[lat lon]
    elif lats[ll]==-16:
        p1 = [-16.50, -39.00] #[lat lon]
        p2 = [-15.00, -36.00] #[lat lon]
    elif lats[ll]==-17:
        p1 = [-16.35, -39.20] #[lat lon]
        p2 = [-18.20, -36.00] #[lat lon]
    elif lats[ll]==-18:
       p1 = [-18.90, -39.80] #[lat lon]
       p2 = [-17.00, -35.50] #[lat lon]
    elif lats[ll]==-19:
       p1 = [-18.30, -39.80] #[lat lon]
       p2 = [-19.60, -35.50] #[lat lon]
    elif lats[ll]==-21:
       p1 = [-20.75, -40.70] #[lat lon]
       p2 = [-21.20, -39.00] #[lat lon]
    elif lats[ll]==-22:
       p1 = [-22.00, -41.10] #[lat lon]
       p2 = [-22.00, -39.00] #[lat lon]
    elif lats[ll]==-23:
       p1 = [-22.25, -42.00] #[lat lon]
       p2 = [-23.70, -39.00] #[lat lon]
    elif lats[ll]==-24:
       p1 = [-22.80, -42.25] #[lat lon]
       p2 = [-25.50, -42.25] #[lat lon]
    elif lats[ll]==-25:
       p1 = [-24.00, -46.00] #[lat lon]
       p2 = [-26.40, -43.00] #[lat lon]
    elif lats[ll]==-26:
       p1 = [-24.60, -47.50] #[lat lon]
       p2 = [-26.70, -44.50] #[lat lon]
    elif lats[ll]==-27:
       p1 = [-26.00, -48.50] #[lat lon]
       p2 = [-28.10, -44.50] #[lat lon]
    elif lats[ll]==-28:
       p1 = [-27.25, -48.70] #[lat lon]
       p2 = [-28.55, -46.00] #[lat lon]
    elif lats[ll]==-29:
       p1 = [-28.35, -49.00] #[lat lon]
       p2 = [-29.50, -46.00] #[lat lon]
    elif lats[ll]==-30:
       p1 = [-29.30, -50.00] #[lat lon]
       p2 = [-30.70, -46.00] #[lat lon]
    elif lats[ll]==-32:
       p1 = [-32.00, -51.00] #[lat lon]
       p2 = [-32.00, -48.00] #[lat lon]
    elif lats[ll]==-44:
       p1 = [-23.00, -44.00] #[lat lon]
       p2 = [-28.00, -44.00] #[lat lon]
       
        

    if p1[0]<0:
        s_lat_p1 = 'S'
    else:
        s_lat_p1 = 'N'

    if p2[0]<0:
        s_lat_p2 = 'S'
    else:
        s_lat_p2 = 'N'

    if p1[1]<0:
        s_lon_p1 = 'W'
    else:
        s_lon_p1 = 'E'

    if p2[1]<0:
        s_lon_p2 = 'W'
    else:
        s_lon_p2 = 'E'


    min_lat = min(p1[0],p2[0])-2
    max_lat = max(p1[0],p2[0])+2
    min_lon = min(p1[1],p2[1])-2
    max_lon = max(p1[1],p2[1])+2
    for i in range(1,13,1):
        if i==1:
            t_woa_arq = sobs+'/woa23_decav91C0_t01_04_seasonal.nc'
            nc_file = nc.Dataset(t_woa_arq, 'r')
            lat_woa = nc_file.variables['lat'][:]
            lon_woa = nc_file.variables['lon'][:]

            #ind_min_lat = min(range(len(lat_woa)), key=lambda x:abs(lat_woa[x]-min_lat))
            #ind_max_lat = min(range(len(lat_woa)), key=lambda x:abs(lat_woa[x]-max_lat))
            #ind_min_lon = min(range(len(lon_woa)), key=lambda x:abs(lon_woa[x]-min_lon))
            #ind_max_lon = min(range(len(lon_woa)), key=lambda x:abs(lon_woa[x]-max_lon))

            cond_lat = (lat_woa<=max_lat) & (lat_woa>=min_lat)
            cond_lon = (lon_woa<=max_lon) & (lon_woa>=min_lon)

            lat_woa = lat_woa[cond_lat]
            lon_woa = lon_woa[cond_lon]

            depth_woa_seas = nc_file.variables['depth'][:]
            nc_file.close()
            temp_woa = np.ma.array(np.zeros(len(depth_woa_seas)*len(lat_woa)*len(lon_woa)*12), mask=True)
            temp_woa = temp_woa.reshape(12,len(lat_woa),len(lon_woa),len(depth_woa_seas))
            sal_woa = np.ma.array(np.zeros(len(depth_woa_seas)*len(lat_woa)*len(lon_woa)*12), mask=True)
            sal_woa = sal_woa.reshape(12,len(lat_woa),len(lon_woa),len(depth_woa_seas))
            t_woa_arq = sobs+'/woa23_decav91C0_t01_04.nc'
            nc_file = nc.Dataset(t_woa_arq, 'r')
            depth_woa = nc_file.variables['depth'][:]
            nc_file.close()

        t_woa_arq = sobs+'/woa23_decav91C0_t'+str(i).zfill(2)+'_04.nc'
        nc_file = nc.Dataset(t_woa_arq, 'r')
        t_an = np.squeeze(nc_file.variables['t_an'][:])
        t_an = np.moveaxis(t_an,0,-1)
        woa = t_an[cond_lat,:,:]
        temp_woa[i-1,:,:,0:len(depth_woa)] = woa[:,cond_lon,:]
        nc_file.close()
        del t_an, t_woa_arq, woa

        t_woa_arq = sobs+'/woa23_decav91C0_t'+str(i).zfill(2)+'_04_seasonal.nc'
        nc_file = nc.Dataset(t_woa_arq, 'r')
        t_an = np.squeeze(nc_file.variables['t_an'][:])
        t_an = np.moveaxis(t_an,0,-1)
        woa = t_an[cond_lat,:,len(depth_woa):-1]
        temp_woa[i-1,:,:,len(depth_woa):-1] = woa[:,cond_lon,:]
        nc_file.close()
        del t_an, t_woa_arq, woa

        s_woa_arq = sobs+'/woa23_decav91C0_s'+str(i).zfill(2)+'_04.nc'
        nc_file = nc.Dataset(s_woa_arq, 'r')
        s_an = np.squeeze(nc_file.variables['s_an'][:])
        s_an = np.moveaxis(s_an,0,-1)
        woa = s_an[cond_lat,:,:]
        sal_woa[i-1,:,:,0:len(depth_woa)] = woa[:,cond_lon,:]
        nc_file.close()
        del s_an, s_woa_arq, woa

        s_woa_arq = sobs+'/woa23_decav91C0_s'+str(i).zfill(2)+'_04_seasonal.nc'
        nc_file = nc.Dataset(s_woa_arq, 'r')
        s_an = np.squeeze(nc_file.variables['s_an'][:])
        s_an = np.moveaxis(s_an,0,-1)
        woa = s_an[cond_lat,:,len(depth_woa):-1]
        sal_woa[i-1,:,:,len(depth_woa):-1] = woa[:,cond_lon,:]
        nc_file.close()
        del s_an, s_woa_arq, woa

        print('WOA MES: '+str(i).zfill(2)+' ----OK')

    for rod in range(0,len(rodada),1):
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
            print(rodada[rod]+" DIA: "+current_data.strftime("%d-%m-%Y")+" LAT: "+str(lats[ll]))

            arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                            +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                            +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'
            if os.path.isfile(arq_dia_hycom):
               arq_dia_hycom2 = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                                +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3di.nc'

               arq_dia_hycom3 = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-1)).strftime("%Y%m%d")+'00/archv.' \
                                +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
            else:
               for i in range(2,int(inc_assim)+1,1):
                   arq_dia_hycom = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'00/archv.' \
                                   +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                   +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zu.nc'
                   if (os.path.isfile(arq_dia_hycom)):
                      arq_dia_hycom2 = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'00/archv.' \
                                       +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                       +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3di.nc'
                      arq_dia_hycom3 = dir_expt+'/output/Ncdf/'+(current_data + datetime.timedelta(days=-i)).strftime("%Y%m%d")+'00/archv.' \
                                       +str((current_data).timetuple().tm_year).zfill(4)+'_' \
                                       +str((current_data).timetuple().tm_yday).zfill(3)+'_00_3zt.nc'
                      break

            if not (os.path.isfile(arq_dia_hycom)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq_dia_hycom)
               exit()

            counter = counter+1
            if (counter==0):
               mes =  current_data.month
               counter_mes = -1
            if (current_data==data_inicial):
                nc_file = nc.Dataset(arq_dia_hycom, 'r')
                lon_hycom = nc_file.variables['Longitude'][:]
                lat_hycom = nc_file.variables['Latitude'][:]
                depth_hycom = nc_file.variables['Depth'][:]
                nc_file.close()
                del nc_file
                nc_file = nc.Dataset(arq_dia_hycom2, 'r')
                layers = nc_file.variables['Layer'][:]
                nc_file.close()
                del nc_file

                cond_lat = (lat_hycom<=max_lat) & (lat_hycom>=min_lat)
                cond_lon = (lon_hycom<=max_lon) & (lon_hycom>=min_lon)

                lat_hycom = lat_hycom[cond_lat]
                lon_hycom = lon_hycom[cond_lon]

                dist_pts = geopy.distance.geodesic(p1, p2).km
                res_mod = (lat_hycom[1]-lat_hycom[0])*111.2;
                if res_mod==0:
                   res_mod = (lon_hycom[1]-lon_hycom[0])*111.2;
                res_pts = round(dist_pts/(res_mod/2.0));

                inc_lat = (p2[0]-p1[0])/res_pts
                inc_lon = (p2[1]-p1[1])/res_pts
                lat_sec = np.ma.array(np.zeros(res_pts), mask=False)
                lon_sec = np.ma.array(np.zeros(res_pts), mask=False)
                lat_sec[0] = p1[0]
                lon_sec[0] = p1[1]
                for i in range(0,res_pts-1,1):
                    lat_sec[i+1] = lat_sec[i]+inc_lat
                    lon_sec[i+1] = lon_sec[i]+inc_lon
                dy = geopy.distance.geodesic([p1[0],p1[1]], [p2[0],p1[1]]).km
                dx = geopy.distance.geodesic([p1[0],p1[1]], [p1[0],p2[1]]).km
                if dx==0:
                   ang = 0 
                else:
                   ang = 180*np.arctan(dy/dx)/np.pi
                theta = np.deg2rad(ang)

                rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
                
                temp = np.ma.array(np.zeros(len(depth_woa_seas)*res_pts*12), mask=True)
                temp = temp.reshape(12,res_pts,len(depth_woa_seas))
                temp_sec_woa_mes = np.ma.array(np.zeros(len(depth_hycom)*res_pts*12), mask=True)
                temp_sec_woa_mes = temp_sec_woa_mes.reshape(12,res_pts,len(depth_hycom))
                temp_sec_woa = np.ma.array(np.zeros(len(depth_hycom)*res_pts), mask=True)
                temp_sec_woa = temp_sec_woa.reshape(res_pts,len(depth_hycom))

                sal = np.ma.array(np.zeros(len(depth_woa_seas)*res_pts*12), mask=True)
                sal = sal.reshape(12,res_pts,len(depth_woa_seas))
                sal_sec_woa_mes = np.ma.array(np.zeros(len(depth_hycom)*res_pts*12), mask=True)
                sal_sec_woa_mes = sal_sec_woa_mes.reshape(12,res_pts,len(depth_hycom))
                sal_sec_woa = np.ma.array(np.zeros(len(depth_hycom)*res_pts), mask=True)
                sal_sec_woa = sal_sec_woa.reshape(res_pts,len(depth_hycom))

                v_vel_reg = np.ma.array(np.zeros(days*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                v_vel_reg = v_vel_reg.reshape(days,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                v_vel_reg_mes = np.ma.array(np.zeros(months*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                v_vel_reg_mes = v_vel_reg_mes.reshape(months,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_v = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_v = s_v.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                v_vel_sec_mes = np.ma.array(np.zeros(months*res_pts*len(depth_hycom)), mask=True)
                v_vel_sec_mes = v_vel_sec_mes.reshape(months,res_pts,len(depth_hycom))
                v_vel_rot = np.ma.array(np.zeros(res_pts*len(depth_hycom)), mask=True)
                v_vel_rot = v_vel_rot.reshape(res_pts,len(depth_hycom))

                u_vel_reg = np.ma.array(np.zeros(days*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                u_vel_reg = u_vel_reg.reshape(days,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                u_vel_reg_mes = np.ma.array(np.zeros(months*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                u_vel_reg_mes = u_vel_reg_mes.reshape(months,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_u = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_u = s_u.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                u_vel_sec_mes = np.ma.array(np.zeros(months*res_pts*len(depth_hycom)), mask=True)
                u_vel_sec_mes = u_vel_sec_mes.reshape(months,res_pts,len(depth_hycom))
                u_vel_rot = np.ma.array(np.zeros(res_pts*len(depth_hycom)), mask=True)
                u_vel_rot = u_vel_rot.reshape(res_pts,len(depth_hycom))

                int_depth_reg = np.ma.array(np.zeros(days*len(lat_hycom)*len(lon_hycom)*(len(depth_hycom)-1)), mask=True)
                int_depth_reg = int_depth_reg.reshape(days,len(lat_hycom),len(lon_hycom),len(depth_hycom)-1)
                int_depth_reg_mes = np.ma.array(np.zeros(months*len(lat_hycom)*len(lon_hycom)*(len(depth_hycom)-1)), mask=True)
                int_depth_reg_mes = int_depth_reg_mes.reshape(months,len(lat_hycom),len(lon_hycom),len(depth_hycom)-1)
                s_int_depth = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*(len(depth_hycom)-1)), mask=False)
                s_int_depth = s_int_depth.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom)-1)
                int_depth_sec_mes = np.ma.array(np.zeros(months*res_pts*(len(depth_hycom)-1)), mask=True)
                int_depth_sec_mes = int_depth_sec_mes.reshape(months,res_pts,len(depth_hycom)-1)
                int_depth_sec = np.ma.array(np.zeros(res_pts*(len(depth_hycom)-1)), mask=True)
                int_depth_sec = int_depth_sec.reshape(res_pts,len(depth_hycom)-1)

                temp_reg = np.ma.array(np.zeros(days*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                temp_reg = temp_reg.reshape(days,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                temp_reg_mes = np.ma.array(np.zeros(months*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                temp_reg_mes = temp_reg_mes.reshape(months,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_temp = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_temp = s_temp.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                temp_sec_mes = np.ma.array(np.zeros(months*res_pts*len(depth_hycom)), mask=True)
                temp_sec_mes = temp_sec_mes.reshape(months,res_pts,len(depth_hycom))
                temp_sec = np.ma.array(np.zeros(res_pts*len(depth_hycom)), mask=True)
                temp_sec = temp_sec.reshape(res_pts,len(depth_hycom))

                sal_reg = np.ma.array(np.zeros(days*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                sal_reg = sal_reg.reshape(days,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                sal_reg_mes = np.ma.array(np.zeros(months*len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=True)
                sal_reg_mes = sal_reg_mes.reshape(months,len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_sal = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_sal = s_sal.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                sal_sec_mes = np.ma.array(np.zeros(months*res_pts*len(depth_hycom)), mask=True)
                sal_sec_mes = sal_sec_mes.reshape(months,res_pts,len(depth_hycom))
                sal_sec = np.ma.array(np.zeros(res_pts*len(depth_hycom)), mask=True)
                sal_sec = sal_sec.reshape(res_pts,len(depth_hycom))

                q = 0

            nc_file = nc.Dataset(arq_dia_hycom, 'r')
            u = nc_file.variables['u'][0,:,:,:] #(33, 564, 502)
            u = np.moveaxis(u,0,-1)
            v = nc_file.variables['v'][0,:,:,:]
            v = np.moveaxis(v,0,-1)
            nc_file.close()
            del nc_file
            nc_file = nc.Dataset(arq_dia_hycom2, 'r')
            int_depth = nc_file.variables['interface_depth'][0,:,:,:]
            int_depth = np.moveaxis(int_depth,0,-1)
            nc_file.close()
            del nc_file
            nc_file = nc.Dataset(arq_dia_hycom3, 'r')
            temp_mc = nc_file.variables['temperature'][0,:,:,:]
            temp_mc = np.moveaxis(temp_mc,0,-1)
            sal_mc = nc_file.variables['salinity'][0,:,:,:]
            sal_mc = np.moveaxis(sal_mc,0,-1)
            nc_file.close()
            del nc_file

            u = u[cond_lat,:,:][:,cond_lon,:]
            v = v[cond_lat,:,:][:,cond_lon,:]
            int_depth = int_depth[cond_lat,:,:][:,cond_lon,:]
            temp_mc = temp_mc[cond_lat,:,:][:,cond_lon,:]
            sal_mc = sal_mc[cond_lat,:,:][:,cond_lon,:]

            v_vel_reg[counter,:,:,:] = v
            u_vel_reg[counter,:,:,:] = u
            if 'int_depth' in locals():
                int_depth_reg[counter,:,:,:] = int_depth
            temp_reg[counter,:,:,:] = temp_mc
            sal_reg[counter,:,:,:] = sal_mc

            if (mes==current_data.month):
                s_v = s_v + np.squeeze(v_vel_reg[counter,:,:,:])
                s_u = s_u + np.squeeze(u_vel_reg[counter,:,:,:])
                if 'int_depth' in locals():
                    s_int_depth = s_int_depth + np.squeeze(int_depth_reg[counter,:,:,:])
                s_temp = s_temp + np.squeeze(temp_reg[counter,:,:,:])
                s_sal = s_sal + np.squeeze(sal_reg[counter,:,:,:])
                q = q + 1;
            else:
                counter_mes = counter_mes+1
                mes = current_data.month
                v_vel_reg_mes[counter_mes,:,:,:] = s_v/q
                u_vel_reg_mes[counter_mes,:,:,:] = s_u/q
                if 'int_depth' in locals():
                    int_depth_reg_mes[counter_mes,:,:,:] = s_int_depth/q
                temp_reg_mes[counter_mes,:,:,:] = s_temp/q
                sal_reg_mes[counter_mes,:,:,:] = s_sal/q

                s_v = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_v = s_v.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_u = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_u = s_u.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                if 'int_depth' in locals():
                    s_int_depth = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*(len(depth_hycom)-1)), mask=False)
                    s_int_depth = s_int_depth.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom)-1)
                s_temp = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_temp = s_temp.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                s_sal = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*len(depth_hycom)), mask=False)
                s_sal = s_sal.reshape(len(lat_hycom),len(lon_hycom),len(depth_hycom))
                q = 1;
                s_v = s_v + np.squeeze(v_vel_reg[counter,:,:,:]);
                s_u = s_u + np.squeeze(u_vel_reg[counter,:,:,:]);
                if 'int_depth' in locals():
                    s_int_depth = s_int_depth + np.squeeze(int_depth_reg[counter,:,:,:])
                s_temp = s_temp + np.squeeze(temp_reg[counter,:,:,:])
                s_sal = s_sal + np.squeeze(sal_reg[counter,:,:,:])

            current_data = current_data + datetime.timedelta(days=inc_tempo)

        if q>20:
            counter_mes = counter_mes+1
            v_vel_reg_mes[counter_mes,:,:,:] = s_v/q
            u_vel_reg_mes[counter_mes,:,:,:] = s_u/q
            if 'int_depth' in locals():
                int_depth_reg_mes[counter_mes,:,:,:] = s_int_depth/q
            temp_reg_mes[counter_mes,:,:,:] = s_temp/q
            sal_reg_mes[counter_mes,:,:,:] = s_sal/q

        for niv in range(0,len(depth_hycom),1):
            for mes in range(0,counter_mes+1,1):
                   
                for i in range(0,res_pts,1):
                    I = interpolate.interp2d(lon_hycom,lat_hycom,np.squeeze(v_vel_reg_mes[mes,:,:,niv].filled()),kind='linear')
                    v_vel_sec_mes[mes,i,niv] = I(lon_sec[i],lat_sec[i])
                    del I

                    I = interpolate.interp2d(lon_hycom,lat_hycom,np.squeeze(u_vel_reg_mes[mes,:,:,niv]).filled(),kind='linear')
                    u_vel_sec_mes[mes,i,niv] = I(lon_sec[i],lat_sec[i])
                    del I

                    I = interpolate.interp2d(lon_hycom,lat_hycom,np.squeeze(temp_reg_mes[mes,:,:,niv]).filled(),kind='linear')
                    temp_sec_mes[mes,i,niv] = I(lon_sec[i],lat_sec[i])
                    del I

                    I = interpolate.interp2d(lon_hycom,lat_hycom,np.squeeze(sal_reg_mes[mes,:,:,niv]).filled(),kind='linear')
                    sal_sec_mes[mes,i,niv] = I(lon_sec[i],lat_sec[i])
                    del I

                    if 's_int_depth' in locals():
                       if niv<len(layers):
                          I = interpolate.interp2d(lon_hycom,lat_hycom,np.squeeze(int_depth_reg_mes[mes,:,:,niv]).filled(),kind='linear')
                          int_depth_sec_mes[mes,i,niv] = I(lon_sec[i],lat_sec[i])
                          del I

                cond = abs(v_vel_sec_mes)>99999
                v_vel_sec_mes = np.ma.masked_where(cond,v_vel_sec_mes)
                u_vel_sec_mes = np.ma.masked_where(cond,u_vel_sec_mes)
                temp_sec_mes = np.ma.masked_where(cond,temp_sec_mes)
                sal_sec_mes = np.ma.masked_where(cond,sal_sec_mes)
                [u_vel_sec_mes[mes,:,niv], v_vel_sec_mes[mes,:,niv]] = np.dot(rot, [u_vel_sec_mes[mes,:,niv], v_vel_sec_mes[mes,:,niv]])

                cond = abs(int_depth_sec_mes)>99999
                int_depth_sec_mes = np.ma.masked_where(cond,int_depth_sec_mes)

        u_rot = np.mean(u_vel_sec_mes,axis=0)
        v_rot = np.mean(v_vel_sec_mes,axis=0)
        temp_sec = np.mean(temp_sec_mes,axis=0)
        sal_sec = np.mean(sal_sec_mes,axis=0)
        int_depth_sec = np.mean(int_depth_sec_mes,axis=0)

        ##### INTERPOLACAO WOA #########
        for mes in range(0,1,1):
            for i in range(0,res_pts,1):
                for niv in range(0,len(depth_woa_seas),1):
                    I = interpolate.interp2d(lon_woa,lat_woa,np.squeeze(temp_woa[mes,:,:,niv]),kind='linear')
                    temp[mes,i,niv] = I(lat_sec[i],lon_sec[i])

                    I = interpolate.interp2d(lon_woa,lat_woa,np.squeeze(sal_woa[mes,:,:,niv]),kind='linear')
                    sal[mes,i,niv] = I(lat_sec[i],lon_sec[i])

                temp_interp = np.squeeze(temp[mes,i,:])
                cond = temp_interp.mask==False
                I = interpolate.interp1d(depth_woa_seas[cond],temp_interp[cond].data)
                cond = depth_hycom>=max(depth_woa_seas[cond])
                temp_sec_woa_mes[mes,i,:] = I(depth_hycom[cond].data)

                sal_interp = np.squeeze(sal[mes,i,:])
                cond = sal_interp.mask==False
                I = interpolate.interp1d(depth_woa_seas[cond],sal_interp[cond].data)
                cond = depth_hycom>=max(depth_woa_seas[cond])
                sal_sec_woa_mes[mes,i,:] = I(depth_hycom[cond].data)
        temp_sec_woa[:,:] = np.mean(temp_sec_woa_mes,axis=0)
        sal_sec_woa[:,:] = np.mean(sal_sec_woa_mes,axis=0)

        filename = "VEL_ROT_"+rodada[rod]+"_"+data_inicial.strftime("%Y%m%d")+"-"+data_final.strftime("%Y%m%d")+"_"+str(abs(p1[0]))+s_lat_p1+ \
                   str(abs(p1[1]))+s_lon_p1+'-'+str(abs(p2[0]))+s_lat_p2+'_'+str(abs(p2[1]))+s_lon_p2
        if not (os.path.isdir(output_dir+"/"+rodada[rod])):
           os.system("mkdir -p "+output_dir+"/"+rodada[rod])
        scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'u_rot_'+rodada[rod]: u_rot, 'v_rot_'+rodada[rod]: v_rot, \
                         'int_depth_'+rodada[rod]: int_depth_sec, 'temp_'+rodada[rod]: temp_sec, 'sal_'+rodada[rod]: sal_sec, \
                         'temp_sec_woa_'+rodada[rod]: temp_sec_woa, 'sal_sec_woa_'+rodada[rod]: sal_sec_woa, 'lat_sec': lat_sec, \
                         'lon_sec': lon_sec, 'depth_hycom': depth_hycom})

        ano = data_inicial.year
        mes = data_inicial.month
        for mm in range(0,months,1):
            u_rot = np.squeeze(u_vel_sec_mes[mm,:,:])
            v_rot = np.squeeze(v_vel_sec_mes[mm,:,:])
            temp_sec = np.squeeze(temp_sec_mes[mm,:,:])
            sal_sec = np.squeeze(sal_sec_mes[mm,:,:])
            int_depth_sec = np.squeeze(int_depth_sec_mes[mm,:,:])

            temp_sec_woa = np.squeeze(temp_sec_woa_mes[mm,:,:])
            sal_sec_woa = np.squeeze(sal_sec_woa_mes[mm,:,:])

            filename = "VEL_ROT_"+rodada[rod]+"_"+str(ano)+str(mes).zfill(2)+"01-"+str(ano)+str(mes).zfill(2)+"30_"+str(abs(p1[0]))+ \
                       s_lat_p1+str(abs(p1[1]))+s_lon_p1+'-'+str(abs(p2[0]))+s_lat_p2+'_'+str(abs(p2[1]))+s_lon_p2
            scipy.io.savemat(output_dir+"/"+rodada[rod]+"/"+filename+".mat", mdict={'u_rot_'+rodada[rod]: u_rot, 'v_rot_'+rodada[rod]: v_rot, \
                             'int_depth_'+rodada[rod]: int_depth_sec, 'temp_'+rodada[rod]: temp_sec, 'sal_'+rodada[rod]: sal_sec, \
                             'temp_sec_woa_'+rodada[rod]: temp_sec_woa, 'sal_sec_woa_'+rodada[rod]: sal_sec_woa, 'lat_sec': lat_sec, \
                             'lon_sec': lon_sec, 'depth_hycom': depth_hycom})
            mes = mes+1
            if mes>12:
               mes = 1
               ano = ano+1
        #del v_vel_reg_mes, u_vel_reg_mes, temp_reg_mes, sal_reg_mes, int_depth_reg_mes, s_int_depth 


        print()
        print("EXPT: "+rodada[rod]+" LAT:"+str(lats[ll])+" ---- OK")
        print()


sys.exit()

