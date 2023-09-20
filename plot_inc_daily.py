import remo
import sys
import numpy as np
import datetime
import time
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
jet = matplotlib.cm.get_cmap('jet')

restar_type = True #True==restart / False==archv
ssh_archv = True

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))
inc_assim = float(content_list[6].rstrip('\n'))

dir_hycom = '/home/filipe.costa/previsao/hycom_2_2_18/proc'
fig_output_dir = '/home/filipe.costa/resultados/figs'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       lat_step_plot = 4
       lon_step_plot = 5
       resolution_plot = 150
       if restar_type:
          record_u = 0
          record_v = 64
          record_dp = 128
          record_temp = 192
          record_saln = 256
       dir_ATL = dir_hycom+'/ATLj0.04/'

    elif '008i' in rodada[rod]:
       lat_step_plot = 5
       lon_step_plot = 5
       resolution_plot = 150

       min_temp = 15
       max_temp = 30
       min_temp_diff = -4
       max_temp_diff = 4

       min_ssh = -0.4
       max_ssh = 1.2
       min_ssh_diff = -0.4
       max_ssh_diff = 0.4

       IDM = 628
       JDM = 780
       if restar_type:
          record_u = 0
          record_v = 64
          record_dp = 128
          record_temp = 192
          record_saln = 256
       dir_ATL = dir_hycom+'/ATLi0.08/'

    elif '008d' in rodada[rod]:
       lat_step_plot = 15
       lon_step_plot = 15
       resolution_plot = 150

       min_temp = -2
       max_temp = 30
       min_temp_diff = -4
       max_temp_diff = 4

       min_ssh = -0.4
       max_ssh = 1.2
       min_ssh_diff = -0.4
       max_ssh_diff = 0.4

       IDM = 1717
       JDM = 2345
       if restar_type:
          record_u = 0
          record_v = 64
          record_dp = 128
          record_temp = 192
          record_saln = 256
       dir_ATL = dir_hycom+'/ATLd0.08/'
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

        arq_back = dir_expt+'/data/restart_files/restart_'+str((current_data).timetuple().tm_year).zfill(4)+'d' \
                        +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'

        arq_ana = dir_expt+'/data/restart_assim/restart_'+str((current_data).timetuple().tm_year).zfill(4)+'d' \
                        +str((current_data).timetuple().tm_yday).zfill(3)+'h00.a'

        if ssh_archv:
           arq_back_ssh = dir_expt+'/output/ab/'+(current_data - datetime.timedelta(days=1)).strftime("%Y%m%d")+ \
                          '/archv.'+current_data.strftime("%Y")+'_' \
                          +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
           if not (os.path.isfile(arq_back_ssh)):
              arq_back_ssh = dir_expt+'/output/ab/'+(current_data - datetime.timedelta(days=inc_assim)).strftime("%Y%m%d")+ \
                             '/archv.'+current_data.strftime("%Y")+'_' \
                             +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
           arq_ana_ssh = dir_expt+'/output/ab/'+current_data.strftime("%Y%m%d")+'/archv.'+current_data.strftime("%Y")+'_' \
                         +str((current_data).timetuple().tm_yday).zfill(3)+'_00.a'
        else:
           arq_ana_ssh = dir_expt+'/assim/slasst_gridded_map_error_or_argoZ_stats/Check/'+current_data.strftime("%Y%m%d")+'00/'+ \
                         'increment_adt.ascii'
        if not (os.path.isfile(arq_back)):
           print('HYCOM BACK FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_back)
           current_data = current_data + datetime.timedelta(days=inc_assim)
           continue
        if not (os.path.isfile(arq_ana)):
           print('HYCOM ANA FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
           print(arq_ana)
           current_data = current_data + datetime.timedelta(days=inc_assim)
           continue

        f_back = open(arq_back,'rb')
        f_ana  = open(arq_ana,'rb')

        f_back.seek(record_temp*4*(IJDM+npad))
        field_back = np.fromfile(f_back,dtype='>f4',count=IDM*JDM)
        field_back = np.reshape(field_back,(JDM,IDM))
        cond = abs(field_back) > 999999
        field_back = np.ma.masked_where(cond,field_back)

        f_ana.seek(record_temp*4*(IJDM+npad))
        field_ana = np.fromfile(f_ana,dtype='>f4',count=IDM*JDM)
        field_ana = np.reshape(field_ana,(JDM,IDM))
        cond = abs(field_ana) > 999999
        field_ana = np.ma.masked_where(cond,field_ana)

        out_dir = fig_output_dir+"/"+rodada[rod]+"/TOTAL_INC/MAPS/DAILY"
        if not (os.path.isdir(out_dir)):
           os.system("mkdir -p "+out_dir)

        remo.plot_map(lat_hycom, lon_hycom, field_back, jet,shad='gouraud',cbmin=min_temp,cbmax=max_temp,um='\N{DEGREE SIGN}C', \
                      lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' BACK', \
                      filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_BACK.png')

        remo.plot_map(lat_hycom, lon_hycom, field_ana, jet,shad='gouraud',cbmin=min_temp,cbmax=max_temp,um='\N{DEGREE SIGN}C', \
                      lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' ANA', \
                      filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_ANA.png')

        remo.plot_map(lat_hycom, lon_hycom, field_ana-field_back, jet,shad='gouraud',cbmin=min_temp_diff,cbmax=max_temp_diff,um='\N{DEGREE SIGN}C', \
                      lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' INC', \
                      filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_INC.png')

        f_ana.close()
        f_back.close()
        f.close()

        if ssh_archv:
           if not (os.path.isfile(arq_back_ssh)):
              print('HYCOM BACK FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
              print(arq_back_ssh)
              current_data = current_data + datetime.timedelta(days=inc_assim)
              continue
           if not (os.path.isfile(arq_ana_ssh)):
              print('HYCOM ANA FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
              print(arq_ana_ssh)
              current_data = current_data + datetime.timedelta(days=inc_assim)
              continue

           f_back = open(arq_back_ssh,'rb')
           f_ana  = open(arq_ana_ssh,'rb')

           f_back.seek(1*4*(IJDM+npad))
           field_back = np.fromfile(f_back,dtype='>f4',count=IDM*JDM)
           field_back = np.reshape(field_back,(JDM,IDM))/9.806
           cond = abs(field_back) > 999999
           field_back = np.ma.masked_where(cond,field_back)

           f_ana.seek(1*4*(IJDM+npad))
           field_ana = np.fromfile(f_ana,dtype='>f4',count=IDM*JDM)
           field_ana = np.reshape(field_ana,(JDM,IDM))/9.806
           cond = abs(field_ana) > 999999
           field_ana = np.ma.masked_where(cond,field_ana)

           out_dir = fig_output_dir+"/"+rodada[rod]+"/TOTAL_INC/MAPS/DAILY"
           if not (os.path.isdir(out_dir)):
              os.system("mkdir -p "+out_dir)

           remo.plot_map(lat_hycom, lon_hycom, field_back, jet,shad='gouraud',cbmin=min_ssh,cbmax=max_ssh,um='m', \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' BACK', \
                         filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_BACK.png')

           remo.plot_map(lat_hycom, lon_hycom, field_ana, jet,shad='gouraud',cbmin=min_ssh,cbmax=max_ssh,um='m', \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' ANA', \
                         filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_ANA.png')

           remo.plot_map(lat_hycom, lon_hycom, field_ana-field_back, jet,shad='gouraud',cbmin=min_ssh_diff,cbmax=max_ssh_diff,um='m', \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot,subtitle=rodada[rod]+' INC', \
                         filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_INC.png')

           f_ana.close()
           f_back.close()

        else:
           ssh = pd.read_table(arq_ana_ssh, header=None, sep='\s+')
           cond = abs(ssh) > 999999
           ssh = np.ma.masked_where(cond,ssh)

           out_dir = fig_output_dir+"/"+rodada[rod]+"/TOTAL_INC/MAPS/DAILY"
           if not (os.path.isdir(out_dir)):
              os.system("mkdir -p "+out_dir)

           remo.plot_map(lat_hycom, lon_hycom, ssh, jet,shad='gouraud',cbmin=min_ssh_diff,cbmax=max_ssh_diff, \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                         filename=out_dir+'/INC_ADT_HYCOM_'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_ANA.png')


        current_data = current_data + datetime.timedelta(days=inc_assim)
