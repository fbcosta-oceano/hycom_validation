import remo
import sys
import numpy as np
import datetime
import time
import os
import matplotlib
import matplotlib.pyplot as plt
jet = matplotlib.cm.get_cmap('jet')

obs_types = ['adt_track','sal_roms','tem_roms','sst_roms','swot_err_big','swot_err_med','swot_err_sma']
#obs_types = ['swot_err_big','swot_err_med','swot_err_sma']
obs_types = ['adt_track','sst','tem_argoZ','sal_argoZ','tem_mrbZ','sal_mrbZ','tem_gldZ','sal_gldZ']

lat_step_plot = 4
lon_step_plot = 5
resolution_plot = 150

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
       record_u = 0
       record_v = 32
       record_dp = 64
       record_temp = 98
       record_saln = 128
       record_ssh = 160
       dir_ATL = dir_hycom+'/ATLj0.04/'
    elif '008i' in rodada[rod]:
       IDM = 628
       JDM = 780
       record_u = 0
       record_v = 32
       record_dp = 64
       record_temp = 98
       record_saln = 128
       record_ssh = 160
       dir_ATL = dir_hycom+'/ATLi0.08/'
    else:
       print('CASE '+rodada[rod]+' NOT DEFINED')
       continue

    dir_expt = dir_ATL+expt[rod]
    IJDM=IDM*JDM
    npad=4096-(IJDM%4096)
    counter = -1

    out_dir = fig_output_dir+'/'+rodada[rod]+'/INC_PER_OBS'
    if not (os.path.isdir(out_dir)):
       os.system("mkdir -p "+out_dir)
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

        for obs in range(0,len(obs_types),1):
            arq = dir_expt+'/assim/slasst_gridded_map_error_or_argoZ_stats/Check/'+current_data.strftime("%Y%m%d")+'00/inc_'+obs_types[obs]+'.a'
            print(arq)

            if not (os.path.isfile(arq)):
               print('HYCOM FILE FOR DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
               print(arq)
               continue

            f = open(arq,'rb')

            f.seek(record_temp*4*(IJDM+npad))
            field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
            field = np.reshape(field,(JDM,IDM))

            if obs==0:
               inc_swot_sst = np.ma.array(np.zeros(field.shape), mask=False)
               inc_swot_ssh = np.ma.array(np.zeros(field.shape), mask=False)
               cond_swot = np.ma.array(np.zeros(field.shape), mask=True)

               inc_argo_sst = np.ma.array(np.zeros(field.shape), mask=False)
               inc_argo_ssh = np.ma.array(np.zeros(field.shape), mask=False)
               cond_argo = np.ma.array(np.zeros(field.shape), mask=True)

            cond = abs(field) < 999
            if 'swot' in obs_types[obs]:
               inc_swot_sst[cond] = inc_swot_sst[cond] + field[cond]
               cond_swot.mask[cond] = False
            if 'tem' in obs_types[obs] or 'sal' in obs_types[obs]:
               inc_argo_sst[cond] = inc_argo_sst[cond] + field[cond]
               cond_argo.mask[cond] = False

            cond = abs(field) > 999
            field = np.ma.masked_where(cond,field)
            #print(field)

            remo.plot_map(lat_hycom, lon_hycom, field, jet,shad='gouraud',cbmin=-4,cbmax=4, \
                          lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                          title='INC ON SST FROM '+str(obs_types[obs]).upper(), \
                          filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_'+str(obs_types[obs]).upper()+'.png')


            del cond,field

            f.seek(record_ssh*4*(IJDM+npad))
            field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
            field = np.reshape(field,(JDM,IDM))
            field = field/9.806
            cond = abs(field) < 999
            if 'swot' in obs_types[obs]:
               inc_swot_ssh[cond] = inc_swot_ssh[cond] + field[cond]
            if 'tem' in obs_types[obs] or 'sal' in obs_types[obs]:
               inc_argo_ssh[cond] = inc_argo_ssh[cond] + field[cond]

            cond = abs(field) > 999
            field = np.ma.masked_where(cond,field)
            #print(field)

            remo.plot_map(lat_hycom, lon_hycom, field, jet,shad='gouraud',cbmin=-0.2,cbmax=0.2, \
                          lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                          title='INC ON SSH FROM '+str(obs_types[obs]).upper(), \
                          filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_'+str(obs_types[obs]).upper()+'.png')

            f.close()
            del field, cond

        if 'swot' in str(obs_types):
           inc_swot_sst = np.ma.masked_where(cond_swot,inc_swot_sst)
           remo.plot_map(lat_hycom, lon_hycom, inc_swot_sst, jet,shad='gouraud',cbmin=-4,cbmax=4, \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                         title='INC ON SST FROM SWOT', \
                         filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_SWOT.png')

           inc_swot_ssh = np.ma.masked_where(cond_swot,inc_swot_ssh)
           remo.plot_map(lat_hycom, lon_hycom, inc_swot_ssh, jet,shad='gouraud',cbmin=-0.2,cbmax=0.2, \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                         title='INC ON SSH FROM SWOT', \
                         filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_SWOT.png')

        if 'tem' in str(obs_types) or 'sal' in str(obs_types):
           inc_argo_sst = np.ma.masked_where(cond_argo,inc_argo_sst)
           remo.plot_map(lat_hycom, lon_hycom, inc_argo_sst, jet,shad='gouraud',cbmin=-4,cbmax=4, \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                         title='INC ON SST FROM ARGO', \
                         filename=out_dir+'/SST_HYCOM_'+rodada[rod]+'_ARGO.png')

           inc_argo_ssh = np.ma.masked_where(cond_argo,inc_argo_ssh)
           remo.plot_map(lat_hycom, lon_hycom, inc_argo_ssh, jet,shad='gouraud',cbmin=-0.2,cbmax=0.2, \
                         lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                         title='INC ON SSH FROM ARGO', \
                         filename=out_dir+'/SSH_HYCOM_'+rodada[rod]+'_ARGO.png')

        current_data = current_data + datetime.timedelta(days=inc_assim)
