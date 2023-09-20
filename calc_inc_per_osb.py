import remo
import sys
import numpy as np
import datetime
import time
import os
import copy
import matplotlib
import matplotlib.pyplot as plt
cmap = copy.copy(matplotlib.cm.get_cmap('jet'))
#jet = matplotlib.cm.get_cmap('jet')

#obs_types = ['adt_track','sal_roms','tem_roms','sst_roms','swot_err_big','swot_err_med','swot_err_sma']
obs_types = ['adt_track','sst','tem_argoZ','sal_argoZ','tem_mrbZ','sal_mrbZ','tem_gldZ','sal_gldZ']
model_vars = ['ssh','sst']
#model_vars = ['ssh']

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
save_daily_figs = True
#save_daily_figs = False
#save_mean_figs = True
save_mean_figs = False

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])

for rod in range(0,len(rodada),1):
    current_data = data_inicial
    if '004j' in rodada[rod]:
       IDM = 502
       JDM = 564
       dir_ATL = dir_hycom+'/ATLj0.04/'
    elif '008i' in rodada[rod]:
       IDM = 628
       JDM = 780
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
        counter = counter +1

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

           mean_inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*(len(obs_types)+1)*len(model_vars)),mask=True)
           mean_inc = mean_inc.reshape(len(lat_hycom),len(lon_hycom),len(obs_types)+1,len(model_vars))
           cond_mean_inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*(len(obs_types)+1)*len(model_vars)),mask=True)
           cond_mean_inc = cond_mean_inc.reshape(len(lat_hycom),len(lon_hycom),len(obs_types)+1,len(model_vars))

           inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*int((data_final - data_inicial).days/inc_assim+1)* \
                 (len(obs_types)+1)*len(model_vars)),mask=False)
           inc = inc.reshape(len(lat_hycom),len(lon_hycom),int((data_final - data_inicial).days/inc_assim+1),len(obs_types)+1,len(model_vars))

           cond_inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)*int((data_final - data_inicial).days/inc_assim+1)* \
                      (len(obs_types)+1)),mask=True)
           cond_inc = cond_inc.reshape(len(lat_hycom),len(lon_hycom),int((data_final - data_inicial).days/inc_assim+1),len(obs_types)+1)



        for var in range(0,len(model_vars),1):
            #record_u = 0
            #record_v = 32
            #record_dp = 64
            #record_saln = 128
            if 'ssh' in model_vars[var]:
               record = 160
            elif 'sst' in model_vars[var]:
               record = 98
            else:
               print(model_vars[var]+' not found in record options.')

            tot_inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=False)
            tot_inc =  tot_inc.reshape(len(lat_hycom),len(lon_hycom))
            cond_tot_inc = np.ma.array(np.zeros(len(lat_hycom)*len(lon_hycom)), mask=True)
            cond_tot_inc =  cond_tot_inc.reshape(len(lat_hycom),len(lon_hycom))

            for obs in range(0,len(obs_types),1):  #FAZENDO INC TOTAL
                arq = dir_expt+'/assim/slasst_gridded_map_error_or_argoZ_stats/Check/'+current_data.strftime("%Y%m%d")+'00/inc_'+ \
                      obs_types[obs]+'.a'

                if not (os.path.isfile(arq)):
                   print('INC FOR OBS: '+obs_types[obs]+' DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
                   print(arq)
                   continue

                f = open(arq,'rb')

                f.seek(record*4*(IJDM+npad))
                field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
                field = np.reshape(field,(JDM,IDM))

                cond = abs(field) < 999

                tot_inc[cond] = tot_inc[cond] + abs(field[cond])
                cond_tot_inc.mask[cond] = False

                f.close()
                del field, cond


            for obs in range(0,len(obs_types),1):  #FAZENDO INC/TOT_INC

                arq = dir_expt+'/assim/slasst_gridded_map_error_or_argoZ_stats/Check/'+current_data.strftime("%Y%m%d")+'00/inc_'+ \
                      obs_types[obs]+'.a'

                if not (os.path.isfile(arq)):
                   print('INC FOR OBS: '+obs_types[obs]+' DAY: '+current_data.strftime("%Y%m%d")+' DOES NOT EXIST')
                   print(arq)
                   continue

                f = open(arq,'rb')

                f.seek(record*4*(IJDM+npad))
                field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
                field = np.reshape(field,(JDM,IDM))

                cond = abs(field) < 999
                inc[cond,counter,obs,var] = abs(field[cond])/tot_inc[cond]*100
                cond_inc.mask[cond,counter,obs] = False
                cond_mean_inc.mask[cond,obs,var] = False
                if 'swot' in obs_types[obs]:
                   inc[cond,counter,len(obs_types),var] = inc[cond,counter,len(obs_types),var] + inc[cond,counter,obs,var]
                   cond_inc.mask[cond,counter,len(obs_types)] = False
                   cond_mean_inc.mask[cond,len(obs_types),var] = False

                inc[:,:,counter,obs,var] = np.ma.masked_where(np.ma.getmask(cond_inc[:,:,counter,obs]),inc[:,:,counter,obs,var])

                f.close()
                del field, cond

                if save_daily_figs:
                   out_dir = fig_output_dir+"/"+rodada[rod]+"/INC_PER_OBS/MAPS/DAILY"
                   if not (os.path.isdir(out_dir)):
                      os.system("mkdir -p "+out_dir)

                   remo.plot_map(lat_hycom, lon_hycom, inc[:,:,counter,obs,var], cmap,shad='gouraud',cbmin=0,cbmax=100, \
                                 lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                                 title='CONTRIBUTION OF '+str(obs_types[obs]).upper()+'\n ON '+str(model_vars[var]).upper()+' ('+ \
                                 current_data.strftime("%d/%m/%Y")+') '+rodada[rod],um='%', \
                                 filename=out_dir+'/'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_'+str(model_vars[var]).upper()+'_'+ \
                                 str(obs_types[obs]).upper()+'.png')

            inc[:,:,counter,len(obs_types),var] = np.ma.masked_where(np.ma.getmask(cond_inc[:,:,counter,len(obs_types)]), \
                                                  inc[:,:,counter,len(obs_types),var])

            if save_daily_figs and 'swot' in obs_types[:]:
               out_dir = fig_output_dir+"/"+rodada[rod]+"/INC_PER_OBS/MAPS/DAILY"
               remo.plot_map(lat_hycom, lon_hycom, inc[:,:,counter,len(obs_types),var], cmap,shad='gouraud',cbmin=0,cbmax=100, \
                             lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                             title='CONTRIBUTION OF SWOT \n ON '+str(model_vars[var]).upper()+' ('+ \
                             current_data.strftime("%d/%m/%Y")+') '+rodada[rod],um='%', \
                             filename=out_dir+'/'+rodada[rod]+'_'+current_data.strftime("%Y%m%d")+'_'+str(model_vars[var]).upper()+'_SWOT.png')

        current_data = current_data + datetime.timedelta(days=inc_assim)

    if save_mean_figs:
       mean_inc = np.nanmean(inc,2)
       mean_inc = np.ma.masked_where(np.ma.getmask(cond_mean_inc),mean_inc)
       for var in range(0,len(model_vars),1):
           for obs in range(0,len(obs_types),1):

               out_dir = fig_output_dir+"/"+rodada[rod]+"/INC_PER_OBS/MAPS"
               if not (os.path.isdir(out_dir)):
                  os.system("mkdir -p "+out_dir)

               remo.plot_map(lat_hycom, lon_hycom, mean_inc[:,:,obs,var], cmap,shad='gouraud',cbmin=0,cbmax=100, \
                             lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                             title='CONTRIBUTION OF '+str(obs_types[obs]).upper()+'\n ON '+str(model_vars[var]).upper()+' ('+ \
                             data_inicial.strftime("%d/%m/%Y")+' - '+data_final.strftime("%d/%m/%Y")+') '+rodada[rod],um='%', \
                             filename=out_dir+'/'+rodada[rod]+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+'_'+ \
                             str(model_vars[var]).upper()+'_'+str(obs_types[obs]).upper()+'.png')

           if 'swot' in obs_types or 'swot_err_sma' in obs_types:
              remo.plot_map(lat_hycom, lon_hycom, mean_inc[:,:,len(obs_types),var], cmap,shad='gouraud',cbmin=0,cbmax=100, \
                            lonstep=lon_step_plot,latstep=lat_step_plot, resolution=resolution_plot, \
                            title='CONTRIBUTION OF SWOT \n ON '+str(model_vars[var]).upper()+' ('+ \
                            data_inicial.strftime("%d/%m/%Y")+' - '+data_final.strftime("%d/%m/%Y")+') '+rodada[rod],um='%', \
                            filename=out_dir+'/'+rodada[rod]+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+'_'+ \
                            str(model_vars[var]).upper()+'_SWOT.png')




