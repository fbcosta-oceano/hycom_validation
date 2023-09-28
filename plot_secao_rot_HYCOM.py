from shiftedColorMap import shiftedColorMap
import remo
import datetime
import time
import os
import sys
import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
from calendar import monthrange
warnings.filterwarnings("ignore")

months = [1,2,3,4,5,6,7,8,9,10,11,12]
#months = [1,2,3,4,5,6,7,8]

max_prof_sup = 500
max_prof_inf = 1500

cfg_file=sys.argv[1]
opt=open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip('\n').split(','))
rodada = list(content_list[1].rstrip('\n').split(','))
legenda = list(content_list[2].rstrip('\n').split(','))
data_inicial = str(content_list[3].rstrip('\n'))
data_final = str(content_list[4].rstrip('\n'))
inc_tempo = float(content_list[5].rstrip('\n'))

input_dir = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/resultados'
fig_output_dir = '/mnt/nfs/dpns33/data1/home_dpns31/fbcosta/resultados/figs'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
total_months = (data_final.year - data_inicial.year)*12 + data_final.month - data_inicial.month +1

font_size = 8
resolution_plot = 150
#under='#000033'
under='#0000c2'
over='#330000'
bad='#000000'
cmap=matplotlib.cm.get_cmap('bwr')
cmap.set_under(under)
cmap.set_over(over)
cmap.set_bad(color='k', alpha=1)

#lats = np.array([-22])
lats = np.array([-22, -23, -24, -25, -26, -27, -28, -29, -30, -32])

for ll in range(0,len(lats),1):
    if lats[ll]==-9:
       p1 = [-08.80, -35.20] #[lat lon]
       p2 = [-09.50, -32.50] #[lat lon]

    elif lats[ll]==-10:
       p1 = [-09.60, -35.50] #[lat lon]
       p2 = [-11.10, -33.00] #[lat lon]    

    elif lats[ll]==-22:
       p1 = [-22.00, -41.10] #[lat lon]
       p2 = [-22.00, -39.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625
    elif lats[ll]==-23:
       p1 = [-22.25, -42.00] #[lat lon]
       p2 = [-23.70, -39.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-24:
       p1 = [-22.80, -42.25] #[lat lon]
       p2 = [-25.50, -42.25] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-25:
       p1 = [-24.00, -46.00] #[lat lon]
       p2 = [-26.40, -43.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-26:
       p1 = [-24.60, -47.50] #[lat lon]
       p2 = [-26.70, -44.50] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-27:
       p1 = [-26.00, -48.50] #[lat lon]
       p2 = [-28.10, -44.50] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-28:
       p1 = [-27.25, -48.70] #[lat lon]
       p2 = [-28.55, -46.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-29:
       p1 = [-28.35, -49.00] #[lat lon]
       p2 = [-29.50, -46.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-30:
       p1 = [-29.30, -50.00] #[lat lon]
       p2 = [-30.70, -46.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

    elif lats[ll]==-32:
       p1 = [-32.00, -51.00] #[lat lon]
       p2 = [-32.00, -48.00] #[lat lon]
       prof_plot_sup = 400
       prof_plot_inf = 1500
       min_mean = -0.8
       max_mean = 0.5
       inc_mean = 0.05
       midpoint_colorbar=0.625

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

    for rod in range(0,len(rodada),1):
        current_data = data_inicial
        #y = data_inicial.year
        counter = -1

        while current_data<=data_final:
            while current_data.month != months[0]:
               days = monthrange(current_data.year,current_data.month)[1]
               current_data = current_data + datetime.timedelta(days=days)
            for m in list(months):
                days = monthrange(current_data.year,current_data.month)[1]
                filename = input_dir+"/"+rodada[rod]+"/VEL_ROT_"+rodada[rod]+"_"+str(current_data.year).zfill(4)+ \
                           str(m).zfill(2)+"01-"+str(current_data.year).zfill(4)+str(m).zfill(2)+"30_"+str(abs(p1[0]))+s_lat_p1+ \
                           str(abs(p1[1]))+s_lon_p1+'-'+str(abs(p2[0]))+s_lat_p2+'_'+str(abs(p2[1]))+s_lon_p2+".mat"
                if not (os.path.isfile(filename)):
                   print(filename, ' NAO ENCONTRADO')
                   exit()
                counter = counter+1
                mat_contents = scipy.io.loadmat(filename)

                if not 'lat_sec' in locals():
                   lat_sec = np.squeeze(mat_contents['lat_sec'])
                   lon_sec = np.squeeze(mat_contents['lon_sec'])
                   depth_hycom = np.squeeze(mat_contents['depth_hycom'])
                if not 'v_rot_all' in locals():
                   u_rot_all = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                   u_rot_all = u_rot_all.reshape(total_months,len(lat_sec),len(depth_hycom))
                   v_rot_all = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                   v_rot_all = v_rot_all.reshape(total_months,len(lat_sec),len(depth_hycom))
                   int_depth_all = np.ma.array(np.zeros(len(lat_sec)*(len(depth_hycom)-1)*total_months), mask=True)
                   int_depth_all = int_depth_all.reshape(total_months,len(lat_sec),len(depth_hycom)-1)
                   temp_all = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                   temp_all = temp_all.reshape(total_months,len(lat_sec),len(depth_hycom))
                   sal_all = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                   sal_all = sal_all.reshape(total_months,len(lat_sec),len(depth_hycom))

                   if rod==0:
                      temp_all_woa = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                      temp_all_woa = temp_all_woa.reshape(total_months,len(lat_sec),len(depth_hycom))
                      sal_all_woa = np.ma.array(np.zeros(len(lat_sec)*len(depth_hycom)*total_months), mask=True)
                      sal_all_woa = sal_all_woa.reshape(total_months,len(lat_sec),len(depth_hycom))

                u_rot_all[counter,:,:] = mat_contents['u_rot_'+rodada[rod]]
                v_rot_all[counter,:,:] = mat_contents['v_rot_'+rodada[rod]]
                int_depth_all[counter,:,:] = mat_contents['int_depth_'+rodada[rod]]
                temp_all[counter,:,:] = mat_contents['temp_'+rodada[rod]]
                sal_all[counter,:,:] = mat_contents['sal_'+rodada[rod]]

                if rod==0:
                   temp_all_woa[counter,:,:] = mat_contents['temp_sec_woa_'+rodada[rod]]
                   sal_all_woa[counter,:,:] = mat_contents['sal_sec_woa_'+rodada[rod]]


                print("EXPT: "+rodada[rod]+" LAT: "+str(lats[ll])+" MES: "+str(m).zfill(2)+" ANO: "+str(current_data.year)+" ---- OK")
                current_data = current_data + datetime.timedelta(days=days)

        cond = abs(u_rot_all)>99999
        u_rot_all = np.ma.masked_where(cond,u_rot_all)
        v_rot_all = np.ma.masked_where(cond,v_rot_all)
        temp_all = np.ma.masked_where(cond,temp_all)
        sal_all = np.ma.masked_where(cond,sal_all)
        cond = abs(int_depth_all)>99999
        int_depth_all = np.ma.masked_where(cond,int_depth_all)
        locals()['u_rot_'+rodada[rod]] = np.squeeze(np.mean(u_rot_all,axis=0))
        locals()['v_rot_'+rodada[rod]] = np.squeeze(np.mean(v_rot_all,axis=0))
        locals()['int_depth_'+rodada[rod]] = np.squeeze(np.mean(int_depth_all,axis=0))
        locals()['temp_'+rodada[rod]] = np.squeeze(np.mean(temp_all,axis=0))
        locals()['sal_'+rodada[rod]] = np.squeeze(np.mean(sal_all,axis=0))
        if rod==0:
           cond = abs(temp_all_woa)>99999
           temp_all_woa = np.ma.masked_where(cond,temp_all_woa)
           sal_all_woa = np.ma.masked_where(cond,sal_all_woa)
           temp_woa = np.squeeze(np.mean(temp_all_woa,axis=0)) 
           sal_woa = np.squeeze(np.mean(sal_all_woa,axis=0)) 

    all_rod_save = rodada[0]
    print(locals()['int_depth_'+rodada[0]].shape)
    for rod in range(1,len(rodada),1):
        all_rod_save = all_rod_save+"_"+rodada[rod]
    out_dir = fig_output_dir+"/"+all_rod_save+"/V_VEL_ROT/"+str(abs(lats[ll]))+s_lat_p2
    if not (os.path.isdir(out_dir)):
       os.system("mkdir -p "+out_dir)

    new_cmap = shiftedColorMap(cmap, midpoint=midpoint_colorbar, name='shifted')

    fig = plt.figure(figsize=(6.4, 4.8), constrained_layout=True)
    outer = gridspec.GridSpec(len(rodada), 1)
    for rod in range(0,len(rodada),1):

        data = np.transpose(locals()['v_rot_'+rodada[rod]])
        ind_x = np.ma.notmasked_edges(data[0,:], axis=None)[0]+5

        gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[rod], hspace = 0.05)
        ax = plt.subplot(gs1[0])

        int_depth = np.transpose(locals()['int_depth_'+rodada[rod]])
        int_depth = np.ma.masked_where(data[0:int_depth.shape[0],:].mask,int_depth)
        plt.plot(lon_sec[ind_x:-1],int_depth[9,ind_x:-1], color='k', linestyle='--', linewidth=.5, alpha=.5)
        plt.plot(lon_sec[ind_x:-1],int_depth[14,ind_x:-1], color='k', linestyle='--', linewidth=.5, alpha=.5)

        ind_y = depth_hycom<=prof_plot_sup
        plt.gca().patch.set_color('.9')
        plt.gca().patch.set_alpha(0.9)

        if lats[ll]==-24:
           im = ax.contourf(lat_sec[ind_x:-1], depth_hycom[ind_y], data[ind_y,ind_x:-1],np.arange(min_mean,max_mean+inc_mean,inc_mean), \
                   cmap=new_cmap,vmin=min_mean,vmax=max_mean,linestyles='dashed',extend='both')
        else:
           im = ax.contourf(lon_sec[ind_x:-1], depth_hycom[ind_y], data[ind_y,ind_x:-1],np.arange(min_mean,max_mean+inc_mean,inc_mean), \
                   cmap=new_cmap,vmin=min_mean,vmax=max_mean,linestyles='dashed',extend='both')

        #levels = np.ma.array(np.zeros(len(data[0,ind_x:-1])), mask=True, dtype=np.int16)
        #for i in range(ind_x,len(data[0,ind_x:-1]),1):
        #    levels[i] = depth_hycom[np.ma.notmasked_edges(data[:,i], axis=None)[-1]]
        #plt.plot(lon_sec[ind_x:-1],levels, color='k', linestyle='-', linewidth=.5)
        
        ax.set_ylim(0, prof_plot_sup)
        ax.invert_yaxis()
        #ax.set_yticks(np.arange(0,prof_plot_sup+100,100)) 
        ##ax.set_yticks(np.arange(0,prof_plot_sup+100,100),minor=True) 
        #ax.grid(axis='y',which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
        ax.axes.get_xaxis().set_visible(False)

        ax = plt.subplot(gs1[1])

        int_depth = np.transpose(locals()['int_depth_'+rodada[rod]])
        int_depth = np.ma.masked_where(data[0:int_depth.shape[0],:].mask,int_depth)
        plt.plot(lon_sec[ind_x:-1],int_depth[14,ind_x:-1], color='k', linestyle='--', linewidth=.5, alpha=.5)
        plt.plot(lon_sec[ind_x:-1],int_depth[20,ind_x:-1], color='k', linestyle='--', linewidth=.5, alpha=.5)

        ind_y = (depth_hycom>=prof_plot_sup) & (depth_hycom<=prof_plot_inf)
        plt.gca().patch.set_color('.9')
        plt.gca().patch.set_alpha(0.9)
        im = ax.contourf(lon_sec[ind_x:-1], depth_hycom[ind_y], data[ind_y,ind_x:-1],np.arange(min_mean,max_mean+inc_mean,inc_mean), \
                cmap=new_cmap,vmin=min_mean,vmax=max_mean,linestyles='dashed',extend='both')
        #im = ax.pcolormesh(lon_sec[ind_x:-1], depth_hycom[ind_y], data[ind_y,ind_x:-1], shading='interp', \
        #        cmap=new_cmap,vmin=min_mean,vmax=max_mean)

        ax.set_ylim(prof_plot_sup, prof_plot_inf)
        ax.invert_yaxis()

        #ax.set_yticks(np.arange(prof_plot_sup+200,prof_plot_inf+50,400)) 
        #ax.set_yticks(np.arange(prof_plot_sup+200,prof_plot_inf+50,200),minor=True) 
        #ax.grid(axis='y',which='major', color='k', linestyle='--', linewidth=.5, alpha=.5)
        leg = str(abs(lats[ll]))+s_lat_p2+' - '+legenda[rod]
        ax.text(lon_sec[ind_x+1], prof_plot_inf-400, leg, horizontalalignment='left', fontsize=font_size, fontweight='bold')

    #cbar_ax = fig.add_axes([0.92, 0.010, 0.02, 0.85])
    cbar_ax = fig.add_axes([0.91, 0.110, 0.02, 0.77])
    cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')        
    if len(months)==12:
       fig.savefig(out_dir+'/V_VEL_ROT_'+all_rod_save+'_'+data_inicial.strftime("%Y%m%d")+'-'+data_final.strftime("%Y%m%d")+ \
                   '_'+str(abs(lats[ll]))+s_lat_p2+'.png',dpi=resolution_plot,transparent=False,bbox_inches='tight',pad_inches=0.05)
    else:
       fig.savefig(out_dir+'/V_VEL_ROT_'+all_rod_save+'_'+'_months_'+str(months[0]).zfill(2)+'-'+str(months[-1]).zfill(2)+'_YEARS_'+ \
                   str(data_inicial.year)+'-'+str(data_final.year)+'_'+str(abs(lats[ll]))+s_lat_p2+'.png',dpi=resolution_plot, \
                   transparent=False,bbox_inches='tight',pad_inches=0.05)

    del u_rot_all, v_rot_all, int_depth_all, temp_all, sal_all, lat_sec, lon_sec, depth_hycom
exit()
    

