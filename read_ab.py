import remo
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings("ignore")
jet = matplotlib.cm.get_cmap('jet')

IDM=1717
JDM=2345

#IDM=628
#JDM=780

#IDM=502
#JDM=564

IJDM=IDM*JDM
npad=4096-(IJDM%4096)

arch_name = 'archv.2012_086_00'
arch_dir = 'expt_02.2/output/ab/20120325'
arch_name = 'restart_2015d001h00'
arch_dir = 'expt_02.2/output/ab/20120325'
#arch_name = 'depth_ATLj0.04_01'
#arch_dir = 'topo'

f = open('/home/filipe.costa/previsao/hycom_2_2_18/proc/ATLd0.08/topo/regional.grid.a','rb')
f.seek(0*4*(npad+IJDM))
field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
lon = np.reshape(field,(JDM,IDM))
lon = lon[0,:]
f.seek(1*4*(npad+IJDM))
field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
lat = np.reshape(field,(JDM,IDM))
lat = lat[:,0]
f.close()

#lon, lat = np.meshgrid(lon, lat)

f = open('/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc/ATLd0.08/'+arch_dir+'/'+arch_name+'.a','rb')
#f.seek(11*4*(npad+IJDM)) #TEMP
f.seek(1*4*(IJDM+npad))
#field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
field = np.fromfile(f,dtype='>f4',count=IDM*JDM)
cond = abs(field) > 999999
field = np.ma.masked_where(cond,field)
field = np.reshape(field,(JDM,IDM))
field = field/9.806

#field=field.transpose()
#print(field[570:578,520:528])
f.close()

#remo.plot_map(lat, lon, field,jet,cbmin=15,cbmax=30,lonstep=3.0,latstep=2.0,filename='TEMP_archv.2009_009_00.png') 
#remo.plot_map(lat, lon, field,jet,cbmin=5,cbmax=30,lonstep=3.0,latstep=2.0,show=True,filename='TEMP_archv.2009_004_00_ATLd.png')
remo.plot_map(lat, lon, field,jet,cbmin=.2,cbmax=1.2,lonstep=7.0,latstep=4.0,filename='SSH_'+arch_name+'_new.png')
#remo.plot_map(lat, lon, field,jet,cbmin=.0,cbmax=5500,lonstep=7.0,latstep=4.0,filename=arch_name+'.png')
#fig = plt.figure(figsize=(10, 8))
#m = Basemap(projection='merc', resolution='l',llcrnrlon=-54,llcrnrlat=-34,urcrnrlon=-32,urcrnrlat=-12)
##m.shadedrelief(scale=0.5)
##m.pcolormesh(lon, lat, field, latlon=True)
#m.pcolormesh(lon, lat, field, shading='flat', cmap='jet')
#plt.savefig('test')

