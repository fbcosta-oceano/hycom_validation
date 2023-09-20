from netCDF4 import Dataset
import remo
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
jet = matplotlib.cm.get_cmap('jet')


#model_file='/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc/ATLi0.08/expt_01.1/archv.2009_004_00_fsd.nc'
model_file='/scratch/rodas/filipe.costa/previsao/hycom_2_2_18/proc/ATLj0.04/expt_00.1/output/Ncdf/2009010800/archv.2009_009_00_fsd.nc'

nc = Dataset(model_file, 'r') #abre o netcdf
latmod=nc.variables['Latitude'][:]
lonmod=nc.variables['Longitude'][:]
sshmod=nc.variables['ssh'][:]
#sshmod=sshmod.transpose()
print(latmod.shape)
print(lonmod.shape)
print(sshmod.shape)

remo.plot_map(latmod, lonmod, sshmod[0,:,:],jet,cbmin=-.5,cbmax=1,lonstep=7.0,latstep=4.0,filename='SSH_archv.2009_009_00.png')
