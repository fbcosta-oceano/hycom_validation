#!/usr/bin/python

from datetime import *

def datespan(startDate, endDate, delta=timedelta(days=1)):
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta

def rmsd(a,b): #a:sample b:estimated value
 from numpy import mean, sqrt, arange, isnan, logical_not, isfinite

 a=a.flatten()
 b=b.flatten()
 mask = isfinite(a) & isfinite(b)
 a=a[mask]
 b=b[mask]
 diff=a-b
 value = sqrt(mean(diff**2))
 return value

def nanmean(a):
 from numpy import mean, isfinite
 a=a.flatten()
 mask = isfinite(a)
 a=a[mask]
 value = mean(a)
 return value

def interp2d(X,Y,data,XI,YI):
 from scipy import interpolate

 I = interpolate.RectBivariateSpline(X,Y,data,kx=1,ky=1)
 interp_data = I(XI,YI)

 return interp_data

def get_data_model(path,varname,lat,lon,date,level=0):
  from scipy import interpolate
  from netCDF4 import Dataset
  import numpy as np
  import gc
  import os.path

  year=date.strftime("%Y")
  jul=date.strftime("%j")
  sufix=""
  data=[[1e30]]
  if(varname=='temperature' or varname=='salinity'):
    sufix="3zt"
  elif(varname=='ssh'):
    sufix="fsd"
  elif(varname=='u' or varname=='v' or varname=='speed'):
    sufix="3zu"
  else:
    return data

  ncfile="%s/archv.%s_%s_00_%s.nc" % (path,year,jul,sufix)
  if(not os.path.isfile(ncfile)):
    return data

  model = Dataset(ncfile)
  lons = model.variables['Longitude'][:]
  lats = model.variables['Latitude'][:]
  depths = model.variables['Depth'][:]
  idz=np.nonzero(depths == level)[0][0]
  if(varname=='temperature' or varname=='salinity'):
    var=model.variables[varname][0,idz,:,:]
  elif(varname=='ssh'):
    var=model.variables[varname][0,:,:]
  elif(varname=='u' or varname=='v' or varname=='speed'):
    var=model.variables[varname][0,idz,:,:]
  model.close()

  I = interpolate.RectBivariateSpline(lats,lons,var,kx=1,ky=1)
  data = I(lat,lon)
  del I, lons, lats, var, sufix, idz, year, jul
  gc.collect()
  return data

def plot_section(X,Y,data,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',interp='linear',origin='upper',show=True,
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',aspect='auto',extend='both',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0):
  import gc
  import os
  import matplotlib.pyplot as plot
  import numpy as np

  ext=[np.min(X),np.max(X),np.min(Y),np.max(Y)]
  CF = plot.contourf(X, Y, data, contour_levels,origin=origin,extend=extend,cmap=cmap,vmin=cbmin,vmax=cbmax,extent=ext)
  CF.cmap.set_under(under)
  CF.cmap.set_over(over)
  CF.cmap.set_bad(bad)
  CL = plot.contour(X,Y,data,contour_levels,colors=('k',),linewidths=(1,),origin=origin,extent=ext,linestyles='dashed')
  cbar = plot.colorbar(CF,extend=extend,)
  cbar.set_label(um)
  cbar.ax.tick_params(labelsize=20)
  plot.clabel(CL, fmt = '%2.1f', colors = 'k', fontsize=18)
  plot.title(title,fontsize=25)
  plot.gca().invert_yaxis()
  plot.rc('xtick', labelsize=20)
  plot.rc('ytick', labelsize=20)
  if(show):
    plot.show()
  else:
    os.system('rm -f '+filename+' 2> /dev/null')
    plot.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  plot.close()
  del CF, CL, X, Y, data
  gc.collect()

def plot_points(dataX,dataY,vecdata,coordsX,coordsY,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',interp='linear',origin='upper',show=True,
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',aspect='auto',extend='both',plot_points=False,render='PS',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0):
  import gc
  import os
  import matplotlib.pyplot as plot
  #matplotlib.use(render)
  import numpy as np
  from matplotlib.mlab import griddata

  X = coordsX 
  Y = coordsY
  data = griddata(dataX, dataY, vecdata, X, Y, interp='linear')

  ext=[np.min(X),np.max(X),np.min(Y),np.max(Y)]
  CF = plot.contourf(X, Y, data, contour_levels,origin=origin,extend=extend,cmap=cmap,vmin=cbmin,vmax=cbmax,extent=ext)
  CF.cmap.set_under(under)
  CF.cmap.set_over(over)
  CF.cmap.set_bad(bad)
  CL = plot.contour(X,Y,data,contour_levels,colors=('k',),linewidths=(1,),origin=origin,extent=ext,linestyles='dashed')
  cbar = plot.colorbar(CF,extend=extend,)
  cbar.set_label(um)
  #plot.clabel(CL, fmt = '%2.1f', colors = 'k', fontsize=10)
  cbar.ax.tick_params(labelsize=20) # fontsize colorbar
  plot.clabel(CL, fmt = '%2.1f', colors = 'k', fontsize=18) 
  #plot.title(title,fontsize=25)
  plot.rc('xtick', labelsize=20)
  plot.rc('ytick', labelsize=20)

  if(plot_points):
    plot.scatter(dataX, dataY, marker='*', c='k', s=2, zorder=10)
  plot.title(title,fontsize=25)
  plot.gca().invert_yaxis()
  if(show):
    plot.show()
  else:
    os.system('rm -f '+filename+' 2> /dev/null')
    plot.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  plot.close()
  del CF, CL, X, Y, data
  gc.collect()

def pcolor(lon,lat,data):
  
  #import matplotlib
  #from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plt
  import numpy as np

  ax=plt.figure
  plt.pcolormesh(lon,lat,data);
  ax.show()

def plot_map(lat,lon,data,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat',render='PS',
             filename='figure.png',resolution=600,colormap='jet',under='#000033',show=False,
             over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0,
             pointsX=[],pointsY=[],pointSize=2,pointColor='k',pointMarker='*',subtitle=''):

  import gc
  import os
  import matplotlib
  #matplotlib.use(render)
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  fig = plot.figure()

  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat)+0.1,latstep)
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize=10,linewidth=1)
  meridians = np.arange(int(minlon),int(maxlon)+0.1,lonstep)
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize=10,linewidth=1)
  map_ax.drawcountries()

  x, y = map_ax(*np.meshgrid(lon, lat))
  cmap.set_under(under)
  cmap.set_over(over)
  cmap.set_bad(bad)
  if(contour):
     CS = map_ax.contour(x,y,data,30,levels=contour_levels,linewidths=1,colors='k')
     plot.clabel(CS, fontsize=4, inline=1)
  im = map_ax.pcolormesh(x,y,data,shading=shad,cmap=cmap,vmin=cbmin,vmax=cbmax)
  cbar = map_ax.colorbar(im,pad='5%',extend='both',)
  cbar.set_label(um,fontsize=10)
  plot.title(title,fontsize=10)
  px, py = map_ax(pointsX, pointsY)
  map_ax.plot(px, py, marker=pointMarker, c=pointColor, markersize=pointSize)
  plot.annotate(subtitle, xy=(0.08, 0.9), xycoords='axes fraction',fontsize=10)
  cbar.ax.tick_params(labelsize=10)
  if(show):
     plot.show()
  else:
     os.system('rm -f '+filename+' 2> /dev/null')
     fig.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  fig.clf()
  plot.close()
  del map_ax, im, lat, lon, data
  gc.collect()
                                                                                             

def plot_map_show(lat,lon,data,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat'):
  
  import matplotlib
  #matplotlib.use('TkAgg',warn=False)

  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plt
  import numpy as np

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  #fig = plot.figure()
  fig = plt.figure()
  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat),20.)
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  meridians = np.arange(int(minlon),int(maxlon),20.)
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')
  map_ax.drawcountries()

  x, y = map_ax(*np.meshgrid(lon, lat))
  cmap=plt.cm.get_cmap('jet',20)
  cmap.set_under('#000033')
  cmap.set_over('#330000')
  cmap.set_bad('#FFFFFF')
  im = map_ax.pcolormesh(x,y,data,shading=shad,cmap=cmap,vmin=cbmin,vmax=cbmax)
  cbar = map_ax.colorbar(im,pad='5%',extend='both')
  cbar.set_label(um)
  #mng = plt.get_current_fig_manager()
  #mng.resize(*mng.window.maxsize())
  fig.suptitle(title)
  fig.show()

def plot_map_save(lat,lon,data,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat',
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0):
  import gc
  import os
  import matplotlib
  matplotlib.use('TkAgg',warn=False)

  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  #matplotlib.rc('font', **{'sans-serif' : 'Arial','family' : 'sans-serif'})
  #matplotlib.rcParams.update({'font.size': 18})

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  fig = plot.figure()

  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat)+0.1,latstep)
  #map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize=10,linewidth=1)
  meridians = np.arange(int(minlon),int(maxlon)+0.1,lonstep)
  #map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize=10,linewidth=1)
  map_ax.drawcountries()

  x, y = map_ax(*np.meshgrid(lon, lat))
  #cmap=plot.cm.get_cmap(colormap,cint)
  cmap.set_under(under)
  cmap.set_over(over)
  cmap.set_bad(bad)
  #if(contour):
  CS = map_ax.contour(x,y,data,[27],linewidths=1,colors='k')
  #plot.clabel(CS, fontsize=4, inline=1)
  im = map_ax.pcolormesh(x,y,data,shading=shad,cmap=cmap,vmin=cbmin,vmax=cbmax)
  cbar = map_ax.colorbar(im,pad='5%',extend='both',)
  cbar.set_label(um)
  cbar.ax.tick_params(labelsize=13)
  plot.title(title,fontsize=25)
  #mng = plot.get_current_fig_manager()
  #mng.resize(*mng.window.maxsize())
  os.system('rm -f '+filename+' 2> /dev/null')
  fig.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  fig.clf()
  plot.close()
  del map_ax, im, lat, lon, data
  gc.collect()

def satelite_view(lat,lon,proj='cyl'):
  import matplotlib
  matplotlib.use('TkAgg',warn=False)
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)
  
  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.bluemarble()
  parallels = np.arange(int(minlat),int(maxlat),20.)
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  meridians = np.arange(int(minlon),int(maxlon),20.)
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')

  plot.show()
  plot.savefig('teste.png',dpi=resolution)
  plot.clf()

def plot_map_stream_save(lat,lon,U,V,field,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat',
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False):
  import gc
  import os
  import matplotlib
  matplotlib.use('TkAgg',warn=False)

  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  fig = plot.figure()

  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat),20.)
  #map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize=18,linewidth=0)
  meridians = np.arange(int(minlon),int(maxlon),30.)
  #map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize=18,linewidth=0)
  map_ax.drawcountries()

  x, y = map_ax(*np.meshgrid(lon, lat))
  #cmap=plot.cm.get_cmap(colormap,cint)
  cmap.set_under(under)
  cmap.set_over(over)
  cmap.set_bad(bad)

  map_ax.streamplot(x,y,U,V,color=field,linewidth=2,density=2,cmap=cmap)

  cbar = map_ax.colorbar(im,pad='5%',extend='both',)
  cbar.set_label(um)
  plot.title(title,fontsize=25)
  os.system('rm -f '+filename+' 2> /dev/null')
  fig.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  fig.clf()
  plot.close()
  del map_ax, im, lat, lon, U, V, field
  gc.collect()


def plot_map_speed_save(lat,lon,U,V,speed,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat',
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0):
  
  import gc
  import os
  import matplotlib
  matplotlib.use('TkAgg',warn=False)

  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  #matplotlib.rc('font', **{'sans-serif' : 'Arial','family' : 'sans-serif'})
  #matplotlib.rcParams.update({'font.size': 18})

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  fig = plot.figure()

  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat)+0.1,latstep)
  #map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize=10,linewidth=1) # 9,0.5
  meridians = np.arange(int(minlon),int(maxlon)+0.1,lonstep)
  #map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize=10,linewidth=1)
  map_ax.drawcountries()
  x, y = map_ax(*np.meshgrid(lon, lat))
  #cmap=plot.cm.get_cmap(colormap,cint)
  cmap.set_under(under)
  cmap.set_over(over)
  cmap.set_bad(bad)
  if(contour):
    CS = map_ax.contour(x,y,speed,30,levels=contour_levels,linewidths=1,colors='k')
    plot.clabel(CS, fontsize=4, inline=1)
  im = map_ax.pcolormesh(x,y,speed,shading=shad,cmap=cmap,vmin=cbmin,vmax=cbmax)


  uproj,vproj,xx,yy = \
  map_ax.transform_vector(U,V,lon,lat,50,80,returnxy=True,masked=True)
  Q = map_ax.quiver(xx[::3],yy[::3],uproj[::3],vproj[::3],scale=3.0,pivot='mid',scale_units='inches')
  # make quiver key.
  qk = plot.quiverkey(Q, 0.1, 0.6, 1.0, '1 m/s', labelpos='S')

  cbar = map_ax.colorbar(im,pad='5%',extend='both',)
  cbar.set_label(um)
  cbar.ax.tick_params(labelsize=13)
  plot.title(title,fontsize=25)
  #mng = plot.get_current_fig_manager()
  #mng.resize(*mng.window.maxsize())
  os.system('rm -f '+filename+' 2> /dev/null')
  fig.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  fig.clf()
  plot.close()
  del map_ax, im, lat, lon, U, V, speed
  gc.collect()


def plot_map_speed_save_line(lat,lon,U,V,speed,pointsX,pointsY,cmap,title='',um='',cbmin=0.0,cbmax=38,proj='cyl',shad='flat',
                  filename='figure.png',resolution=600,colormap='jet',under='#000033',
                  over='#330000',bad='#FFFFFF',cint=20,trans=False,contour=False,contour_levels=[],lonstep=1.0,latstep=1.0):
  
  import gc
  import os
  import matplotlib
  matplotlib.use('TkAgg',warn=False)

  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  #matplotlib.rc('font', **{'sans-serif' : 'Arial','family' : 'sans-serif'})
  #matplotlib.rcParams.update({'font.size': 18})

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)

  fig = plot.figure()

  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.drawcoastlines()
  map_ax.fillcontinents(lake_color='aqua')
  parallels = np.arange(int(minlat),int(maxlat)+0.1,latstep)
  #map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize=10,linewidth=1) # 9,0.5
  meridians = np.arange(int(minlon),int(maxlon)+0.1,lonstep)
  #map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize=10,linewidth=1)
  map_ax.drawcountries()
  x, y = map_ax(*np.meshgrid(lon, lat))
  #cmap=plot.cm.get_cmap(colormap,cint)
  cmap.set_under(under)
  cmap.set_over(over)
  cmap.set_bad(bad)
  #if(contour):
  im = map_ax.pcolormesh(x,y,speed,shading=shad,cmap=cmap,vmin=cbmin,vmax=cbmax)
  CS = map_ax.contour(x,y,speed,[27],linewidths=1,colors='k')
  #plot.clabel(CS, fontsize=4, inline=1)


  uproj,vproj,xx,yy = \
  map_ax.transform_vector(U,V,lon,lat,50,80,returnxy=True,masked=True)
  Q = map_ax.quiver(xx[::3],yy[::3],uproj[::3],vproj[::3],scale=3.0,pivot='mid',scale_units='inches')
  # make quiver key.
  qk = plot.quiverkey(Q, 0.1, 0.6, 1.0, '1 m/s', labelpos='S')

  cbar = map_ax.colorbar(im,pad='5%',extend='both',)
  cbar.set_label(um)
  cbar.ax.tick_params(labelsize=13)
  plot.title(title,fontsize=25)
  px, py = map_ax(pointsX, pointsY)
  #map_ax.plot(px, py, marker='.', c='k', markersize=2,linewidth=2)
  #mng = plot.get_current_fig_manager()
  #mng.resize(*mng.window.maxsize())
  os.system('rm -f '+filename+' 2> /dev/null')
  fig.savefig(filename,dpi=resolution,transparent=trans,bbox_inches='tight',pad_inches=0)
  fig.clf()
  plot.close()
  del map_ax, im, lat, lon, U, V, speed
  gc.collect()


def etopo_view(lat,lon,proj='cyl'):
  import matplotlib
  matplotlib.use('TkAgg',warn=False)
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plot
  import numpy as np

  minlat=min(lat)
  maxlat=max(lat)
  minlon=min(lon)
  maxlon=max(lon)
  
  map_ax = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,projection=proj)
  map_ax.etopo()
  parallels = np.arange(int(minlat),int(maxlat),20.)
  map_ax.drawparallels(parallels,labels=[True,False,False,False],fontsize='small')
  meridians = np.arange(int(minlon),int(maxlon),20.)
  map_ax.drawmeridians(meridians,labels=[False,False,False,True],fontsize='small')

  plot.show()
  plot.savefig('teste.png',dpi=600)
  plot.clf()
