from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib import cm

print 'Please choose between direct, diffuse and global ...'
varchosen = raw_input('--> ')

### Read in data
ncfilnm = "/tigress/xl4/CERES_SFCFLUX/SYN1deg-3H/ceres_syn1deg-annual_2014.nc"
ncfh    = Dataset(ncfilnm,'r')

sfc_direct_all     = ncfh.variables['sfc_comp_sw_direct_all_3h'][0][:][:]
sfc_direct_noaero  = ncfh.variables['sfc_comp_sw_direct_noaero_3h'][0][:][:]
sfc_diffuse_all    = ncfh.variables['sfc_comp_sw_diffuse_all_3h'][0][:][:]
sfc_diffuse_noaero = ncfh.variables['sfc_comp_sw_diffuse_noaero_3h'][0][:][:]

lat = ncfh.variables['lat'][:]
lon = ncfh.variables['lon'][:]

### Calculate Data
diff_direct = sfc_direct_noaero - sfc_direct_all
diff_diffuse = sfc_diffuse_noaero - sfc_diffuse_all
diff_global = diff_direct+diff_diffuse

diff={}
diff['direct'] = diff_direct
diff['diffuse'] = diff_diffuse
diff['global'] = diff_global

del(diff_direct)
del(diff_diffuse)
del(diff_global)

### Plot
fig = plt.figure(figsize=(8.0,5.0))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])

llcrnrlat=-90.; urcrnrlat=90
llcrnrlon=-180; urcrnrlon=180
lat_0 = (llcrnrlat + urcrnrlat) / 2.
lon_0 = (llcrnrlon + urcrnrlon) / 2.

m = Basemap(projection='cea',lat_0=lat_0,lon_0=lon_0,\
            llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,\
            llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,\
            resolution='l',area_thresh=1e6)
m.drawcoastlines(linewidth=0.5)
m.drawcountries()
m.drawparallels(np.arange(-60,90,30),labels=[1,0,0,0])
m.drawmeridians(np.arange(m.lonmin,m.lonmax+30,60),labels=[0,0,0,1])
#m.readshapefile('/home/xl4/maps_masks_shapefiles/world_countries_2008/world_countries_2008','world_countries_2008', 
#                drawbounds = True, linewidth=1.0)

#imshow
diff[varchosen],lon = shiftgrid(180.5,diff[varchosen],lon,start=False)
#diff[varchosen],lon = addcyclic(diff[varchosen],lon)

print lon

lons, lats = np.meshgrid(lon, lat)
x, y = m(lons, lats)

# imshow
#nx = int((m.xmax-m.xmin)/25000.)+1; ny = int((m.ymax-m.ymin)/25000.)+1
#print m.xmax,m.xmin,m.ymax,m.ymin
#print nx,ny
#diffplot = m.transform_scalar(diff[varchosen],lon,lat,nx,ny)

# draw filled contours
if varchosen == 'global':
    clevs = np.linspace(0,50,11,endpoint=True)
elif varchosen == 'direct':
    clevs = np.linspace(0,150,16,endpoint=True)
elif varchosen == 'diffuse':
    clevs = np.linspace(-100,0,11,endpoint=True)

if varchosen == 'diffuse':
    #cs = m.contourf(x,y,diff[varchosen],clevs,cmap=cm.YlGnBu_r,extend="both") # IMPORTANT: extend="both" enables contour for data out of range of the colorbar
    #cs = m.imshow(diffplot,vmin=clevs[0],vmax=clevs[-1],cmap=cm.YlGnBu_r)
    cs = m.pcolormesh(x,y,diff[varchosen],vmin=clevs[0],vmax=clevs[-1],cmap=cm.YlGnBu_r)

else:
    #cs = m.contourf(x,y,diff[varchosen],clevs,cmap=cm.YlOrRd,extend="both") # IMPORTANT: extend="both" enables contour for data out of range of the colorbar
    #cs = m.imshow(diffplot,vmin=clevs[0],vmax=clevs[-1],cmap=cm.YlOrRd)
    cs = m.pcolormesh(x,y,diff[varchosen],vmin=clevs[0],vmax=clevs[-1],cmap=cm.YlOrRd)

# add colorbar
cbar = m.colorbar(cs,location='bottom',pad="12%")
cbar.set_label('$W/m^2$')

# add title
plt.title('$F_{'+varchosen+',noaero} - F_{'+varchosen+',all}$')
plt.savefig("diff_"+varchosen+"_2014_wld.ps",format="ps")

del(diff)

plt.show()
