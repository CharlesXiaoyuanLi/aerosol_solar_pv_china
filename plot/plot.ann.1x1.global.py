from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
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

m = Basemap(projection='cyl',lat_0=0.,lon_0=180.,\
            llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=359.9999,\
            resolution='i',area_thresh=1e4)
m.drawcoastlines(linewidth=0.5)
m.drawcountries()

ny = len(lat); nx = len(lon)
lons, lats = m.makegrid(nx,ny)
x, y = m(lons, lats)
# draw filled contours
clevs = np.linspace(-140,140,15,endpoint=True)
#clevs = [-250,-100,-90,-80,-70,-60,-50,-40,-30,-20,\
#         -10,0,10,20,30,40,50,60,70,80,90,100,250]
cs = m.contourf(x,y,diff[varchosen],clevs,cmap=cm.RdYlGn,extend="both") # IMPORTANT: extend="both" enables contour for data out of range of the colorbar
# add colorbar
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('$W/m^2$')

# add title
plt.title('$F_{'+varchosen+',noaero} - F_{'+varchosen+',all}$')
plt.savefig("diff_"+varchosen+"_2014.ps",format="ps")

del(diff)

plt.show()
