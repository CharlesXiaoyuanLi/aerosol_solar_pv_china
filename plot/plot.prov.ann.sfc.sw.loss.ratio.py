from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
import numpy as np
from matplotlib import cm, colors
import matplotlib as mpl
import sys
from netCDF4 import Dataset


fig, ax = plt.subplots(frameon=False)
#ax = fig.add_subplot(111)

map = Basemap(llcrnrlon=80.,llcrnrlat=12.5,urcrnrlon=140.,urcrnrlat=52.5,resolution='i',projection='laea',lat_0=35., lon_0=110.)

#map.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='#ddaa66') #lake_color='aqua'
#map.drawcoastlines()

map.readshapefile('/home/xl4/maps_masks_shapefiles/china_provinces/china_provinces00/china_provinces00','cnprovinces', drawbounds = True)

patches = {}
subpatches = []
count = 0
for info, shape in zip(map.cnprovinces_info, map.cnprovinces):
    if count > 31:
        break
    if info['RINGNUM'] == 1:
        if info['SHAPENUM'] != 1: 
            patches[pronm] = subpatches
            subpatches = []
        pronm = info['EPROV']
        count += 1
    subpatches.append( Polygon(np.array(shape), True) )

print count

### China Province Namelist
pro_nm = ["Beijing","Tianjin","Hebei","Shanxi","Neimenggu","Liaoning","Jilin",\
           "Heilongjiang","Shanghai","Jiangsu","Zhejiang","Anhui","Fujian",\
           "Jiangxi","Shandong","Henan","Hubei","Hunan","Guangdong","Guangxi",\
           "Hainan","Chongqing","Sichuan","Guizhou","Yunnan","Xizang","Shaanxi",\
           "Gansu","Qinghai","Ningxia","Xinjiang"]

proids = range(1,32) 

### Read in surface sw radiation data
nc_file = "/home/xl4/aerosol_solar_pv_china/annual/prov_annual_sfc_sw_31.nc"

fh = Dataset(nc_file, 'r')

sfc_direct_all = fh.variables['sfc_direct_all'][:]
sfc_direct_clr = fh.variables['sfc_direct_clr'][:]
sfc_direct_pri = fh.variables['sfc_direct_pri'][:]
sfc_direct_noaero = fh.variables['sfc_direct_noaero'][:]

sfc_diffuse_all = fh.variables['sfc_diffuse_all'][:]
sfc_diffuse_clr = fh.variables['sfc_diffuse_clr'][:]
sfc_diffuse_pri = fh.variables['sfc_diffuse_pri'][:]
sfc_diffuse_noaero = fh.variables['sfc_diffuse_noaero'][:]

### Calculate Ratio
sfc_direct_cld_ratio = 100. * (sfc_direct_clr - sfc_direct_all) / sfc_direct_all
sfc_direct_aer_ratio = 100. *(sfc_direct_noaero - sfc_direct_all) / sfc_direct_all

sfc_diffuse_cld_ratio = 100. * (sfc_diffuse_clr - sfc_diffuse_all) / sfc_diffuse_all
sfc_diffuse_aer_ratio = 100. * (sfc_diffuse_noaero - sfc_diffuse_all) / sfc_diffuse_all

sfc_global_cld_ratio = 100. * (sfc_direct_clr + sfc_diffuse_clr - sfc_direct_all -
                        sfc_diffuse_all) / (sfc_direct_all + sfc_diffuse_all)
sfc_global_aer_ratio = 100. * (sfc_direct_noaero + sfc_diffuse_noaero - sfc_direct_all -
                        sfc_diffuse_all) / (sfc_direct_all + sfc_diffuse_all)

### Create color dictionary that maps values of each province to colorbar

color={}

for i in range(len(pro_nm)):
    color[pro_nm[i]] = sfc_direct_aer_ratio[i]
color['Shanghai'] = color['Jiangsu']

#cmap=cm.YlOrRd
cmap=cm.YlGnBu_r

#bounds = np.linspace(0,140,8)
bounds = np.linspace(-30,0,7)
norm = colors.BoundaryNorm(bounds, cmap.N)

#bounds = np.linspace(0,120,21)
#norm = colors.BoundaryNorm(bounds, ncolors=256)
#print norm(100)

#sys.exit()

for key, value in patches.iteritems():
    print key+" "+str(color[key])
#    ax.add_collection(PatchCollection(value,\
#        facecolor=cm.YlOrRd(norm(color[key])),\
#        cmap=cm.YlOrRd, edgecolor='k', linewidth=1., zorder=2))
    ax.add_collection(PatchCollection(value,\
        facecolor=cm.YlGnBu_r(norm(color[key])),\
        cmap=cm.YlGnBu_r, edgecolor='k', linewidth=1., zorder=2))



#    ax.add_collection(PatchCollection(value, facecolors=color[key], cmap=cm.RdYlBu, edgecolor='k',
#        linewidth=1., zorder=2, norm=norm))


### Plot colorbar
#set room for colorbar
fig.subplots_adjust(right = 0.85)


axcb = fig.add_axes([0.87,0.15,0.03,0.7])

cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, boundaries=bounds,
        orientation='vertical', ticks=bounds)

cb.set_label('%')

plt.savefig('prov_sfc_sw_diffuse_loss_rmvaero_ratio.ps',format='ps')

plt.show()
