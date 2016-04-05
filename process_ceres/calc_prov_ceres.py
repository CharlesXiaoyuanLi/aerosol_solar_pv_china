#!/bin/python2

### Calculate surface albedo

import numpy as np
from netCDF4 import Dataset
import sys


### Define constants

re  = 6.37122E+6
d2r = np.pi/180.


### Read data

nc_file = "/storage/CERES/CERES_SFCFLUX/SYN1deg-3H/ceres_syn1deg-3h_2014.nc"

fh = Dataset(nc_file, 'r')

sfc_direct_all  = fh.variables['sfc_comp_sw_direct_all_3h'][:][:][:]
sfc_direct_clr  = fh.variables['sfc_comp_sw_direct_clr_3h'][:][:][:]
sfc_direct_pri  = fh.variables['sfc_comp_sw_direct_pri_3h'][:][:][:]
sfc_direct_noaero  = fh.variables['sfc_comp_sw_direct_noaero_3h'][:][:][:]

sfc_diffuse_all  =fh.variables['sfc_comp_sw_diffuse_all_3h'][:][:][:]
sfc_diffuse_clr  =fh.variables['sfc_comp_sw_diffuse_clr_3h'][:][:][:]
sfc_diffuse_pri  =fh.variables['sfc_comp_sw_diffuse_pri_3h'][:][:][:]
sfc_diffuse_noaero  =fh.variables['sfc_comp_sw_diffuse_noaero_3h'][:][:][:]

lat  = fh.variables['lat'][:]
lon  = fh.variables['lon'][:]
time = fh.variables['time'][:]

nlat = len(lat)
nlon = len(lon)
ntime = len(time)

fh.close()

latu = lat + 0.5
latl = lat - 0.5

### Calculate area of each grid
weight = np.sin(latu*d2r) - np.sin(latl*d2r)
area = 2. * np.pi * np.power(re,2) * weight / nlon
area_globe = np.array([ area for x in range(nlon) ]).transpose()

regions = ["Beijing","Tianjin","Hebei","Shanxi","Neimenggu","Liaoning","Jilin",\
           "Heilongjiang","Shanghai","Jiangsu","Zhejiang","Anhui","Fujian",\
           "Jiangxi","Shandong","Henan","Hubei","Hunan","Guangdong","Guangxi",\
           "Hainan","Chongqing","Sichuan","Guizhou","Yunnan","Xizang","Shaanxi",\
           "Gansu","Qinghai","Ningxia","Xinjiang"]

regionids = range(1,32) 

### Read mask ncfile
mask_file = "../../maps_masks_shapefiles/china_provinces_"+str(nlat)+"x"+str(nlon)+".nc"
fh_mask = Dataset(mask_file, 'r')
rgmask = fh_mask.variables['mask'][:][:]

fh_mask.close()

prov_sfc_direct_all = np.zeros((len(regions),ntime))
prov_sfc_direct_clr = np.zeros((len(regions),ntime))
prov_sfc_direct_pri = np.zeros((len(regions),ntime))
prov_sfc_direct_noaero = np.zeros((len(regions),ntime))

prov_sfc_diffuse_all = np.zeros((len(regions),ntime))
prov_sfc_diffuse_clr = np.zeros((len(regions),ntime))
prov_sfc_diffuse_pri = np.zeros((len(regions),ntime))
prov_sfc_diffuse_noaero = np.zeros((len(regions),ntime))

### Calculate average for each region
for region, provid in zip(regions, regionids):

    print region
    rgind = regions.index(region)

    ilats, ilons = np.where( rgmask == provid )

    print ilats, ilons

    prov_sfc_direct_all_t = np.array([ sfc_direct_all[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_direct_clr_t = np.array([ sfc_direct_clr[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_direct_pri_t = np.array([ sfc_direct_pri[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_direct_noaero_t = np.array([ sfc_direct_noaero[:,x,y] for x,y in zip(ilats,ilons) ])

    prov_sfc_diffuse_all_t = np.array([ sfc_diffuse_all[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_diffuse_clr_t = np.array([ sfc_diffuse_clr[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_diffuse_pri_t = np.array([ sfc_diffuse_pri[:,x,y] for x,y in zip(ilats,ilons) ])
    prov_sfc_diffuse_noaero_t = np.array([ sfc_diffuse_noaero[:,x,y] for x,y in zip(ilats,ilons) ])

    print prov_sfc_direct_all_t

    prov_area_array = np.array([ area_globe[x,y] for x,y in zip(ilats,ilons) ])

    if (region == 'Shanghai'):
        continue

    prov_sfc_direct_all[provid-1,:] = np.average(prov_sfc_direct_all_t, axis=0, weights=prov_area_array)
    prov_sfc_direct_clr[provid-1,:] = np.average(prov_sfc_direct_clr_t, axis=0, weights=prov_area_array)
    prov_sfc_direct_pri[provid-1,:] = np.average(prov_sfc_direct_pri_t, axis=0, weights=prov_area_array)
    prov_sfc_direct_noaero[provid-1,:] = np.average(prov_sfc_direct_noaero_t, axis=0, weights=prov_area_array)

    prov_sfc_diffuse_all[provid-1,:] = np.average(prov_sfc_diffuse_all_t,\
                                                  axis=0, weights=prov_area_array)
    prov_sfc_diffuse_clr[provid-1,:] = np.average(prov_sfc_diffuse_clr_t,\
                                                  axis=0, weights=prov_area_array)
    prov_sfc_diffuse_pri[provid-1,:] = np.average(prov_sfc_diffuse_pri_t,\
                                                  axis=0, weights=prov_area_array)
    prov_sfc_diffuse_noaero[provid-1,:] = np.average(prov_sfc_diffuse_noaero_t,\
                                                  axis=0, weights=prov_area_array)
 
    print prov_sfc_direct_all

### Write data to ncfile
outfile = Dataset('prov_sfc_sw_3H_2014.nc','w',format='NETCDF4')

outfile.createDimension('province',len(regions))
outfile.createDimension('time',ntime)

time_dim = outfile.createVariable('time',np.float32,('time',))

region_dim = outfile.createVariable('province',np.int32,('province',))
description = ''
for region, rgid in zip(regions, regionids):
    description += region + ': ' + str(rgid) + ',  '
region_dim.description =description

prov_sfc_direct_all_out = outfile.createVariable('sfc_direct_all',np.float32,('province','time'))
prov_sfc_direct_clr_out = outfile.createVariable('sfc_direct_clr',np.float32,('province','time'))
prov_sfc_direct_pri_out = outfile.createVariable('sfc_direct_pri',np.float32,('province','time'))
prov_sfc_direct_noaero_out = outfile.createVariable('sfc_direct_noaero',np.float32,('province','time'))

prov_sfc_diffuse_all_out = outfile.createVariable('sfc_diffuse_all',np.float32,('province','time'))
prov_sfc_diffuse_clr_out = outfile.createVariable('sfc_diffuse_clr',np.float32,('province','time'))
prov_sfc_diffuse_pri_out = outfile.createVariable('sfc_diffuse_pri',np.float32,('province','time'))
prov_sfc_diffuse_noaero_out = outfile.createVariable('sfc_diffuse_noaero',np.float32,('province','time'))

region_dim[:]   = regionids
time_dim[:]     = time

prov_sfc_direct_all_out[:] = prov_sfc_direct_all 
prov_sfc_direct_clr_out[:] = prov_sfc_direct_clr 
prov_sfc_direct_pri_out[:] = prov_sfc_direct_pri 
prov_sfc_direct_noaero_out[:] = prov_sfc_direct_noaero 

prov_sfc_diffuse_all_out[:] = prov_sfc_diffuse_all 
prov_sfc_diffuse_clr_out[:] = prov_sfc_diffuse_clr 
prov_sfc_diffuse_pri_out[:] = prov_sfc_diffuse_pri 
prov_sfc_diffuse_noaero_out[:] = prov_sfc_diffuse_noaero 

outfile.close()


