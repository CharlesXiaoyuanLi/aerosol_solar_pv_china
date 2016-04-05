#!/bin/python2

### Calculate surface albedo

import numpy as np
from netCDF4 import Dataset
import sys


### Define constants

re  = 6.37122E+6
d2r = np.pi/180.


### Prepare for reading data

mons = np.arange(1,13,1)

for mon in mons:

    ### Read data
    
    nc_file = "/storage/MERRA-2/tavg3H_1x1/merra2.tavg3h.1x1.2014{:0>2d}.nc".format(mon)
    
    fh = Dataset(nc_file, 'r')
    
    T2M    = fh.variables['T2M'][:][:][:]
    T10M   = fh.variables['T10M'][:][:][:]
    T2MDEW = fh.variables['T2MDEW'][:][:][:]
    T2MWET = fh.variables['T2MWET'][:][:][:]
    TS     = fh.variables['TS'][:][:][:]
    U2M    = fh.variables['U2M'][:][:][:]
    U10M   = fh.variables['U10M'][:][:][:]
    V2M    = fh.variables['V2M'][:][:][:]
    V10M   = fh.variables['V10M'][:][:][:]

    lat  = fh.variables['lat'][:]
    lon  = fh.variables['lon'][:]
    
    nlat = len(lat)
    nlon = len(lon)
    ntime = T2M.shape[0]
    
    #time = fh.variables['time'][:]
    
    fh.close()
    
    if mon == 1:
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
    
    prov_T2M    = np.zeros((len(regions),ntime))
    prov_T10M   = np.zeros((len(regions),ntime))
    prov_T2MDEW = np.zeros((len(regions),ntime))
    prov_T2MWET = np.zeros((len(regions),ntime))
    prov_TS     = np.zeros((len(regions),ntime))
    prov_U2M    = np.zeros((len(regions),ntime))
    prov_U10M   = np.zeros((len(regions),ntime))
    prov_V2M    = np.zeros((len(regions),ntime))
    prov_V10M   = np.zeros((len(regions),ntime))

    ### Calculate average for each region
    for region, provid in zip(regions, regionids):
    
        print region
        rgind = regions.index(region)
    
        ilats, ilons = np.where( rgmask == provid )
    
        print ilats, ilons
    
        prov_T2M_t    = np.array([ T2M[:,x,y] for x,y in zip(ilats,ilons) ])   
        prov_T10M_t   = np.array([ T10M[:,x,y] for x,y in zip(ilats,ilons) ])  
        prov_T2MDEW_t = np.array([ T2MDEW[:,x,y] for x,y in zip(ilats,ilons) ])
        prov_T2MWET_t = np.array([ T2MWET[:,x,y] for x,y in zip(ilats,ilons) ])
        prov_TS_t     = np.array([ TS[:,x,y] for x,y in zip(ilats,ilons) ])    
        prov_U2M_t    = np.array([ U2M[:,x,y] for x,y in zip(ilats,ilons) ])   
        prov_U10M_t   = np.array([ U10M[:,x,y] for x,y in zip(ilats,ilons) ])  
        prov_V2M_t    = np.array([ V2M[:,x,y] for x,y in zip(ilats,ilons) ])   
        prov_V10M_t   = np.array([ V10M[:,x,y] for x,y in zip(ilats,ilons) ])  

        prov_area_array = np.array([ area_globe[x,y] for x,y in zip(ilats,ilons) ])
    
        if (region == 'Shanghai'):
            continue
    
        prov_T2M[provid-1,:]    = np.average( prov_T2M_t , axis=0, weights=prov_area_array)  
        prov_T10M[provid-1,:]   = np.average( prov_T10M_t , axis=0, weights=prov_area_array) 
        prov_T2MDEW[provid-1,:] = np.average( prov_T2MDEW_t, axis=0, weights=prov_area_array)
        prov_T2MWET[provid-1,:] = np.average( prov_T2MWET_t, axis=0, weights=prov_area_array)
        prov_TS[provid-1,:]     = np.average( prov_TS_t , axis=0, weights=prov_area_array)   
        prov_U2M[provid-1,:]    = np.average( prov_U2M_t , axis=0, weights=prov_area_array)  
        prov_U10M[provid-1,:]   = np.average( prov_U10M_t , axis=0, weights=prov_area_array) 
        prov_V2M[provid-1,:]    = np.average( prov_V2M_t , axis=0, weights=prov_area_array)  
        prov_V10M[provid-1,:]   = np.average( prov_V10M_t , axis=0, weights=prov_area_array) 
    
    if mon == 1:
        prov_T2M_all    = prov_T2M
        prov_T10M_all   = prov_T10M  
        prov_T2MDEW_all = prov_T2MDEW
        prov_T2MWET_all = prov_T2MWET
        prov_TS_all     = prov_TS    
        prov_U2M_all    = prov_U2M   
        prov_U10M_all   = prov_U10M  
        prov_V2M_all    = prov_V2M   
        prov_V10M_all   = prov_V10M
    else:
        prov_T2M_all    = np.concatenate(( prov_T2M_all, prov_T2M), axis=1)
        prov_T10M_all   = np.concatenate(( prov_T10M_all, prov_T10M), axis=1)  
        prov_T2MDEW_all = np.concatenate(( prov_T2MDEW_all, prov_T2MDEW), axis=1)
        prov_T2MWET_all = np.concatenate(( prov_T2MWET_all, prov_T2MWET), axis=1)
        prov_TS_all     = np.concatenate(( prov_TS_all, prov_TS), axis=1)    
        prov_U2M_all    = np.concatenate(( prov_U2M_all, prov_U2M), axis=1)   
        prov_U10M_all   = np.concatenate(( prov_U10M_all, prov_U10M), axis=1)  
        prov_V2M_all    = np.concatenate(( prov_V2M_all, prov_V2M), axis=1)   
        prov_V10M_all   = np.concatenate(( prov_V10M_all, prov_V10M), axis=1)

    del(prov_T2M)
    del(prov_T10M)  
    del(prov_T2MDEW)
    del(prov_T2MWET)
    del(prov_TS)    
    del(prov_U2M)   
    del(prov_U10M)  
    del(prov_V2M)   
    del(prov_V10M)  

print prov_T2M_all.shape

### Create new time dimension variable
nc_file = "/storage/CERES/CERES_SFCFLUX/SYN1deg-3H/ceres_syn1deg-3h_2014.nc"
fh = Dataset(nc_file, 'r')
time = fh.variables['time'][:]
fh.close()

ntime = len(time)

if ntime != prov_T2M_all.shape[1]:
    print "shape does not match"
    sys.exit()

### Write data to ncfile
outfile = Dataset('prov_met_3H_2014.nc','w',format='NETCDF4')

outfile.createDimension('province',len(regions))
outfile.createDimension('time',ntime)

time_dim = outfile.createVariable('time',np.float32,('time',))

region_dim = outfile.createVariable('province',np.int32,('province',))
description = ''
for region, rgid in zip(regions, regionids):
    description += region + ': ' + str(rgid) + ',  '
region_dim.description =description

prov_T2M_out = outfile.createVariable('T2M',np.float32,('province','time'))
prov_T10M_out = outfile.createVariable('T10M',np.float32,('province','time'))
prov_T2MDEW_out = outfile.createVariable('T2MDEW',np.float32,('province','time'))
prov_T2MWET_out = outfile.createVariable('T2MWET',np.float32,('province','time'))
prov_TS_out = outfile.createVariable('TS',np.float32,('province','time'))
prov_U2M_out = outfile.createVariable('U2M',np.float32,('province','time'))
prov_U10M_out = outfile.createVariable('U10M',np.float32,('province','time'))
prov_V2M_out = outfile.createVariable('V2M',np.float32,('province','time'))
prov_V10M_out = outfile.createVariable('V10M',np.float32,('province','time'))

region_dim[:]   = regionids
time_dim[:]     = time

prov_T2M_out[:]    =  prov_T2M_all     
prov_T10M_out[:]   =  prov_T10M_all  
prov_T2MDEW_out[:] =  prov_T2MDEW_all
prov_T2MWET_out[:] =  prov_T2MWET_all
prov_TS_out[:]     =  prov_TS_all    
prov_U2M_out[:]    =  prov_U2M_all   
prov_U10M_out[:]   =  prov_U10M_all  
prov_V2M_out[:]    =  prov_V2M_all   
prov_V10M_out[:]   =  prov_V10M_all  

outfile.close()


