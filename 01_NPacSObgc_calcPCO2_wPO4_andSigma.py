# -*- coding: utf-8 -*-
"""
Created on Fri Sept 22 14:52:52 2023

@author: mgs23

2025-06-26
"""



# %% Set up

# Import modules

# Note, might also need scipy (v 1.11.3 for this code), netCDF4 (1.6.2), 
# and pandas (2.1.4) libraries installed in environment for xarray 
# to open .nc files properly
import xarray as xr
import numpy as np
import PyCO2SYS as pyco2
import gsw.conversions as gsw_c
import gsw.density as gsw_d

xr.set_options(keep_attrs=True)



# %% (1) Read in UVic Output

# %%% Read in UVic dic, alk, po4, T, S

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
x = xr.open_dataset(pathname+filename, decode_times=False) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years

dic_UV_pmoc  = x.O_dic.isel(time=-1)
alk_UV_pmoc  = x.O_alk.isel(time=-1)
po4_UV_pmoc  = x.O_po4.isel(time=-1)
potT_UV_pmoc = x.O_temp.isel(time=-1)
sal_UV_pmoc  = x.O_sal.isel(time=-1)

dic_UV_timeEvol  = x.O_dic
alk_UV_timeEvol  = x.O_alk
po4_UV_timeEvol  = x.O_po4
potT_UV_timeEvol = x.O_temp
sal_UV_timeEvol  = x.O_sal


pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg.02001.01.01.nc'
x = xr.open_dataset(pathname+filename, decode_times=False)
# time=1. array([2000.5]) (the end of the LGM simulation - our control)

dic_UV_ctrl  = x.O_dic.isel(time=0)
alk_UV_ctrl  = x.O_alk.isel(time=0)
po4_UV_ctrl  = x.O_po4.isel(time=0)
potT_UV_ctrl = x.O_temp.isel(time=0)
sal_UV_ctrl  = x.O_sal.isel(time=0)


# %%% Calc rho_UV 

# https://teos-10.github.io/GSW-Python/index.html
# gsw = Python implementation of the Thermodynamic Equation of Seawater 2010 (TEOS-10) 

# Will need abs salinity (good) & conservative temp (from my pot T) to calc
# sigma0 (Pot dens ref 0dbar)

import gsw

# (1) Convert Pot Temp to Conservative Temp

conT_UV_ctrl = gsw.conversions.CT_from_pt(sal_UV_ctrl, potT_UV_ctrl)
conT_UV_pmoc = gsw.conversions.CT_from_pt(sal_UV_pmoc, potT_UV_pmoc)
conT_UV_timeEvol = gsw.conversions.CT_from_pt(sal_UV_timeEvol, potT_UV_timeEvol)


# Re-assign attributes (long_name, standard_name, units)
conT_UV_ctrl.attrs['long_name'] = 'conservative temperature'
conT_UV_ctrl.attrs['standard_name'] = 'conservative temperature (TEOS-10)'
conT_UV_ctrl.attrs['units'] = 'degC'
conT_UV_pmoc.attrs['long_name'] = 'conservative temperature'
conT_UV_pmoc.attrs['standard_name'] = 'conservative temperature (TEOS-10)'
conT_UV_pmoc.attrs['units'] = 'degC'
conT_UV_timeEvol.attrs['long_name'] = 'conservative temperature'
conT_UV_timeEvol.attrs['standard_name'] = 'conservative temperature (TEOS-10)'
conT_UV_timeEvol.attrs['units'] = 'degC'
# Remove 'valid_range' attribute (from sal, n/a for conT)
del conT_UV_ctrl.attrs['valid_range']
del conT_UV_pmoc.attrs['valid_range']
del conT_UV_timeEvol.attrs['valid_range']


# (2) Calc rho
# gsw: Calculates potential density anomaly with reference pressure of 0 dbar, 
#  this being this particular potential density minus 1000 kg/m^3.
rho_UV_ctrl = gsw.density.sigma0(sal_UV_ctrl, conT_UV_ctrl)
rho_UV_pmoc = gsw.density.sigma0(sal_UV_pmoc, conT_UV_pmoc)
rho_UV_timeEvol = gsw.density.sigma0(sal_UV_timeEvol, conT_UV_timeEvol)

# Add 1000 to make it not an anomaly
rho_UV_ctrl = rho_UV_ctrl + 1000
rho_UV_pmoc = rho_UV_pmoc + 1000
rho_UV_timeEvol = rho_UV_timeEvol + 1000

# Re-assign attributes (long_name, units)
rho_UV_ctrl.attrs['long_name'] = 'potential density 0dbar'
rho_UV_ctrl.attrs['standard_name'] = 'potential density 0dbar'
rho_UV_ctrl.attrs['units'] = 'kg m-3'
rho_UV_pmoc.attrs['long_name'] = 'potential density 0dbar'
rho_UV_pmoc.attrs['standard_name'] = 'potential density 0dbar'
rho_UV_pmoc.attrs['units'] = 'kg m-3'
rho_UV_timeEvol.attrs['long_name'] = 'potential density 0dbar'
rho_UV_timeEvol.attrs['standard_name'] = 'potential density 0dbar'
rho_UV_timeEvol.attrs['units'] = 'kg m-3'
# Remove 'valid_range' attribute (from sal, n/a for conT)
del rho_UV_ctrl.attrs['valid_range']
del rho_UV_pmoc.attrs['valid_range']
del rho_UV_timeEvol.attrs['valid_range']



# %%% Use rho_UV to convert alk, dic, po4 mol m-3 > µmol kg-1.

dic_UV_ctrl = dic_UV_ctrl / rho_UV_ctrl *(10**6)
dic_UV_pmoc = dic_UV_pmoc / rho_UV_pmoc *(10**6)
dic_UV_timeEvol = dic_UV_timeEvol / rho_UV_timeEvol *(10**6)

alk_UV_ctrl = alk_UV_ctrl / rho_UV_ctrl *(10**6)
alk_UV_pmoc = alk_UV_pmoc / rho_UV_pmoc *(10**6)
alk_UV_timeEvol = alk_UV_timeEvol / rho_UV_timeEvol *(10**6)

po4_UV_ctrl = po4_UV_ctrl / rho_UV_ctrl *(10**6)
po4_UV_pmoc = po4_UV_pmoc / rho_UV_pmoc *(10**6)
po4_UV_timeEvol = po4_UV_timeEvol / rho_UV_timeEvol *(10**6)


dic_UV_ctrl.attrs['units'] = 'µmol kg-1'
dic_UV_pmoc.attrs['units'] = 'µmol kg-1'
dic_UV_timeEvol.attrs['units'] = 'µmol kg-1'

alk_UV_ctrl.attrs['units'] = 'µmol kg-1'
alk_UV_pmoc.attrs['units'] = 'µmol kg-1'
alk_UV_timeEvol.attrs['units'] = 'µmol kg-1'

po4_UV_ctrl.attrs['units'] = 'µmol kg-1'
po4_UV_pmoc.attrs['units'] = 'µmol kg-1'
po4_UV_timeEvol.attrs['units'] = 'µmol kg-1'


# %%% Rename all dims to "depth", "lat", "lon"

dic_UV_ctrl = dic_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
dic_UV_pmoc = dic_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
dic_UV_timeEvol = dic_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})

alk_UV_ctrl = alk_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
alk_UV_pmoc = alk_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
alk_UV_timeEvol = alk_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})

po4_UV_ctrl = po4_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
po4_UV_pmoc = po4_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
po4_UV_timeEvol = po4_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})

potT_UV_ctrl = potT_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
potT_UV_pmoc = potT_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
potT_UV_timeEvol = potT_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})

sal_UV_ctrl = sal_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
sal_UV_pmoc = sal_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
sal_UV_timeEvol = sal_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})

rho_UV_ctrl = rho_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
rho_UV_pmoc = rho_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
rho_UV_timeEvol = rho_UV_timeEvol.rename({'longitude':'lon', 'latitude':'lat'})





# Probably no longer needed, probably for time-evolving DIC anomaly depth 
# transect  figure (Supp Fig 5)

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
x = xr.open_dataset(pathname+filename, decode_times=False) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# XXX - Take last timestep at the end of the 1000 years
dic_UV_pmoc_full  = x.O_dic# .isel(time=-1)    # time: 10

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg.02001.01.01.nc'
x = xr.open_dataset(pathname+filename, decode_times=False)
# time=1. array([2000.5]) (the end of the LGM simulation - our control)
dic_UV_ctrl_full  = x.O_dic.isel(time=0)  # only 1 time-pt anyway (2e+05)





# %% (2) Read in LOVECLIM Output

# %%% Read in LOVECLIM dic, alk, po4, T, S, rho

pathname = '../data/Output_2014_LOVECLIM_L-fNA/';
# T10 = 11. array([0., 200., 400., 600., 800., 1000., 1200., 1400., 1600., 1800., 2000.])
# first time step = initial starting/ctrl point. Then n=10 10yr avgs at end of every 200 years
dic_LC_ctrl = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['CT'].isel(T10=0)
alk_LC_ctrl = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['AT'].isel(T10=0)
po4_LC_ctrl = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['PO4'].isel(T10=0)

dic_LC_pmoc = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['CT'].sel(T10=1000.0)
alk_LC_pmoc = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['AT'].sel(T10=1000.0)
po4_LC_pmoc = xr.open_dataset(pathname+'LOVECLIM-BGCnoComp.nc')['PO4'].sel(T10=1000.0)


# temp
# time = 2000 (annual resolution). Avg first and last 10 yrs of simulation for your control and experiment. Units: Kelvin, convert later
potT_LC_ctrl = xr.open_dataset(pathname+'temp-19-172.nc')['temp'].isel(time=slice(0,10)).mean(dim='time')
potT_LC_pmoc = xr.open_dataset(pathname+'temp-19-172.nc')['temp'].isel(time=slice(990,1000)).mean(dim='time')

# sal
# time = 2000 (annual resolution). units: psu
sal_LC_ctrl = xr.open_dataset(pathname+'sal-19-172.nc')['sal'].isel(time=slice(0,10)).mean(dim='time')
sal_LC_pmoc = xr.open_dataset(pathname+'sal-19-172.nc')['sal'].isel(time=slice(990,1000)).mean(dim='time')

# time = 2000 (annual resolution). units: kg/m^3
rho_LC_ctrl = xr.open_dataset(pathname+'rho-19-172.nc')['rho'].isel(time=slice(0,10)).mean(dim='time')
rho_LC_pmoc = xr.open_dataset(pathname+'rho-19-172.nc')['rho'].isel(time=slice(990,1000)).mean(dim='time')



# %%% Rename all dims to "depth", "lat", "lon"

dic_LC_ctrl = dic_LC_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
dic_LC_pmoc = dic_LC_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

alk_LC_ctrl = alk_LC_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
alk_LC_pmoc = alk_LC_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

po4_LC_ctrl = po4_LC_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
po4_LC_pmoc = po4_LC_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

potT_LC_ctrl  = potT_LC_ctrl.rename({'z':'depth'})
potT_LC_pmoc  = potT_LC_pmoc.rename({'z':'depth'})

sal_LC_ctrl  = sal_LC_ctrl.rename({'z':'depth'})
sal_LC_pmoc  = sal_LC_pmoc.rename({'z':'depth'})

rho_LC_ctrl  = rho_LC_ctrl.rename({'z':'depth'})
rho_LC_pmoc  = rho_LC_pmoc.rename({'z':'depth'})


# %%% Fix BGC lon (27.5-384.5)

# LC output has lon that goes [27.5:384.5]. Fix these to wrap back around to 0
#  at 360, and re-order to go monotonically up [0:360]
# This fixes values >360: ds.coords[lon] = np.mod(LC14_ctrl.lon.values,360)
# This re-orders: ds.reindex({ lon : np.sort(ds[lon])})
# Source: https://github.com/pydata/xarray/issues/577

# Re-number

dic_LC_ctrl.coords['lon'] = np.mod(dic_LC_ctrl['lon'], 360)
dic_LC_pmoc.coords['lon'] = np.mod(dic_LC_pmoc['lon'], 360)

alk_LC_ctrl.coords['lon'] = np.mod(alk_LC_ctrl['lon'], 360)
alk_LC_pmoc.coords['lon'] = np.mod(alk_LC_pmoc['lon'], 360)

po4_LC_ctrl.coords['lon'] = np.mod(po4_LC_ctrl['lon'], 360)
po4_LC_pmoc.coords['lon'] = np.mod(po4_LC_pmoc['lon'], 360)


# Re-order

dic_LC_ctrl = dic_LC_ctrl.reindex({ 'lon': np.sort(dic_LC_ctrl['lon'])}) 
dic_LC_pmoc = dic_LC_pmoc.reindex({ 'lon': np.sort(dic_LC_pmoc['lon'])}) 

alk_LC_ctrl = alk_LC_ctrl.reindex({ 'lon': np.sort(alk_LC_ctrl['lon'])}) 
alk_LC_pmoc = alk_LC_pmoc.reindex({ 'lon': np.sort(alk_LC_pmoc['lon'])}) 

po4_LC_ctrl = po4_LC_ctrl.reindex({ 'lon': np.sort(po4_LC_ctrl['lon'])}) 
po4_LC_pmoc = po4_LC_pmoc.reindex({ 'lon': np.sort(po4_LC_pmoc['lon'])}) 



# %%% Fix PHYS depth (make go surf-down)

potT_LC_ctrl = potT_LC_ctrl.reindex(depth = potT_LC_ctrl.depth[::-1]) 
potT_LC_pmoc = potT_LC_pmoc.reindex(depth = potT_LC_pmoc.depth[::-1])  

sal_LC_ctrl = sal_LC_ctrl.reindex(depth = sal_LC_ctrl.depth[::-1]) 
sal_LC_pmoc = sal_LC_pmoc.reindex(depth = sal_LC_pmoc.depth[::-1])  

rho_LC_ctrl = rho_LC_ctrl.reindex(depth = rho_LC_ctrl.depth[::-1]) 
rho_LC_pmoc = rho_LC_pmoc.reindex(depth = rho_LC_pmoc.depth[::-1])  



# %%% Convert T from Kelvin to degC

# Note from Filestor: "data in degK, please make sure you define values below 271K as bad."

# Only return numbers >= 271 (replaces with NaNs everywhere else)
potT_LC_ctrl = potT_LC_ctrl.where(potT_LC_ctrl > 271)
potT_LC_pmoc = potT_LC_pmoc.where(potT_LC_pmoc > 271)

# Then convert units
potT_LC_ctrl = potT_LC_ctrl - 273.15
potT_LC_ctrl.attrs['units'] = 'degC'

potT_LC_pmoc = potT_LC_pmoc - 273.15
potT_LC_pmoc.attrs['units'] = 'degC'



# %%% Interp PHYS vars onto BGC coords

potT_LC_ctrl = potT_LC_ctrl.interp_like(dic_LC_ctrl)
potT_LC_pmoc = potT_LC_pmoc.interp_like(dic_LC_ctrl)

sal_LC_ctrl = sal_LC_ctrl.interp_like(dic_LC_ctrl)
sal_LC_pmoc = sal_LC_pmoc.interp_like(dic_LC_ctrl)

rho_LC_ctrl = rho_LC_ctrl.interp_like(dic_LC_ctrl)
rho_LC_pmoc = rho_LC_pmoc.interp_like(dic_LC_ctrl)



# %% (3) Read in LC-LGM Output

# %%% Read in LC-LGM dic, alk, po4, T, S, rho

# NO time dimension for bgc output. Just a  ctrl and a pmoc simulation, 1 data-time-point each

# V3L - "North Atlantic strong" - Ctrl.
pathname = '../data/Output_2016_LOVECLIM-LGM/V3L/';
dic_LGM_ctrl = xr.open_dataset(pathname+'dic.nc')['CT']
alk_LGM_ctrl = xr.open_dataset(pathname+'alk.nc')['AT']
po4_LGM_ctrl = xr.open_dataset(pathname+'V3L-PO4-18250.nc')['PO']

# cresta is 200yrs in annual res, yr.800-1000. Take average of last 10yr
potT_LGM_ctrl = xr.open_dataset(pathname+'cresta18250.nc')['temp'].sel(time=slice(191,200)).mean(dim='time')
sal_LGM_ctrl = xr.open_dataset(pathname+'cresta18250.nc')['sal'].sel(time=slice(191,200)).mean(dim='time')
rho_LGM_ctrl = xr.open_dataset(pathname+'cresta18250.nc')['rho'].sel(time=slice(191,200)).mean(dim='time')


# V3LNAw - "North Atlantic weak" - PMOC. 0.05Sv freshwater to N.Atl  
pathname = '../data/Output_2016_LOVECLIM-LGM/V3LNAw/';
dic_LGM_pmoc = xr.open_dataset(pathname+'dic.nc')['CT']
alk_LGM_pmoc = xr.open_dataset(pathname+'alk.nc')['AT']
po4_LGM_pmoc = xr.open_dataset(pathname+'V3LNAw-PO4.nc')['PO']

# cresta is 200yrs in annual res, yr.800-1000. Take average of last 10yr
potT_LGM_pmoc = xr.open_dataset(pathname+'cresta18250.nc')['temp'].sel(time=slice(191,200)).mean(dim='time')
sal_LGM_pmoc = xr.open_dataset(pathname+'cresta18250.nc')['sal'].sel(time=slice(191,200)).mean(dim='time')
rho_LGM_pmoc = xr.open_dataset(pathname+'cresta18250.nc')['rho'].sel(time=slice(191,200)).mean(dim='time')



# %%% Rename all dims to "depth", "lat", "lon"

dic_LGM_ctrl = dic_LGM_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
dic_LGM_pmoc = dic_LGM_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

alk_LGM_ctrl = alk_LGM_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
alk_LGM_pmoc = alk_LGM_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

po4_LGM_ctrl = po4_LGM_ctrl.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})
po4_LGM_pmoc = po4_LGM_pmoc.rename({'Z':'depth', 'Y10':'lat', 'X10':'lon'})

potT_LGM_ctrl  = potT_LGM_ctrl.rename({'z':'depth'})
potT_LGM_pmoc  = potT_LGM_pmoc.rename({'z':'depth'})

sal_LGM_ctrl  = sal_LGM_ctrl.rename({'z':'depth'})
sal_LGM_pmoc  = sal_LGM_pmoc.rename({'z':'depth'})

rho_LGM_ctrl  = rho_LGM_ctrl.rename({'z':'depth'})
rho_LGM_pmoc  = rho_LGM_pmoc.rename({'z':'depth'})


# %%% Fix BGC lon (27.5-384.5)

# LC output has lon that goes [27.5:384.5]. Fix these to wrap back around to 0
#  at 360, and re-order to go monotonically up [0:360]
# This fixes values >360: ds.coords[lon] = np.mod(LC14_ctrl.lon.values,360)
# This re-orders: ds.reindex({ lon : np.sort(ds[lon])})
# Source: https://github.com/pydata/xarray/issues/577

# Re-number

dic_LGM_ctrl.coords['lon'] = np.mod(dic_LGM_ctrl['lon'], 360)
dic_LGM_pmoc.coords['lon'] = np.mod(dic_LGM_pmoc['lon'], 360)

alk_LGM_ctrl.coords['lon'] = np.mod(alk_LGM_ctrl['lon'], 360)
alk_LGM_pmoc.coords['lon'] = np.mod(alk_LGM_pmoc['lon'], 360)

po4_LGM_ctrl.coords['lon'] = np.mod(po4_LGM_ctrl['lon'], 360)
po4_LGM_pmoc.coords['lon'] = np.mod(po4_LGM_pmoc['lon'], 360)



# Re-order

dic_LGM_ctrl = dic_LGM_ctrl.reindex({ 'lon': np.sort(dic_LGM_ctrl['lon'])}) 
dic_LGM_pmoc = dic_LGM_pmoc.reindex({ 'lon': np.sort(dic_LGM_pmoc['lon'])}) 

alk_LGM_ctrl = alk_LGM_ctrl.reindex({ 'lon': np.sort(alk_LGM_ctrl['lon'])}) 
alk_LGM_pmoc = alk_LGM_pmoc.reindex({ 'lon': np.sort(alk_LGM_pmoc['lon'])}) 

po4_LGM_ctrl = po4_LGM_ctrl.reindex({ 'lon': np.sort(po4_LGM_ctrl['lon'])}) 
po4_LGM_pmoc = po4_LGM_pmoc.reindex({ 'lon': np.sort(po4_LGM_pmoc['lon'])}) 



# %%% Fix PHYS depth (make go surf-down)

potT_LGM_ctrl = potT_LGM_ctrl.reindex(depth = potT_LGM_ctrl.depth[::-1]) 
potT_LGM_pmoc = potT_LGM_pmoc.reindex(depth = potT_LGM_pmoc.depth[::-1])  

sal_LGM_ctrl = sal_LGM_ctrl.reindex(depth = sal_LGM_ctrl.depth[::-1]) 
sal_LGM_pmoc = sal_LGM_pmoc.reindex(depth = sal_LGM_pmoc.depth[::-1])  

rho_LGM_ctrl = rho_LGM_ctrl.reindex(depth = rho_LGM_ctrl.depth[::-1]) 
rho_LGM_pmoc = rho_LGM_pmoc.reindex(depth = rho_LGM_pmoc.depth[::-1])  



# %%% Convert T from Kelvin to degC

# First replace 0s with NaNs - nah, new temp output from Laurie for LC she said 
# all numbers <271 are bad.

potT_LGM_ctrl = potT_LGM_ctrl.where(potT_LGM_ctrl > 271)
potT_LGM_pmoc = potT_LGM_pmoc.where(potT_LGM_pmoc > 271)


# Then convert units
potT_LGM_ctrl = potT_LGM_ctrl - 273.15
potT_LGM_ctrl.attrs['units'] = 'degC'

potT_LGM_pmoc = potT_LGM_pmoc - 273.15
potT_LGM_pmoc.attrs['units'] = 'degC'



# %%% Interp PHYS vars onto BGC coords

potT_LGM_ctrl = potT_LGM_ctrl.interp_like(dic_LGM_ctrl)
potT_LGM_pmoc = potT_LGM_pmoc.interp_like(dic_LGM_ctrl)

sal_LGM_ctrl = sal_LGM_ctrl.interp_like(dic_LGM_ctrl)
sal_LGM_pmoc = sal_LGM_pmoc.interp_like(dic_LGM_ctrl)

rho_LGM_ctrl = rho_LGM_ctrl.interp_like(dic_LGM_ctrl)
rho_LGM_pmoc = rho_LGM_pmoc.interp_like(dic_LGM_ctrl)



# %% (4) Read in GENIE Output

# %%% Read in in GENIE dic, alk, po4, T, S

# Get data from:
pathname = '../data/Output_2020_cGENIE_Rae/';

vars = ['ocn_DIC_Snorm', 'ocn_ALK_Snorm', 'ocn_PO4_Snorm', 'ocn_temp', 'ocn_sal']
# note the ctrl simulation has rho (phys_ocn_rho) but not the experiment simulation

rawdata_GENIE_ctrl = xr.open_dataset(pathname+'Spin_LGM__SPIN_worjh2_Fe14C_preAge_Dye_LGMsave99/'+'fields_biogem_3d.nc')[vars]
rawdata_GENIE_exp = xr.open_dataset(pathname+'PA-15_negpt28Sv/'+'fields_biogem_3d.nc')[vars]
# PA35_negpt12Sv WEAKEST pmoc
# PA15_negpt19Sv MOST LGM-like (Rae 2020 Sci)
# PA-15_negpt28Sv MOST extm

# Same for both ctl & exp
    # ocn_DIC_Snorm     DIC normalized by salinity [mol kg-1]
    # ocn_ALK_Snorm     ALK normalized by salinity [mol kg-1]  
    # ocn_PO4_Snorm     PO4 normalized by salinity [mol kg-1]   
    # ocn_temp          temperature [degrees C]
    # ocn_sal           salinity [PSU]
    #
    # (time: 1, zt: 16, lat: 36, lon: 36) (Note, time exp: 12 ( ))
    # * time     (time) float64 999.5 // or, exp, n=12: [0.5, 1.5, 4.5, 9.5, 19.5, 49.5, 99.5, 199.5, 499.5, 999.5, 1999.5, 4999.5]
    # * lon      (lon) float64 -255.0 -245.0 -235.0 -225.0 ... 65.0 75.0 85.0 95.0
    # * lat      (lat) float64 -76.46 -66.44 -59.44 -53.66 ... 59.44 66.44 76.46
    # * zt       (zt) float64 40.42 127.6 228.8 ... 3.283e+03 3.894e+03 4.604e+03


# - - - - - Rename variables to short names - - - - - #

rawdata_GENIE_ctrl = rawdata_GENIE_ctrl.rename_vars({"ocn_DIC_Snorm": "dic", "ocn_ALK_Snorm": "alk", \
                                         "ocn_PO4_Snorm": "po4", "ocn_temp"     : "T", \
                                         "ocn_sal"      : "sal"})
rawdata_GENIE_exp = rawdata_GENIE_exp.rename_vars({"ocn_DIC_Snorm": "dic", "ocn_ALK_Snorm": "alk", \
                                         "ocn_PO4_Snorm": "po4", "ocn_temp"     : "T", \
                                         "ocn_sal"      : "sal"})




# %%% Rename all dims to "depth", "lat", "lon"

rawdata_GENIE_ctrl = rawdata_GENIE_ctrl.rename({'zt':'depth'})
rawdata_GENIE_ctrl = rawdata_GENIE_ctrl.isel(time=-1)

rawdata_GENIE_exp = rawdata_GENIE_exp.rename({'zt':'depth'})
rawdata_GENIE_exp = rawdata_GENIE_exp.isel(time=-1)



# %%% Convert units on BGC

# Convert DIC, alk, po4 from [mol kg-1] to [umol kg-1] for PyCO2SYS
# temp in degC (good) and sal in psu (good)

xr.set_options(keep_attrs=True)
convert_vars = ['dic', 'alk', 'po4']

for var_str in convert_vars:
    # Ctrl
    rawdata_GENIE_ctrl[var_str] = rawdata_GENIE_ctrl[var_str]*(10**6)
    rawdata_GENIE_ctrl[var_str] = rawdata_GENIE_ctrl[var_str].assign_attrs(units='µmol kg-1')
    # Exp
    rawdata_GENIE_exp[var_str] = rawdata_GENIE_exp[var_str]*(10**6)
    rawdata_GENIE_exp[var_str] = rawdata_GENIE_exp[var_str].assign_attrs(units='µmol kg-1')




# %%% Calc GENIE rho



# (1) Convert depth to sea pressure ('p' [dbar], absolute pressure - 10.1325 dbar)
# depth (z) is positive up [m], so need to *-1

# Make a 3D arrays for depth, lat
depth_zlatlon = rawdata_GENIE_ctrl.alk.depth.expand_dims(dim={"lat": np.size(rawdata_GENIE_ctrl.alk.lat), "lon": np.size(rawdata_GENIE_ctrl.alk.lon)})
depth_zlatlon = depth_zlatlon.assign_coords(lat=rawdata_GENIE_ctrl.alk.lat, lon=rawdata_GENIE_ctrl.alk.lon)
depth_zlatlon = depth_zlatlon*-1

lat_zlatlon = rawdata_GENIE_ctrl.alk.lat.expand_dims(dim={"depth": np.size(rawdata_GENIE_ctrl.alk.depth), "lon": np.size(rawdata_GENIE_ctrl.alk.lon)})
lat_zlatlon = lat_zlatlon.assign_coords(depth=rawdata_GENIE_ctrl.alk.depth, lon=rawdata_GENIE_ctrl.alk.lon)

lon_zlatlon = rawdata_GENIE_ctrl.alk.lon.expand_dims(dim={"depth": np.size(rawdata_GENIE_ctrl.alk.depth), "lat": np.size(rawdata_GENIE_ctrl.alk.lat)})
lon_zlatlon = lon_zlatlon.assign_coords(depth=rawdata_GENIE_ctrl.alk.depth, lat=rawdata_GENIE_ctrl.alk.lat)


abs_press = gsw_c.p_from_z(depth_zlatlon, lat_zlatlon)
abs_press = abs_press.transpose('depth', 'lat', 'lon') # re-order dimensions
# Re-assign units and long_name, and add note, delete other atts
abs_press = abs_press.assign_attrs(units='dbar', long_name='pressure', note='sea pressure ( i.e. absolute pressure - 10.1325 dbar )')
del abs_press.attrs['axis']; del abs_press.attrs['edges']; del abs_press.attrs['standard_name']


# (2) Convert Practical Salinity to Absolute Salinity [g kg-1]

abssal_exp = gsw_c.SA_from_SP(rawdata_GENIE_exp.sal, abs_press, lon_zlatlon, lat_zlatlon)
abssal_exp = abssal_exp.assign_attrs(units='g kg-1', long_name='absolute salinity')

abssal_ctrl = gsw_c.SA_from_SP(rawdata_GENIE_ctrl.sal, abs_press, lon_zlatlon, lat_zlatlon)
abssal_ctrl = abssal_ctrl.assign_attrs(units='g kg-1', long_name='absolute salinity')


# (3) Get Conservative Temperature from potential T

contT_exp_potT  = gsw_c.CT_from_pt(abssal_exp, rawdata_GENIE_exp.T)
contT_exp_potT = contT_exp_potT.assign_attrs(units='degrees C', long_name='conservative temperature', valid_range='[-9.999 99.999]')

contT_ctrl_potT  = gsw_c.CT_from_pt(abssal_ctrl, rawdata_GENIE_ctrl.T)
contT_ctrl_potT = contT_ctrl_potT.assign_attrs(units='degrees C', long_name='conservative temperature', valid_range='[-9.999 99.999]')


# (4) sigma0 from CT_from_insitu and CT_from_potT

rho_GENIE_ctrl = gsw_d.sigma0(abssal_ctrl, contT_ctrl_potT)
rho_GENIE_exp  = gsw_d.sigma0(abssal_exp, contT_exp_potT)

rho_GENIE_ctrl  = rho_GENIE_ctrl.assign_attrs(units='kg m-3', long_name='sigma 0dbar', note='potential density')
del rho_GENIE_ctrl.attrs['valid_range']    
rho_GENIE_exp  = rho_GENIE_exp.assign_attrs(units='kg m-3', long_name='sigma 0dbar', note='potential density')
del rho_GENIE_exp.attrs['valid_range']    
    

rawdata_GENIE_ctrl['rho'] = rho_GENIE_ctrl + 1000 # to make not anom
rawdata_GENIE_exp['rho'] = rho_GENIE_ctrl + 1000 # to make not anom



# %% (5) Calc PCO2 on all


# %%% Info on PyCO2SYS parameters

# Following Chen et al 2022 (following Wiliams et al 2017)
# Potential pCO2 (PCO2): the pCO2 a water parcel would have if brought 
#  adiabatically up to the surface. Corrects for P effects on T and partial P.
# (Chen et al used CO2SYS (van Heuven et al 2011)

# Following UVic's protocol, we will use the...
# dissociation constants of carbonic acid    Dickson & Millero (1987) DSR Part A  
#                   "  " of bisulfate        Dickson (1990) J. Chem. Thermodyn.
#                   "  " of hyd. fluoride    Dickson & Riley (1979) Mar Chem 7  
#  and the boron to salinity ratio           Lee et al (2010) GCA                


# PyCO2SYS v1.8.1 (Humphreys et al., 2022)
# https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
# INPUT
#   - Two carbonate system parameters: type=1 tot ALK in umolkg, type2 DIC in umolkg
#   - opt_pH_scale     = 2   seawater scale
#   - opt_k_carbonic   = 5   Dickson & Millero '87 (H73a, H73b and MCHP73 refit by DM87)
#   - opt_k_bisulfate  = 1   Dickson 1990 (D90a)
#   - opt_k_fluoride   = 1   Dickson & Riley 1979 (DR79)
#   - opt_total_borate = 2   Lee et al 2010 (LKB10) (this is the boron:salinity ratio to estimate total borate)
# Have model output po4, include it:
#   - total_phosphate: total phosphate in μmol·kg−1
# OUTPUT
#   - Results: 'pCO2': seawater partial pressure of CO2 (uatm)

# input_P = 0.0     # dbar
# # input sal and T will be in-situ



# %%% Define a func for calc'ing PCO2

def construct_3D_PCO2_wPO4(alk, dic, sal, T, po4, west_lon, east_lon, P=0):
    """
    Feed in 5 3D xr.DataSets (alk, dic, and po4, all in umol kg-1; plus potT in
    degC and practical salinity in psu) (all with the same coordinates and dims 
    'depth', 'lat', 'lon'). The function  will return a 3D xr.DataArray of 
    PCO2 - calc'd with in-situ sal & T at the specified reference pressure P 
    (dbar) - for all depths and lats of the input datasets and across the lons 
    specified by west_lon and east_lon limits (degE). Note, PyCO2SYS parameters 
    set within the function.
    """
    
    # Snip down to lons of interest
    subset_alk = alk.sel(lon=slice(west_lon, east_lon))
    subset_dic = dic.sel(lon=slice(west_lon, east_lon))
    subset_sal = sal.sel(lon=slice(west_lon, east_lon))
    subset_temp =  T.sel(lon=slice(west_lon, east_lon))
    subset_po4 = po4.sel(lon=slice(west_lon, east_lon))

    # Set PyCO2SYS parameters
    input_P = P
    my_opt_pH_scale     = 2  # 2 = sw
    my_opt_k_carbonic   = 5  # 5 = H73a, H73b and MCHP73 refit by DM87
    my_opt_k_bisulfate  = 1  # 1 = D90a
    my_opt_k_fluoride   = 1  # 1 = DR79
    my_opt_total_borate = 2  # 2 = LKB10
    
    my_opt_buffers_mode = 2 # how to calculate the various buffer factors (or not)
    # 1: using automatic differentiation, which accounts for the effects of all equilibrating solutes (default).
    # 2: using explicit equations reported in the literature, which only account for carbonate, borate and water alkalinity.
    # 0: not at all.


    # Initialize PCO2_3Darray before the loop
    PCO2_3Darray = None

    # Iterate through lons of the dataset, calc'ing 2D (depth vs lat) PCO2
    # Save that 1st 2D (numpy) array, then append others onto it
    counter = 1
    
    for ii in subset_alk.lon.values:
        # Select a 2D slice along a certain line of lon
        input_alk  = subset_alk.sel(lon=ii)
        input_dic  = subset_dic.sel(lon=ii)
        input_sal  = subset_sal.sel(lon=ii)
        input_temp = subset_temp.sel(lon=ii)
        input_po4  = subset_po4.sel(lon=ii)
        # Calculate a 2D array of PCO2 from 'alk' and 'dic'
        temp_pyco2sys = pyco2.sys(par1 = input_alk, par2 = input_dic, par1_type=1, par2_type=2, \
                                  salinity = input_sal, temperature = input_temp, pressure = input_P, \
                                      opt_pH_scale = my_opt_pH_scale, opt_k_carbonic = my_opt_k_carbonic, \
                                          opt_k_bisulfate = my_opt_k_bisulfate, opt_k_fluoride = my_opt_k_fluoride, \
                                              opt_total_borate = my_opt_total_borate, opt_buffers_mode=my_opt_buffers_mode, \
                                                    total_phosphate = input_po4);
        PCO2 = temp_pyco2sys['pCO2']
        # print(PCO2)

        if counter == 1:
            # If it's the 1st 2D array of PCO2, save it as first_array
            first_array = PCO2
        elif counter == 2:
            # second time, stack the two to create a 3rd dimension
            PCO2_3Darray = np.stack([first_array,PCO2], axis=2)
        else:
            # after that, just append to end, extending out 3rd dim
            PCO2_3Darray = np.dstack((PCO2_3Darray, PCO2))
        counter += 1
        print(counter)
        
                
    # Now that that loop's done, should have a 3D array of PCO2 values that is 
    #  the same dims as 'subset_data'.
    # Create an xr.DataArray from the 3D PCO2 you have now, give it dims, and 
    #  assign coords from 'subset_data'
    da = xr.DataArray(
        data=PCO2_3Darray,
        dims=['depth', 'lat', 'lon'],
        attrs=dict(
            description="PCO_2",
            units="µatm",
        ),
    )
    
    # This should work, because if you assign "dimnames" in wrong order above (i.e.
    # in not the same order as x), the below won't work and will raise an 
    # "conflicting sizes" error
    return_array = da.assign_coords(depth=subset_alk.depth, lat=subset_alk.lat, lon=subset_alk.lon)
                 
    return return_array  




# %%% Calc PCO2 across the globe

my_west_lon = 0    # 140-280 for Pac only, 0-360 for global
my_east_lon = 360  



# With PO4

PCO2_UV_ctrl = construct_3D_PCO2_wPO4(alk_UV_ctrl, dic_UV_ctrl, sal_UV_ctrl, potT_UV_ctrl, po4_UV_ctrl, \
                                 my_west_lon, my_east_lon, P=0)
PCO2_UV_pmoc = construct_3D_PCO2_wPO4(alk_UV_pmoc, dic_UV_pmoc, sal_UV_pmoc, potT_UV_pmoc, po4_UV_pmoc, \
                                 my_west_lon, my_east_lon, P=0)

PCO2_LC_ctrl = construct_3D_PCO2_wPO4(alk_LC_ctrl, dic_LC_ctrl, sal_LC_ctrl, potT_LC_ctrl, po4_LC_ctrl, \
                                  my_west_lon, my_east_lon, P=0)
PCO2_LC_pmoc = construct_3D_PCO2_wPO4(alk_LC_pmoc, dic_LC_pmoc, sal_LC_pmoc, potT_LC_pmoc, po4_LC_pmoc, \
                                  my_west_lon, my_east_lon, P=0)    
    
PCO2_LGM_ctrl = construct_3D_PCO2_wPO4(alk_LGM_ctrl, dic_LGM_ctrl, sal_LGM_ctrl, potT_LGM_ctrl, po4_LGM_ctrl, \
                                 my_west_lon, my_east_lon, P=0)
PCO2_LGM_pmoc = construct_3D_PCO2_wPO4(alk_LGM_pmoc, dic_LGM_pmoc, sal_LGM_pmoc, potT_LGM_pmoc, po4_LGM_pmoc, \
                                 my_west_lon, my_east_lon, P=0)

    
my_west_lon = -260    # 140-280 for Pac only, 0-360 for global (Genie different)
my_east_lon = 100  

dataset = rawdata_GENIE_ctrl
PCO2_GENIE_ctrl = construct_3D_PCO2_wPO4(dataset.alk, dataset.dic, dataset.sal, dataset.T, dataset.po4, \
                                  my_west_lon, my_east_lon, P=0)
dataset = rawdata_GENIE_exp
PCO2_GENIE_exp = construct_3D_PCO2_wPO4(dataset.alk, dataset.dic, dataset.sal, dataset.T, dataset.po4, \
                                  my_west_lon, my_east_lon, P=0)
    

      
# %% Calc time-evolving PCO2 over UVic simulation and save

my_west_lon = 0    # 140-280 for Pac only, 0-360 for global
my_east_lon = 360


for idx, ii in enumerate(rho_UV_timeEvol.time):

    alk_2use = alk_UV_timeEvol.sel(time=ii)
    dic_2use = dic_UV_timeEvol.sel(time=ii)
    sal_2use = sal_UV_timeEvol.sel(time=ii)
    potT_2use = potT_UV_timeEvol.sel(time=ii)
    po4_2use = po4_UV_timeEvol.sel(time=ii)

    # Either make new dataarray or concatenate onto old
    if idx == 0:
        PCO2_UV_timeEvol = construct_3D_PCO2_wPO4(alk_2use, dic_2use, sal_2use, potT_2use, po4_2use, \
                                         my_west_lon, my_east_lon, P=0)
        # And expand "time" to make it the 0th dimension
        PCO2_UV_timeEvol.expand_dims(dim='time', axis=0)
    # Then, from now on calc PCO2, expand dim, and concatenate onto the above
    else:
        PCO2 = construct_3D_PCO2_wPO4(alk_2use, dic_2use, sal_2use, potT_2use, po4_2use, \
                                         my_west_lon, my_east_lon, P=0)
        PCO2 = PCO2.expand_dims(dim='time', axis=0)
        PCO2_UV_timeEvol = xr.concat([PCO2_UV_timeEvol, PCO2], dim='time')

    print('Loop idx = '+str(idx)+' done')


# Save time-evolving PCO2
PCO2_UV_timeEvol.to_netcdf('../results/PCO2_UV_timeEvol.nc')




# %% (6) Compile vars into DataSets

UV_ctrl = alk_UV_ctrl.to_dataset(name='alk')
UV_ctrl['dic'] = dic_UV_ctrl 
UV_ctrl['po4'] = po4_UV_ctrl 
UV_ctrl['potT'] = potT_UV_ctrl 
UV_ctrl['sal']  = sal_UV_ctrl 
UV_ctrl['rho']  = rho_UV_ctrl 
UV_ctrl['PCO2'] = PCO2_UV_ctrl 

UV_pmoc = alk_UV_pmoc.to_dataset(name='alk')
UV_pmoc['dic'] = dic_UV_pmoc 
UV_pmoc['po4'] = po4_UV_pmoc
UV_pmoc['potT'] = potT_UV_pmoc 
UV_pmoc['sal']  = sal_UV_pmoc 
UV_pmoc['rho']  = rho_UV_pmoc 
UV_pmoc['PCO2'] = PCO2_UV_pmoc 



# LC

LC_ctrl = alk_LC_ctrl.to_dataset(name='alk')
LC_ctrl['dic'] = dic_LC_ctrl 
LC_ctrl['po4'] = po4_LC_ctrl
LC_ctrl['potT'] = potT_LC_ctrl
LC_ctrl['sal']  = sal_LC_ctrl
LC_ctrl['rho']  = rho_LC_ctrl 
LC_ctrl['PCO2'] = PCO2_LC_ctrl 

LC_pmoc = alk_LC_pmoc.to_dataset(name='alk')
LC_pmoc['dic'] = dic_LC_pmoc 
LC_pmoc['po4'] = po4_LC_pmoc
LC_pmoc['potT'] = potT_LC_pmoc
LC_pmoc['sal']  = sal_LC_pmoc
LC_pmoc['rho']  = rho_LC_pmoc 
LC_pmoc['PCO2'] = PCO2_LC_pmoc 



# LGM

LGM_ctrl = alk_LGM_ctrl.to_dataset(name='alk')
LGM_ctrl['dic'] = dic_LGM_ctrl 
LGM_ctrl['po4'] = po4_LGM_ctrl
LGM_ctrl['potT'] = potT_LGM_ctrl 
LGM_ctrl['sal']  = sal_LGM_ctrl 
LGM_ctrl['rho']  = rho_LGM_ctrl 
LGM_ctrl['PCO2'] = PCO2_LGM_ctrl 

LGM_pmoc = alk_LGM_pmoc.to_dataset(name='alk')
LGM_pmoc['dic'] = dic_LGM_pmoc 
LGM_pmoc['po4'] = po4_LGM_pmoc
LGM_pmoc['potT'] = potT_LGM_pmoc 
LGM_pmoc['sal']  = sal_LGM_pmoc 
LGM_pmoc['rho']  = rho_LGM_pmoc 
LGM_pmoc['PCO2'] = PCO2_LGM_pmoc 
 


# GENIE

GENIE_ctrl = rawdata_GENIE_ctrl.alk.to_dataset(name='alk')
GENIE_ctrl['dic'] = rawdata_GENIE_ctrl.dic 
GENIE_ctrl['po4'] = rawdata_GENIE_ctrl.po4
GENIE_ctrl['potT'] = rawdata_GENIE_ctrl.T
GENIE_ctrl['sal']  = rawdata_GENIE_ctrl.sal
GENIE_ctrl['rho']  = rawdata_GENIE_ctrl.rho 
GENIE_ctrl['PCO2'] = PCO2_GENIE_ctrl 

GENIE_pmoc = rawdata_GENIE_exp.alk.to_dataset(name='alk')
GENIE_pmoc['dic'] = rawdata_GENIE_exp.dic 
GENIE_pmoc['po4'] = rawdata_GENIE_exp.po4
GENIE_pmoc['potT'] = rawdata_GENIE_exp.T
GENIE_pmoc['sal']  = rawdata_GENIE_exp.sal
GENIE_pmoc['rho']  = rawdata_GENIE_exp.rho 
GENIE_pmoc['PCO2'] = PCO2_GENIE_exp 



# %%% Fix all GENIE lon values: make go 0:360 instead of -255:95
# Orig:
# array([-255., -245., -235., -225., -215., -205., -195., -185., -175.,
#        -165., -155., -145., -135., -125., -115., -105.,  -95.,  -85.,
#         -75.,  -65.,  -55.,  -45.,  -35.,  -25.,  -15.,   -5.,    5.,
#          15.,   25.,   35.,   45.,   55.,   65.,   75.,   85.,   95.])
# Now (ok out of order):
# array([105., 115., 125., 135., 145., 155., 165., 175., 185., 195., 205.,
#        215., 225., 235., 245., 255., 265., 275., 285., 295., 305., 315.,
#        325., 335., 345., 355.,   5.,  15.,  25.,  35.,  45.,  55.,  65.,
#         75.,  85.,  95.])
GENIE_ctrl.coords['lon'] = (360 + GENIE_ctrl.lon.values) % 360
GENIE_pmoc.coords['lon'] = (360 + GENIE_pmoc.lon.values) % 360

GENIE_ctrl = GENIE_ctrl.sortby(GENIE_ctrl.lon)
GENIE_pmoc = GENIE_pmoc.sortby(GENIE_pmoc.lon)






# %% SUPP FIG: time-evolving PCO2 and anomaly, UVic

from matplotlib import pyplot as plt

# for ii in PCO2_UV_timeEvol.time:
#     fig, ax = plt.subplots()
#     PCO2_UV_timeEvol.sel(time=ii, lat=slice(-80,60)).interp(lon=200).plot(vmin=175, vmax=400)
#     ax.invert_yaxis()


# for ii in PCO2_UV_timeEvol.time:
#     fig, ax = plt.subplots()
#     var = PCO2_UV_timeEvol.sel(time=ii) - UV_ctrl.PCO2
#     var.sel(lat=slice(-80,60)).interp(lon=200).plot(vmin=-200, vmax=200, cmap='seismic')
#     ax.invert_yaxis()


lon_of_interest = 200


fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10.5,15))

# Left Col: absolute PCO2 evolving over ventilated UVic simulation 
col_idx = 0
var2plot = PCO2_UV_timeEvol
my_vmin = 170; my_vmax= 400
panel_labels = ['a', 'b', 'c', 'd', 'e']

for ii, time_idx in enumerate([0,3,5,7,9]):
    ax=axs[ii,col_idx]    
    var2plot.isel(time=time_idx).sel(lat=slice(-80,60)).interp(lon=lon_of_interest).plot(ax=ax, vmin=my_vmin, vmax=my_vmax) #, cbar_kwargs={'label': dic_diff.units})
    ax.invert_yaxis(); ax.set_facecolor([0.7, 0.7, 0.7])
    
    ax.set_title('%d years' % np.round(var2plot.time[time_idx].values, -2), fontweight='bold')
    ax.set_title(panel_labels[ii], loc='left', fontweight='bold')
    
    if ii != 4:
        ax.set_xlabel('')
    else:
        ax.set_xlabel(r'Latitude [$\degree$N] at '+str(360-lon_of_interest)+r'$\degree$W')

# Right Col: absolute PCO2 evolving over ventilated UVic simulation 
col_idx = 1
my_vmin = -200; my_vmax= 200
panel_labels = ['f', 'g', 'h', 'i', 'j']

for ii, time_idx in enumerate([0,3,5,7,9]):
    ax=axs[ii,col_idx]    
    var2plot = PCO2_UV_timeEvol.isel(time=time_idx) - UV_ctrl.PCO2
    
    var2plot.sel(lat=slice(-80,60)).interp(lon=lon_of_interest).plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap='seismic')
    ax.invert_yaxis(); ax.set_facecolor([0.7, 0.7, 0.7])
    
    ax.set_title('%d years' % np.round(PCO2_UV_timeEvol.time[time_idx].values, -2), fontweight='bold')
    ax.set_title(panel_labels[ii], loc='left', fontweight='bold')
    
    if ii != 4:
        ax.set_xlabel('')
    else:
        ax.set_xlabel(r'Latitude [$\degree$N] at '+str(360-lon_of_interest)+r'$\degree$W')



fig.tight_layout()
fig.suptitle(r'PCO$_2$ and PCO$_2$ Anomaly Over Time [µatm]', fontweight='bold', y=1.01, fontsize=15)


# Save fig
fig.savefig('../results/SuppFig_TimeEvolPCO2andAnom.eps', format='eps', dpi=600, bbox_inches='tight')
fig.savefig('../results/SuppFig_TimeEvolPCO2andAnom.png', dpi=600, bbox_inches='tight')
fig.savefig('../results/SuppFig_TimeEvolPCO2andAnom.pdf', dpi=600, bbox_inches='tight')

plt.show()




# %% (7) Calc sigma1/2 on all    

def calc_sig1_sig2(sal, potT):
    '''
    Feed in 3D arrays of same dims and coords of sal (g/kg) and pot. temperature (degC)

    '''
    
    # https://teos-10.github.io/GSW-Python/index.html
    # gsw = Python implementation of the Thermodynamic Equation of Seawater 2010 (TEOS-10) 
    
    # Will need abs salinity (good) & conservative temp (from my pot T) to calc
    # sigma0 (Pot dens ref 0dbar)
    
    import gsw

    # (1) Convert Pot Temp to Conservative Temp 
    # (takes abs sal (g/kg) and pot T "referenced to a sea pressure" (?) (C))
    
    conT = gsw.conversions.CT_from_pt(sal, potT)


    # (2) Calc sigmas
    
    # # gsw: Calculates potential density anomaly with reference pressure of 0 dbar, 
    # #  this being this particular potential density minus 1000 kg/m^3. Has inputs 
    # #  of abs sal and consv T.

    # Calculates potential density anomaly with reference pressure of 1000 dbar, 
    #  this being this particular potential density minus 1000 kg/m^3. Has inputs 
    #  of abs sal and consv T. (as above)
    sigma1 = gsw.density.sigma1(sal, conT)

    # And 2000dbar
    sigma2 = gsw.density.sigma2(sal, conT)

    # Re-assign attributes (long_name, units)
    sigma1.attrs['long_name'] = 'sigma1 - potential density anomaly 1000dbar'
    sigma1.attrs['standard_name'] = 'sigma1'
    sigma1.attrs['units'] = 'kg m-3'
    sigma2.attrs['long_name'] = 'sigma2 - potential density anomaly 2000dbar'
    sigma2.attrs['standard_name'] = 'sigma2'
    sigma2.attrs['units'] = 'kg m-3'

    return sigma1, sigma2
    
    
    

    
sigma1_UV_ctrl, sigma2_UV_ctrl = calc_sig1_sig2(UV_ctrl.sal, UV_ctrl.potT) 
sigma1_UV_pmoc, sigma2_UV_pmoc = calc_sig1_sig2(UV_pmoc.sal, UV_pmoc.potT) 

sigma1_LC_ctrl, sigma2_LC_ctrl = calc_sig1_sig2(LC_ctrl.sal, LC_ctrl.potT) 
sigma1_LC_pmoc, sigma2_LC_pmoc = calc_sig1_sig2(LC_pmoc.sal, LC_pmoc.potT) 

sigma1_LGM_ctrl, sigma2_LGM_ctrl = calc_sig1_sig2(LGM_ctrl.sal, LGM_ctrl.potT) 
sigma1_LGM_pmoc, sigma2_LGM_pmoc = calc_sig1_sig2(LGM_pmoc.sal, LGM_pmoc.potT) 

sigma1_GENIE_ctrl, sigma2_GENIE_ctrl = calc_sig1_sig2(GENIE_ctrl.sal, GENIE_ctrl.potT) 
sigma1_GENIE_pmoc, sigma2_GENIE_pmoc = calc_sig1_sig2(GENIE_pmoc.sal, GENIE_pmoc.potT) 



# %% (8) FIG3 and SUPP FIG2 : 4-paneled PCO2 and PO4 depth-lat section at 160W

# - - - - - Calc anoms - - - - - - - - - - - - - - - #

PCO2_anom_UV = UV_pmoc.PCO2 - UV_ctrl.PCO2
PCO2_anom_LC = LC_pmoc.PCO2 - LC_ctrl.PCO2
PCO2_anom_LGM   = LGM_pmoc.PCO2 - LGM_ctrl.PCO2
PCO2_anom_GENIE = GENIE_pmoc.PCO2 - GENIE_ctrl.PCO2


po4_anom_UV = UV_pmoc.po4 - UV_ctrl.po4
po4_anom_LC = LC_pmoc.po4 - LC_ctrl.po4
po4_anom_LGM   = LGM_pmoc.po4 - LGM_ctrl.po4
po4_anom_GENIE = GENIE_pmoc.po4 - GENIE_ctrl.po4


# - - - - -  Re-name "description" - - - - - - - - - #
PCO2_anom_UV.attrs['description'] = 'PCO_2 anomaly (perturbed - ctrl)'
PCO2_anom_LC.attrs['description'] = 'PCO_2 anomaly (perturbed - ctrl)'
PCO2_anom_LGM.attrs['description'] = 'PCO_2 anomaly (perturbed - ctrl)'
PCO2_anom_GENIE.attrs['description'] = 'PCO_2 anomaly (perturbed - ctrl)'

PCO2_anom_UV.attrs['model'] = 'UVic'
PCO2_anom_LC.attrs['model'] = 'LOVECLIM'
PCO2_anom_LGM.attrs['model'] = 'LOVECLIM-LGM'
PCO2_anom_GENIE.attrs['model'] = 'c-GENIE'


po4_anom_UV.attrs['description'] = 'po4 anomaly (perturbed - ctrl)'
po4_anom_LC.attrs['description'] = 'po4 anomaly (perturbed - ctrl)'
po4_anom_LGM.attrs['description'] = 'po4 anomaly (perturbed - ctrl)'
po4_anom_GENIE.attrs['description'] = 'po4 anomaly (perturbed - ctrl)'

po4_anom_UV.attrs['model'] = 'UVic'
po4_anom_LC.attrs['model'] = 'LOVECLIM'
po4_anom_LGM.attrs['model'] = 'LOVECLIM-LGM'
po4_anom_GENIE.attrs['model'] = 'c-GENIE'


# %% Make Fig3

lon_of_interest = 200


from matplotlib import pyplot as plt

my_fig_width  = 8.2  # 10.2  #11
my_fig_height = 5.2  # 6.5   # 7

cont_LW = 1.5  # contour line width and font size
cont_FS = 11
cont_accent='goldenrod'

PCO2_vars = [PCO2_anom_UV, PCO2_anom_GENIE, PCO2_anom_LC, PCO2_anom_LGM]
po4_vars = [po4_anom_UV, po4_anom_GENIE, po4_anom_LC, po4_anom_LGM]
sigma1_vars = [sigma1_UV_pmoc, sigma1_GENIE_pmoc, sigma1_LC_pmoc, sigma1_LGM_pmoc]
sigma2_vars = [sigma2_UV_pmoc, sigma2_GENIE_pmoc, sigma2_LC_pmoc, sigma2_LGM_pmoc]
rho_vars = [UV_pmoc.rho, GENIE_pmoc.rho, LC_pmoc.rho, LGM_pmoc.rho]
panel_labels = ['a', 'b', 'c', 'd']



# PCO2 Aesthetics
my_vmin = -200; my_vmax = -my_vmin
my_cmap = 'seismic'
my_cb_label = 'PCO$_2$ anomaly [µatm]'
vars2plot = PCO2_vars


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(my_fig_width, my_fig_height))

# Create a new axis for colorbar
cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Adjust the position and size as needed


# Loop through models (Note, (row*2 + col) goes [0, 1, 2, 3])
for (row, col), ax in np.ndenumerate(axs):  
    
    # # Plot PCO2 or po4
    var = vars2plot[row*2 + col]
    a = var.interp(lon=lon_of_interest).plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, add_colorbar = False)


    # # Plot isopycnals (all slightly different for each model)
    # rho_var = rho_vars[row*2 + col] - 1000
    
    if (row*2 + col) == 0:    # UVic
        red_levels = [36.0, 36.8]; yellow_levels = [36.5]
        
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        

        
    elif (row*2 + col) == 1:  # c-GENIE
        red_levels = [37, 38]; yellow_levels = [37.6]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
               
        
    elif (row*2 + col) == 2:  # LC
        red_levels = [37, 37.8]; yellow_levels = [37.6]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='r', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
                            
        
    elif (row*2 + col) == 3:  # LC-LGM
        red_levels = [33, 34.5]; yellow_levels = [34.0]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
            



    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_xlim((-80,60)); ax.set_ylim((0,5000)); ax.invert_yaxis()
    ax.set_title(var.model, fontweight='bold')
    if ((row*2 + col) == 0) or ((row*2 + col) == 1):
        ax.set_xlabel('')
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r'Latitude [$\degree$N] at '+str(360-lon_of_interest)+r'$\degree$W')
    if ((row*2 + col) == 1) or ((row*2 + col) == 3):
        ax.set_ylabel('')
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel(r'Depth [m]')

    # Add a, b, c, d annotations
    if ((row*2 + col) == 0) or ((row*2 + col) == 2):
        left_just = -0.07
    else:
        left_just = -0.025
    ax.text(left_just, 1.1, panel_labels[row*2 + col], transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')





# Add common colorbar
cbar = plt.colorbar(a, cax=cbar_ax)
cbar.set_label(label=my_cb_label, fontsize=14, weight='bold', rotation=-90, va='bottom', labelpad=1)
cbar.ax.tick_params(labelsize=12)

# Adjust so panels don't overlap with colorbar
bottom, top = 0.1, 0.9
left, right = 0.1, 0.9
fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.25, wspace=0.1)


# Save fig
plt.savefig('../results/Fig3_4panelPCO2sigma2.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig3_4panelPCO2sigma2.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig3_4panelPCO2sigma2.pdf', dpi=600, bbox_inches='tight')

plt.show()


# %% Make SuppFig3

lon_of_interest = 200


from matplotlib import pyplot as plt

my_fig_width  = 8.2  # 10.2  #11
my_fig_height = 5.2  # 6.5   # 7

cont_LW = 1.5  # contour line width and font size
cont_FS = 11
cont_accent = 'goldenrod'

PCO2_vars = [PCO2_anom_UV, PCO2_anom_GENIE, PCO2_anom_LC, PCO2_anom_LGM]
po4_vars = [po4_anom_UV, po4_anom_GENIE, po4_anom_LC, po4_anom_LGM]
sigma1_vars = [sigma1_UV_pmoc, sigma1_GENIE_pmoc, sigma1_LC_pmoc, sigma1_LGM_pmoc]
sigma2_vars = [sigma2_UV_pmoc, sigma2_GENIE_pmoc, sigma2_LC_pmoc, sigma2_LGM_pmoc]
rho_vars = [UV_pmoc.rho, GENIE_pmoc.rho, LC_pmoc.rho, LGM_pmoc.rho]
panel_labels = ['a', 'b', 'c', 'd']


# PO4 Aesthetics
my_vmin = -1.5; my_vmax = -my_vmin
my_cmap = 'seismic'
my_cb_label = 'PO$_4^{3-}$ anomaly [µmol kg$^{-1}$]'
vars2plot = po4_vars

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(my_fig_width, my_fig_height))

# Create a new axis for colorbar
cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Adjust the position and size as needed


# Loop through models (Note, (row*2 + col) goes [0, 1, 2, 3])
for (row, col), ax in np.ndenumerate(axs):  
    
    # # Plot PCO2 or po4
    var = vars2plot[row*2 + col]
    a = var.interp(lon=lon_of_interest).plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, add_colorbar = False)


    # # Plot isopycnals (all slightly different for each model)
    # rho_var = rho_vars[row*2 + col] - 1000
    
    if (row*2 + col) == 0:    # UVic
        red_levels = [36.0, 36.8]; yellow_levels = [36.5]
        
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        

        
    elif (row*2 + col) == 1:  # c-GENIE
        red_levels = [37, 38]; yellow_levels = [37.6]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
               
        
    elif (row*2 + col) == 2:  # LC
        red_levels = [37, 37.8]; yellow_levels = [37.6]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='r', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
                            
        
    elif (row*2 + col) == 3:  # LC-LGM
        red_levels = [33, 34.5]; yellow_levels = [34.0]
         
        # Add sigma2
        sigma_var = sigma2_vars[row*2 + col]
        var1 = sigma_var.interp(lon=lon_of_interest)
        CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
        CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)        
            



    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_xlim((-80,60)); ax.set_ylim((0,5000)); ax.invert_yaxis()
    ax.set_title(var.model, fontweight='bold')
    if ((row*2 + col) == 0) or ((row*2 + col) == 1):
        ax.set_xlabel('')
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r'Latitude [$\degree$N] at '+str(360-lon_of_interest)+r'$\degree$W')
    if ((row*2 + col) == 1) or ((row*2 + col) == 3):
        ax.set_ylabel('')
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel(r'Depth [m]')

    # Add a, b, c, d annotations
    if ((row*2 + col) == 0) or ((row*2 + col) == 2):
        left_just = -0.07
    else:
        left_just = -0.025
    ax.text(left_just, 1.1, panel_labels[row*2 + col], transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')





# Add common colorbar
cbar = plt.colorbar(a, cax=cbar_ax)
cbar.set_label(label=my_cb_label, fontsize=14, weight='bold', rotation=-90, va='bottom', labelpad=1)
cbar.ax.tick_params(labelsize=12)

# Adjust so panels don't overlap with colorbar
bottom, top = 0.1, 0.9
left, right = 0.1, 0.9
fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.25, wspace=0.1)


# Save fig
plt.savefig('../results/SuppFig_4panelPO4sigma2.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_4panelPO4sigma2.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_4panelPO4sigma2.pdf', dpi=600, bbox_inches='tight')


plt.show()


# %% Add sigmas to datasets and write as .nc files 


UV_ctrl['sigma1'] = sigma1_UV_ctrl
UV_ctrl['sigma2'] = sigma2_UV_ctrl 

UV_pmoc['sigma1'] = sigma1_UV_pmoc 
UV_pmoc['sigma2'] = sigma2_UV_pmoc 


LC_ctrl['sigma1'] = sigma1_LC_ctrl
LC_ctrl['sigma2'] = sigma2_LC_ctrl 

LC_pmoc['sigma1'] = sigma1_LC_pmoc 
LC_pmoc['sigma2'] = sigma2_LC_pmoc 


LGM_ctrl['sigma1'] = sigma1_LGM_ctrl
LGM_ctrl['sigma2'] = sigma2_LGM_ctrl 

LGM_pmoc['sigma1'] = sigma1_LGM_pmoc 
LGM_pmoc['sigma2'] = sigma2_LGM_pmoc 


GENIE_ctrl['sigma1'] = sigma1_GENIE_ctrl 
GENIE_ctrl['sigma2'] = sigma2_GENIE_ctrl 

GENIE_pmoc['sigma1'] = sigma1_GENIE_pmoc
GENIE_pmoc['sigma2'] = sigma2_GENIE_pmoc



# Save files

UV_ctrl.to_netcdf('../results/UV_ctrl.nc')
UV_pmoc.to_netcdf('../results/UV_pmoc.nc')

LC_ctrl.to_netcdf('../results/LC_ctrl.nc') # Don't actually need to save this output for anything
LC_pmoc.to_netcdf('../results/LC_pmoc.nc')

LGM_ctrl.to_netcdf('../results/LGM_ctrl.nc') # But do save this output, calc its outgassing for comparison later
LGM_pmoc.to_netcdf('../results/LGM_pmoc.nc')

GENIE_ctrl.to_netcdf('../results/GENIE_ctrl.nc') # Don't actually need to save this output for anything
GENIE_pmoc.to_netcdf('../results/GENIE_pmoc.nc')




