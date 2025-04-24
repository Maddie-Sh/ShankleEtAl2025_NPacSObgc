#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 15:52:35 2025

@author: mgs23
"""



# %% (0) Set up


# Import modules
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import PyCO2SYS as pyco2

xr.set_options(keep_attrs=True)


# Take S.O. averages below:
cutoff_lat = -40




# Note, all vars in UVic will have same dims, (depth: 19, latitude: 100, longitude: 100):
#     * longitude  (longitude) float64 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
#     * latitude   (latitude) float64 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1
#     * depth      (depth) float64 17.5 82.5 177.5 ... 4.658e+03 5.202e+03 5.778e+03



# %% - - - - - - - - - - - - - - - - - - - - - - - - - -

# %% (1) Radiocarbon 


# %% Open UVic D14C output

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
x = xr.open_dataset(pathname+filename, decode_times=False) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years

d14C_UV_pmoc  = x.O_dc14.isel(time=-1)
    # long_name:    delta carbon 14
    # units:        permil



pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg.02001.01.01.nc'
x = xr.open_dataset(pathname+filename, decode_times=False)
# time=1. array([2000.5]) (the end of the LGM simulation - our control)

d14C_UV_ctrl  = x.O_dc14.isel(time=0)



# Rename all dims to "depth", "lat", "lon"

d14C_UV_ctrl = d14C_UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})
d14C_UV_pmoc = d14C_UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})



# Plot to check how it looks 

plt.figure()
d14C_UV_ctrl.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('CONTROL: delta carbon 14')

plt.figure()
d14C_UV_pmoc.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot(vmin=-220, vmax=-80, cmap='Blues_r')
plt.gca().invert_yaxis()
plt.title('VENTILATED: delta carbon 14')

var2plot = (d14C_UV_pmoc -d14C_UV_ctrl)
plt.figure()
var2plot.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('ANOMALY (VENT - CTRL): delta carbon 14')



# %% Calculate radiocarbon age/years in each case

# Calculating absolute age from radiocarbon not possible/difficult as atmosphere
#  and ocean were not fully equilibrated and atmosphere D14C evolves over time 
#  (see in a timeseries of A_dc14 variable in tavg_tot.nc)

# Instead, can calculate radiocarbon years from the radiocarbon decay equation.
#  radiocarbon years = -8033*ln(d14C/1000+1) where 8033 is the mean lifetime of 
#  14C (Stuiver & Polach, 1977)
#    Stuiver, Minze, and Henry A. Polach. "Discussion reporting of 14C data." 
#    Radiocarbon 19.3 (1977): 355-363.

mean_lifetime = 8033

radCyr_UV_ctrl = -mean_lifetime * np.log((d14C_UV_ctrl/1000)+1)
radCyr_UV_pmoc = -mean_lifetime * np.log((d14C_UV_pmoc/1000)+1)
# (depth: 19, lat: 100, lon: 100)


# Re-assign long_name and units

radCyr_UV_ctrl = radCyr_UV_ctrl.assign_attrs(units="years", long_name="radiocarbon years")
radCyr_UV_pmoc = radCyr_UV_pmoc.assign_attrs(units="years", long_name="radiocarbon years")


# Plot to check how it looks 

plt.figure()
radCyr_UV_ctrl.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('CONTROL: radiocarbon years')

plt.figure()
radCyr_UV_pmoc.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('VENTILATED: radiocarbon years')

var2plot = (radCyr_UV_pmoc - radCyr_UV_ctrl)
plt.figure()
var2plot.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('ANOMALY (VENT - CTRL): radiocarbon years')




# %% Calc diff deep-minus-surface

# Then calculate the difference between each grid cell's radiocarbon age and 
# it's equivalent surface radiocarbon age at the same lat and long (do deep 
#  minus surface to make numbers positive). That's equivalent to a benthic 
#  plankton radiocarbon age and so could be compared directly to data 


radCage_minusSurf_UV_ctrl = radCyr_UV_ctrl - radCyr_UV_ctrl.sel(depth=0, method='nearest')
radCage_minusSurf_UV_pmoc = radCyr_UV_pmoc - radCyr_UV_pmoc.sel(depth=0, method='nearest')


# Re-assign attributes
radCage_minusSurf_UV_ctrl = radCage_minusSurf_UV_ctrl.assign_attrs(long_name = 'deep minus surface radiocarbon age')
radCage_minusSurf_UV_pmoc = radCage_minusSurf_UV_pmoc.assign_attrs(long_name = 'deep minus surface radiocarbon age')


# Plot to check how it looks 

plt.figure()
radCage_minusSurf_UV_ctrl.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('CONTROL: Deep-minus-surface radiocarbon years')

plt.figure()
radCage_minusSurf_UV_pmoc.sel(lon=200, method='nearest').sel(lat=slice(-60,60)).plot()
plt.gca().invert_yaxis()
plt.title('VENTILATED: Deep-minus-surface radiocarbon years')



# %% SUPP FIG: Anomaly in radC years (surface-deep)


lon_of_interest = 200    # 200 E = 160 W

my_fig_width  = 5.8
my_fig_height = 9


#  6.6, 10.5


# Aesthetics
my_vmin = -400; my_vmax = 50
my_cmap = 'seismic'
my_cb_label = 'Radiocarbon years anomaly'


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(my_fig_width, my_fig_height))

for ii in np.arange(3):
    ax = axs[ii]
    if ii == 0:
        var = radCage_minusSurf_UV_ctrl
        my_title = 'UVic-ctrl'; panel_label = 'a'
        my_vmin = 0; my_vmax=4000; my_cmap='Reds'
    elif ii == 1:
        var = radCage_minusSurf_UV_pmoc
        my_title = 'UVic-NP'; panel_label = 'b'
        my_vmin = 0; my_vmax=1200; my_cmap='Reds'
    elif ii == 2:
        var = radCage_minusSurf_UV_pmoc - radCage_minusSurf_UV_ctrl
        var = var.assign_attrs(long_name = 'anomaly')
        my_title = 'Anomaly: UVic-NP  – UVic-ctrl'; panel_label = 'c'
        my_vmin = -3000; my_vmax=0; my_cmap='Blues_r'

    var.interp(lon=lon_of_interest).plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap,
                                         cbar_kwargs={'label': '[radiocarbon years]'})

    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_xlim((-80,60)); ax.set_ylim((0,5000)); ax.invert_yaxis()
    ax.set_title(my_title, fontweight='bold')
    ax.set_title(panel_label, fontweight='bold', loc='left')

    ax.set_ylabel('Depth [m]')
    if ii == 2:
        ax.set_xlabel(r'Latitude [$\degree$N] at '+str(360-lon_of_interest)+r'$\degree$W')
    else:
        ax.set_xlabel('')
        
    fig.suptitle('Deep-minus-surface Radiocarbon Age', fontweight='bold', fontsize=14)

plt.tight_layout()




# Save fig
plt.savefig('../results/SuppFig_RadC_age_UVic.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_RadC_age_UVic.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_RadC_age_UVic.pdf', dpi=600, bbox_inches='tight')

plt.show()




# %% Misc: get Atms_dc14


pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tsi_tot.nc'
x = xr.open_dataset(pathname+filename, decode_times=False) 

x.A_dc14.plot()

# x.A_dc14.values[0]
# Out[20]: 384.73187
# # ctrl


# x.A_dc14.values[-1]
# Out[22]: -112.76231
# # exp/ventilated ^




# %% - - - - - - - - - - - - - - - - - - - - - - - - - -
# %% (2) Calc psi in UVic

# %% Load in necessary UVic data


# UVic (U-fwf)
pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'; 
data_UVic = xr.open_dataset(pathname+filename,decode_times=False);


# Need regular and Gent McWilliams meridional velocity
UVic_v_regular  = data_UVic.O_velY
UVic_gmv  = data_UVic.O_gmvelY

# Grid cell sizes
# UVic_dx   calc'd below by hand
UVic_dz   = data_UVic.G_dzt # (no dzu, only t. = "thickness t grid" - think corresponds with depth_w
                            # coordinate = depth. Corresponds to v's depths, not w's depths


# Rename dims to "lat" and "lon"
UVic_v_regular = UVic_v_regular.rename({'latitude_V':'lat', 'longitude_V':'lon'})
UVic_gmv = UVic_gmv.rename({'latitude_V':'lat', 'longitude_V':'lon'})


# Add regular and Gent McWilliams meridional velocities
UVic_v = UVic_v_regular + UVic_gmv




# %% Load in Control sim data

# UVic (U-fwf)
pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg.02001.01.01.nc'; 
data_UVic_ctrl = xr.open_dataset(pathname+filename,decode_times=False);


# Need regular and Gent McWilliams meridional velocity
UVic_v_regular_ctrl    = data_UVic_ctrl.O_velY
UVic_gmv_ctrl    = data_UVic_ctrl.O_gmvelY


# Rename dims to "lat" and "lon"
UVic_v_regular_ctrl = UVic_v_regular_ctrl.rename({'latitude_V':'lat', 'longitude_V':'lon'})
UVic_gmv_ctrl = UVic_gmv_ctrl.rename({'latitude_V':'lat', 'longitude_V':'lon'})


# Add regular and Gent McWilliams meridional velocities
UVic_v_ctrl = UVic_v_regular_ctrl + UVic_gmv_ctrl



# %% UVic_dx

radius_Earth_m = 6371000

num_lons_UVic = np.size(UVic_v.lon.values)
lats_UVic = UVic_v.lat.values


# - - - dx - - - # 
# tot circumference (as a function of lat) / number of cells around

UVic_dx_working = 2*np.pi*radius_Earth_m*np.cos(np.deg2rad(lats_UVic))/num_lons_UVic
# (72,)
UVic_dx = xr.DataArray(
    data=UVic_dx_working,
    dims=["lat"],
    coords=dict(
        lat=(["lat"], UVic_v_regular.lat.values),
    ),
    attrs=dict(
        description="velocity grid cell width dx",
        units="m",
    ),
)


# %% Generate cmpi6 basin masks to get data from just Pac basin

import cartopy.crs as ccrs
from cmip_basins.basins import generate_basin_codes
import cmip_basins.cmip6 as cmip6
# Note, also requires install of "regionmasks" module


grid_UVic = xr.Dataset()
grid_UVic["lon"] = xr.DataArray(UVic_v_regular.lon.values, dims=("lon"))
grid_UVic["lat"] = xr.DataArray(UVic_v_regular.lat.values, dims=("lat"))


codes_UVic = generate_basin_codes(grid_UVic, lon="lon", lat="lat", persian=True, style="cmip6")


UVic_Pac_mask = xr.where(codes_UVic == 3, 1, 0)
UVic_IndoPac_mask = xr.where((codes_UVic == 3) | (codes_UVic == 5), 1, 0)
UVic_Atl_mask = xr.where(codes_UVic == 2, 1, 0)
# ^ These each 0,1 'mask' (lat: 100, lon: 100)


# Add Southern Ocean to each mask
# hard-coded
# Pac: lon is 147.6 - 291.6 degE (indices 40:80), slice 40:81 does it
#      lat is all lats -37.8 (index 28) and below need to be made from 0 to 1
UVic_Pac_mask = xr.where( (UVic_Pac_mask.lon >147) & (UVic_Pac_mask.lon <292) & (UVic_Pac_mask.lat <-36), 1, UVic_Pac_mask)
UVic_Pac_mask = UVic_Pac_mask.transpose("lat", "lon")

# IndoPac: lon is 21.6 - 291.6degE (indices 5:80)
#          lat is all lats -30.6 and below need to be made from 0 to 1 
UVic_IndoPac_mask = xr.where( (UVic_IndoPac_mask.lon >21) & (UVic_IndoPac_mask.lon <292) & (UVic_IndoPac_mask.lat <-30), 1, UVic_IndoPac_mask)
UVic_IndoPac_mask = UVic_IndoPac_mask.transpose("lat", "lon")
# Good.

# Atl: lon is 18degE and below (index 4) and 295.2degE and above (index 81)
#      lat is all lats -34.2 and below need to be made from 0 to 1
UVic_Atl_mask = xr.where( ((UVic_Atl_mask.lon <=18) | (UVic_Atl_mask.lon >295)) & (UVic_Atl_mask.lat <-34), 1, UVic_Atl_mask)
UVic_Atl_mask = UVic_Atl_mask.transpose("lat", "lon")
# Good.



# %% Def Func for Stream Function

def calc_psi(v, dx, dz, west, east):
    """
    v :      4D xarray DataArray w/ dims time, depth (surf down), lat, lon
    dx % dz: 1D xarray DataArrays w/ dims lat or lon and depth, respectively
    west % east: scalar values of the longitudes across which to integrate zonally (degE)
    
    Returns a 3D xarray DataArray "psi" of dims time, depth, lat
    """
    
    # Integrate v zonally across the basin
    # xr.set_options(keep_attrs=True)

    v_dxdz = v*dx*dz
    
    zon_vdxdz = v_dxdz.sel(lon=slice(west,east)).sum(dim='lon')
    # leaving dimensions: time, depth (rows) and lat (cols)
    
    # Now sum up cumulatively from bottom
  
    
    # This should be adding up from bottom, with a row of 0s on the bottom, good.
    temp_psi = zon_vdxdz.sum(dim='depth') - zon_vdxdz.cumsum(dim='depth')
    # Re-order, chop off last row (0s), and re-assign depth coords to be one less
    temp_psi = temp_psi.transpose("time", "depth", "lat")
    temp_psi = temp_psi.isel(depth=slice(0,-1)) 
    temp_psi = temp_psi.assign_coords(depth=zon_vdxdz.depth.values[1:])

    # Take total sum & !give it a depth! of the uppermost depth value, !so that
    #  it can be appended to the above temp_psi
    sum_to_append = zon_vdxdz.sum(dim='depth')
    # Give it a depth dimension and assign it a coordinate (so that can be appended onto temp_psi)
    sum_to_append = sum_to_append.assign_coords(depth=zon_vdxdz.depth.values[0])
    sum_to_append = sum_to_append.expand_dims(dim='depth', axis=1)
    
    # Append the tot sum onto the top of the shortened temp_psi
    temp_psi_full = xr.concat([sum_to_append, temp_psi], dim='depth')
    
    # Convert to Sv (and *-1 as my calc seems inverse of LC output)
    psi_DataArray = temp_psi_full/1000000*-1
    
    psi_DataArray = psi_DataArray.assign_attrs(units="Sv", \
        description="Stream function over lons "+str(west)+":"+str(east)+" degE")    
        
    return psi_DataArray


# %% Calc psi 
# Takes a minute or two

# # # Indo-Pacific
# west = 50; east = 250  # (200E = 160W)  (260 = 100W)
# # Full Pacific
# west = 140; east = 250  # 155: (190E = 170W)  # (250E = 110W)
# # West Pacific
# west = 150; east = 190  # (190E = 170W)  # 150, 190
# # East Pacific
# west = 200; east = 260  # (200E = 160W)  (260 = 100W)

# Apply basin masks to meridional velocity first
# # Add regular and Gent McWilliams meridional velocities
# UVic_v = UVic_v_regular + UVic_gmv


# Pac
UVic_psi_Pac = calc_psi(UVic_v*UVic_Pac_mask, UVic_dx, UVic_dz, 0, 360)
UVic_psi_Pac_ctrl = calc_psi(UVic_v_ctrl*UVic_Pac_mask, UVic_dx, UVic_dz, 0, 360)
    
# IndoPac
UVic_psi_IndoPac = calc_psi(UVic_v*UVic_IndoPac_mask, UVic_dx, UVic_dz, 0, 360)
UVic_psi_IndoPac_ctrl = calc_psi(UVic_v_ctrl*UVic_IndoPac_mask, UVic_dx, UVic_dz, 0, 360)


# Atlantic
UVic_psi_Atl = calc_psi(UVic_v*UVic_Atl_mask, UVic_dx, UVic_dz, 0, 360)
UVic_psi_Atl_ctrl = calc_psi(UVic_v_ctrl*UVic_Atl_mask, UVic_dx, UVic_dz, 0, 360)





# %% Ex. code to plot psi over time

# # Atlantic
# model_to_plot = UVic_psi_Atl; model_name = 'My UVic Atl'; my_levels = [2, 4, 6, 8]; my_vmin=-10
# model_to_plot = UVic_psi_Atl_ctrl; model_name = 'My UVic Atl - CTRL'; my_levels = [2, 4, 6, 8, 10, 12, 14]; my_vmin=-10


# # Pacific
model_to_plot = UVic_psi_Pac; model_name = 'My UVic Pac'; my_levels = [2, 4, 6, 8, 10, 12, 14]; my_vmin=-10
# model_to_plot = UVic_psi_Pac_ctrl; model_name = 'My UVic Pac - CTRL'; my_levels = [2, 4, 6, 8, 10, 12, 14]; my_vmin=-10


# # Indo-Pacific
# model_to_plot = UVic_psi_IndoPac; model_name = 'My UVic IndoPac'; my_levels = [2, 4, 6, 8, 10, 12, 14]; my_vmin=-10
# model_to_plot = UVic_psi_IndoPac_ctrl; model_name = 'My UVic IndoPac - CTRL'; my_levels = [2, 4, 6, 8, 10, 12, 14]; my_vmin=-10


for ii in np.arange(model_to_plot.time.size):
    fig, ax = plt.subplots()
    CS0 = model_to_plot.isel(time=ii).plot.contour(colors=['k'], levels = [0], linewidths=1.0); plt.gca().clabel(CS0, fontsize=8)
    CS = model_to_plot.isel(time=ii).plot.contour(colors=['k'], levels = my_levels, linewidths=0.5); plt.gca().clabel(CS, fontsize=8)
    model_to_plot.isel(time=ii).plot(vmin=my_vmin, vmax=-my_vmin, cmap='bwr')
    
    plt.gca().invert_yaxis(); plt.xlim(-80,60)
    plt.gca().set_facecolor('grey')
    plt.gca().set_title(model_name+' Stream Func; time = '+str(model_to_plot.time.values[ii]))



# %% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# %% (3) Plot streamfunction of all models

# %% Load in LC Streamfunc output

pathname = '../data/Output_2014_LOVECLIM_L-fNA/';
# time = 2000 (annual resolution)


# This code makes a timeseries of the overturning stream function.
# Start by making ctrl pt (take average of first 10 years, annual resolution so avg t=0:9)
# Then build up a time series (taking avg last 10yr every 200 years)

# Atlantic 
LC_SFatl = xr.open_dataset(pathname+'Atl-SF-19-172.nc', decode_times=False)['ATL'].isel(TAXIS=slice(0,10)).mean(dim='TAXIS')
LC_SFatl = LC_SFatl.expand_dims(dim={"time": 1}, axis=0)    
LC_SFatl = LC_SFatl.assign_coords(time=[0.])

# Build LC_SFatl
# print("Building 'LC_SFatl':")
for ii in np.arange(10):
    start_index = ((ii+1)*200)-10; end_index = ((ii+1)*200)-1; time_val = (ii+1)*200;
    # print(time_val)
    # print('start_index:end_index = '+str(start_index)+':'+str(end_index))
    temp_sf = xr.open_dataset(pathname+'Atl-SF-19-172.nc', decode_times=False)['ATL'].isel(TAXIS=slice(start_index,end_index+1)).mean(dim='TAXIS')
    temp_sf = temp_sf.expand_dims(dim={"time": 1}, axis=0)    
    temp_sf = temp_sf.assign_coords(time=[time_val])
    LC_SFatl = xr.concat([LC_SFatl, temp_sf], dim='time')




# And Indo-Pacific
LC_SFpac = xr.open_dataset(pathname+'Pac-Sf-19-172.nc', decode_times=False)['PAC'].isel(TAXIS=slice(0,10)).mean(dim='TAXIS')
LC_SFpac = LC_SFpac.expand_dims(dim={"time": 1}, axis=0)    
LC_SFpac = LC_SFpac.assign_coords(time=[0.])

for ii in np.arange(10):
    start_index = ((ii+1)*200)-10; end_index = ((ii+1)*200)-1; time_val = (ii+1)*200;
    # print(time_val)
    # print('start_index:end_index = '+str(start_index)+':'+str(end_index))
    temp_sf = xr.open_dataset(pathname+'Pac-Sf-19-172.nc', decode_times=False)['PAC'].isel(TAXIS=slice(start_index,end_index+1)).mean(dim='TAXIS')
    temp_sf = temp_sf.expand_dims(dim={"time": 1}, axis=0)    
    temp_sf = temp_sf.assign_coords(time=[time_val])
    LC_SFpac = xr.concat([LC_SFpac, temp_sf], dim='time')
    

# Rename 'z' to 'depth'
LC_SFatl = LC_SFatl.rename({'ZAX':'depth', 'Y_AXIS':'lat'})
LC_SFpac = LC_SFpac.rename({'ZAX':'depth', 'Y_AXIS':'lat'})

# (time: 11, depth: 20, lat: 57)



# %% Load in LC-LGM Streamfunc output

# This code, unlike code for LC output above, doesn't make a timeseries but just
# loads in the ventilated-state overturning

# V3LNAw - "North Atlantic weak" - PMOC. 0.05Sv freshwater to N.Atl  
pathname = '../data/Output_2016_LOVECLIM-LGM/V3LNAw/';

# 200yrs in annual res, yr.800-1000. Take average of last 10yr
LCLGM_SFatl = xr.open_dataset(pathname+'Fatl_stream_inv.ascATL.cdf', decode_times=False)['ATL'].isel(TAXIS=slice(190,200)).mean(dim='TAXIS')
LCLGM_SFpac = xr.open_dataset(pathname+'Fpac_stream_inv.ascPAC.cdf', decode_times=False)['PAC'].isel(TAXIS=slice(190,200)).mean(dim='TAXIS')

# Rename 'z' to 'depth'
LCLGM_SFatl = LCLGM_SFatl.rename({'ZAX':'depth', 'Y_AXIS':'lat'})
LCLGM_SFpac = LCLGM_SFpac.rename({'ZAX':'depth', 'Y_AXIS':'lat'})



# V3L - "North Atlantic strong" - Ctrl.
pathname = '../data/Output_2016_LOVECLIM-LGM/V3L/';

# 200yrs in annual res, yr.800-1000. Take average of last 10yr
LCLGM_SFatl_ctrl = xr.open_dataset(pathname+'Fatl_stream_inv.ascATL.cdf', decode_times=False)['ATL'].isel(TAXIS=slice(190,200)).mean(dim='TAXIS')
LCLGM_SFpac_ctrl = xr.open_dataset(pathname+'Fpac_stream_inv.ascPAC.cdf', decode_times=False)['PAC'].isel(TAXIS=slice(190,200)).mean(dim='TAXIS')

# Rename 'z' to 'depth'
LCLGM_SFatl_ctrl = LCLGM_SFatl_ctrl.rename({'ZAX':'depth', 'Y_AXIS':'lat'})
LCLGM_SFpac_ctrl = LCLGM_SFpac_ctrl.rename({'ZAX':'depth', 'Y_AXIS':'lat'})

# (depth: 20, lat: 57)




# %% Load in GENIE Streamfunc output

# Control: Get data from:
pathname = '../data/Output_2020_cGENIE_Rae/Spin_LGM__SPIN_worjh2_Fe14C_preAge_Dye_LGMsave99/';


GENIE_SFatl_ctrl = xr.open_dataset(pathname+'fields_biogem_2d.nc')['phys_opsia'].isel(time=0)
GENIE_SFpac_ctrl = xr.open_dataset(pathname+'fields_biogem_2d.nc')['phys_opsip'].isel(time=0)


# (time: 1, zt_moc: 17, lat_moc: 37)
# Coordinates:
#     time     float64 999.5
#   * lat_moc  (lat_moc) float64 -90.0 -70.81 -62.73 -56.44 ... 62.73 70.81 90.0
#   * zt_moc   (zt_moc) float64 0.0 80.84 174.8 ... 3.576e+03 4.235e+03 5e+03
# Attributes:
#     long_name:      Atlantic streamfunction
#     standard_name:  ocean_meridional_overturning_streamfunction
#     units:          Sv
# Attributes:
#     long_name:      Pacific streamfunction
#     standard_name:  ocean_meridional_overturning_streamfunction
#     units:          Sv



# Ventilated: Get data from:
pathname = '../data/Output_2020_cGENIE_Rae/PA-15_negpt28Sv/';

GENIE_SFatl_exp = xr.open_dataset(pathname+'fields_biogem_2d.nc')['phys_opsia'].isel(time=-1)
GENIE_SFpac_exp = xr.open_dataset(pathname+'fields_biogem_2d.nc')['phys_opsip'].isel(time=-1)

# (time: 12, zt_moc: 17, lat_moc: 37)
# Coordinates:
#   * time     (time) float64 0.5 1.5 4.5 9.5 19.5 ... 499.5 999.5 2e+03 5e+03
#   * lat_moc  (lat_moc) float64 -90.0 -70.81 -62.73 -56.44 ... 62.73 70.81 90.0
#   * zt_moc   (zt_moc) float64 0.0 80.84 174.8 ... 3.576e+03 4.235e+03 5e+03
# Attributes:
#     long_name:      Atlantic streamfunction
#     standard_name:  ocean_meridional_overturning_streamfunction
#     units:          Sv
# Attributes:
#     long_name:      Pacific streamfunction
#     standard_name:  ocean_meridional_overturning_streamfunction
#     units:          Sv




# Rename 'zt_moc' to 'depth' and 'lat_moc' to 'lat'
GENIE_SFatl_ctrl = GENIE_SFatl_ctrl.rename({'zt_moc':'depth', 'lat_moc':'lat'})
GENIE_SFpac_ctrl = GENIE_SFpac_ctrl.rename({'zt_moc':'depth', 'lat_moc':'lat'})

GENIE_SFatl_exp = GENIE_SFatl_exp.rename({'zt_moc':'depth', 'lat_moc':'lat'})
GENIE_SFpac_exp = GENIE_SFpac_exp.rename({'zt_moc':'depth', 'lat_moc':'lat'})




# %% Before fig, load in PCO2 output (for use for bathy on figure)

# These files generated from script 01 and contain various tracers, alk, dic, PCO2 etc.

# Get data from:
pathname = '../results/';
UV_forBathy = xr.open_dataset(pathname+'UV_ctrl.nc', decode_times=False)['PCO2']
LC_forBathy = xr.open_dataset(pathname+'LC_ctrl.nc', decode_times=False)['PCO2']
LGM_forBathy = xr.open_dataset(pathname+'LGM_ctrl.nc', decode_times=False)['PCO2']
GENIE_forBathy = xr.open_dataset(pathname+'GENIE_ctrl.nc', decode_times=False)['PCO2']


# And interpret onto the coordinates of the psi data
UV_forBathy = UV_forBathy.interp(lat=UVic_psi_Pac.lat, depth=UVic_psi_Pac.depth, kwargs={"fill_value": "extrapolate"})
LC_forBathy = LC_forBathy.interp(lat=LC_SFpac.lat, depth=LC_SFpac.depth, kwargs={"fill_value": "extrapolate"})
LGM_forBathy = LGM_forBathy.interp(lat=LCLGM_SFpac.lat, depth=LCLGM_SFpac.depth, kwargs={"fill_value": "extrapolate"})
GENIE_forBathy = GENIE_forBathy.interp(lat=GENIE_SFpac_exp.lat, depth=GENIE_SFpac_exp.depth, kwargs={"fill_value": "extrapolate"})



# %% SUPP FIG: 8-panel O.T., ctrl and vent


# Note on:   for (col, row), ax in np.ndenumerate(np.transpose(axs))
# If axs=(nrows=2, ncols=4), then this^ will return (row,col) indices that go:
#   1  3    5  7  
#   2  4    6  8
# I.e., (row,col) = (0,0), (1,0), (0,1), (1,1), (0,2), (1,2), (0,3), (1,3)
# 
#  Want lists to therefore be: UVic, LC, c-GENIE, LC-LGM, and then repeated w/ vent sims


lon_of_interest = 200

my_fig_width  = 17    # 22, 7
my_fig_height = 5

my_title = 'Pacific Meridional Overturning Stream Function (Sv)'

# CONTROL - PACIFIC OVERTURNING
psi_vars_ctrl = [UVic_psi_Pac_ctrl.isel(time=0), LC_SFpac.isel(time=0), GENIE_SFpac_ctrl, LCLGM_SFpac_ctrl]
# VENTILATED - PACIFIC OVERTURNING
psi_vars_vent = [UVic_psi_Pac.isel(time=-1), LC_SFpac.isel(time=-1), GENIE_SFpac_exp, LCLGM_SFpac]

# LIST OF 8
psi_vars = psi_vars_ctrl + psi_vars_vent
bathy_vars = [UV_forBathy, LC_forBathy, GENIE_forBathy, LGM_forBathy, UV_forBathy, LC_forBathy, GENIE_forBathy, LGM_forBathy]
titles = ['UVic-ctrl','LOVECLIM','c-GENIE','LOVECLIM-LGM','UVic-NP','LOVECLIM','c-GENIE','LOVECLIM-LGM']
panel_labels = ['a', 'c', 'b', 'd', 'e', 'g', 'f', 'h'] # this is correct


fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(my_fig_width, my_fig_height), gridspec_kw={'hspace': 0.2, 'wspace': 0.08})


for (col, row), ax in np.ndenumerate(np.transpose(axs)):
    # print('row = '+str(row))
    # print('col = '+str(col))
    ii = row + 2*col  # == 0, 1, 2, 3, 4, 5, 6, 7
    # print(ii)
    # print('\n')


    # # # If want to plot representative bathymetry at lon_of_interest, un-comment
    # # # this chunk and also one other bit just below
    # # Plot bathymetry: plot some output that has NaNs for bathymetry, so that
    # #  bathymetry shows up in figure.
    # # Interp to lon of interest and make all values that are not nans 0
    # varB = bathy_vars[ii]
    # varB = varB.where(np.isnan(varB), 0)
    # varB.interp(lon=lon_of_interest).plot(ax=ax, cmap='RdBu', add_colorbar=False)


    # Plot stream function
    var = psi_vars[ii]

    # # # If want to plot representative bathymetry at lon_of_interest, un-comment
    # # # this chunk and also one other bit just above    
    # # First make stream function NaN where bathymetry is, so can see bathymetry
    # var = var.where(~np.isnan(varB.interp(lon=lon_of_interest)), np.nan)
    
    # Now plot
    im = var.sel(lat=slice(-40, 90)).plot(ax=ax, vmin=-20, vmax=20, cmap='RdBu_r', add_colorbar=False)


    # Aesthetics
    ax.set_facecolor('lightgrey')
    ax.set_xlim((-60,60)); ax.set_ylim((0,5000)); ax.invert_yaxis()
    ax.set_title(titles[ii], fontweight='bold')
    # Add a, b, c, d annotations
    ax.set_title(panel_labels[ii], fontweight='bold', loc='left')

    if (ii == 0) or (ii == 2) or (ii == 4) or (ii == 6):
        ax.set_xlabel('')
        ax.tick_params(labelbottom='')
    else:
        ax.set_xlabel(r'Latitude [$\degree$N]')

    if (ii == 0) or (ii == 1):
        ax.set_ylabel(r'Depth [m]')    
    else:
        ax.set_ylabel('')
        ax.tick_params(labelleft=False)


plt.tight_layout()

# Shift panels left and add common colorbar
fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.87, 0.11, 0.013, 0.77])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label("[Sv]", fontsize=13)
cb.ax.tick_params(labelsize=12)


# Now shift left-most four Ctrl panels more to the left, 1, 2, 5,6 ii=0,1,4,5
move_by = 0.02

# Move panel a
pos = axs[0,0].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[0,0].set_position(new_pos)

# Move panel b
pos = axs[0,1].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[0,1].set_position(new_pos)

# Move panel c
pos = axs[1,0].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[1,0].set_position(new_pos)

# Move panel d
pos = axs[1,1].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[1,1].set_position(new_pos)


# Add figure title
fig.suptitle(my_title, fontweight='bold', fontsize=18, y = 1.06)

# Add "Control"/"Ventilated" annotations
yfrac=0.94
ax.annotate('Control Simulations',
            xy=(0.14, yfrac), xycoords='figure fraction', fontsize=16, fontweight='bold') #horizontalalignment='left', verticalalignment='top',
ax.annotate('Ventilated Simulations',
            xy=(0.53, yfrac), xycoords='figure fraction', fontsize=16, fontweight='bold') #horizontalalignment='left', verticalalignment='top',





# Save fig
plt.savefig('../results/SuppFig_PacSF_CtrlAndVent.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_PacSF_CtrlAndVent.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_PacSF_CtrlAndVent.pdf', dpi=600, bbox_inches='tight')

plt.show()



# %% optional fig: Atlantic 8-panel O.T., ctrl and vent


# Note on:   for (col, row), ax in np.ndenumerate(np.transpose(axs))
# If axs=(nrows=2, ncols=4), then this^ will return (row,col) indices that go:
#   1  3    5  7  
#   2  4    6  8
# I.e., (row,col) = (0,0), (1,0), (0,1), (1,1), (0,2), (1,2), (0,3), (1,3)
# 
#  Want lists to therefore be: UVic, LC, c-GENIE, LC-LGM, and then repeated w/ vent sims


lon_of_interest = 200

my_fig_width  = 17    # 22, 7
my_fig_height = 5

my_title = 'Atlantic Meridional Overturning Stream Function (Sv)'

# CONTROL - Atlantic OVERTURNING
psi_vars_ctrl = [UVic_psi_Atl_ctrl.isel(time=0), LC_SFatl.isel(time=0), GENIE_SFatl_ctrl, LCLGM_SFatl_ctrl]
# VENTILATED - PACIFIC OVERTURNING
psi_vars_vent = [UVic_psi_Atl.isel(time=-1), LC_SFatl.isel(time=-1), GENIE_SFatl_exp, LCLGM_SFatl]

# LIST OF 8
psi_vars = psi_vars_ctrl + psi_vars_vent
bathy_vars = [UV_forBathy, LC_forBathy, GENIE_forBathy, LGM_forBathy, UV_forBathy, LC_forBathy, GENIE_forBathy, LGM_forBathy]
titles = ['UVic-ctrl','LOVECLIM','c-GENIE','LOVECLIM-LGM','UVic-NP','LOVECLIM','c-GENIE','LOVECLIM-LGM']
panel_labels = ['a', 'c', 'b', 'd', 'e', 'g', 'f', 'h'] # this is correct


fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(my_fig_width, my_fig_height), gridspec_kw={'hspace': 0.2, 'wspace': 0.08})


for (col, row), ax in np.ndenumerate(np.transpose(axs)):
    # print('row = '+str(row))
    # print('col = '+str(col))
    ii = row + 2*col  # == 0, 1, 2, 3, 4, 5, 6, 7
    # print(ii)
    # print('\n')


    # # # If want to plot representative bathymetry at lon_of_interest, un-comment
    # # # this chunk and also one other bit just below
    # # Plot bathymetry: plot some output that has NaNs for bathymetry, so that
    # #  bathymetry shows up in figure.
    # # Interp to lon of interest and make all values that are not nans 0
    # varB = bathy_vars[ii]
    # varB = varB.where(np.isnan(varB), 0)
    # varB.interp(lon=lon_of_interest).plot(ax=ax, cmap='RdBu', add_colorbar=False)


    # Plot stream function
    var = psi_vars[ii]

    # # # If want to plot representative bathymetry at lon_of_interest, un-comment
    # # # this chunk and also one other bit just above    
    # # First make stream function NaN where bathymetry is, so can see bathymetry
    # var = var.where(~np.isnan(varB.interp(lon=lon_of_interest)), np.nan)
    
    # Now plot
    im = var.sel(lat=slice(-40, 90)).plot(ax=ax, vmin=-22, vmax=22, cmap='RdBu_r', add_colorbar=False)


    # Aesthetics
    ax.set_facecolor('lightgrey')
    ax.set_xlim((-60,60)); ax.set_ylim((0,5000)); ax.invert_yaxis()
    ax.set_title(titles[ii], fontweight='bold')
    # Add a, b, c, d annotations
    ax.set_title(panel_labels[ii], fontweight='bold', loc='left')

    if (ii == 0) or (ii == 2) or (ii == 4) or (ii == 6):
        ax.set_xlabel('')
        ax.tick_params(labelbottom='')
    else:
        ax.set_xlabel(r'Latitude [$\degree$N]')

    if (ii == 0) or (ii == 1):
        ax.set_ylabel(r'Depth [m]')    
    else:
        ax.set_ylabel('')
        ax.tick_params(labelleft=False)


plt.tight_layout()

# Shift panels left and add common colorbar
fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.87, 0.11, 0.013, 0.77])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label("[Sv]", fontsize=13)
cb.ax.tick_params(labelsize=12)


# Now shift left-most four Ctrl panels more to the left, 1, 2, 5,6 ii=0,1,4,5
move_by = 0.02

# Move panel a
pos = axs[0,0].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[0,0].set_position(new_pos)

# Move panel b
pos = axs[0,1].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[0,1].set_position(new_pos)

# Move panel c
pos = axs[1,0].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[1,0].set_position(new_pos)

# Move panel d
pos = axs[1,1].get_position()
new_pos = [pos.x0-move_by, pos.y0, pos.width, pos.height]
axs[1,1].set_position(new_pos)


# Add figure title
fig.suptitle(my_title, fontweight='bold', fontsize=18, y = 1.06)

# Add "Control"/"Ventilated" annotations
yfrac=0.94
ax.annotate('Control Simulations',
            xy=(0.14, yfrac), xycoords='figure fraction', fontsize=16, fontweight='bold') #horizontalalignment='left', verticalalignment='top',
ax.annotate('Ventilated Simulations',
            xy=(0.53, yfrac), xycoords='figure fraction', fontsize=16, fontweight='bold') #horizontalalignment='left', verticalalignment='top',


# %%  - - - - - - - - - - - - - - - - - - 

# %% (4) Calculate max overturning SF each model - Atl and Pac

# Calc max of zonally-integrated Pacific meridional overturning streamfunction

lat_above = 31  # 31 degN cuts out subtropical gyre circ in all 4 simulations
                # 30-vs-31 starts to get into influence of STG in LC sim, but 
                # beyond this, results insensitive to exact lat (30-vs-31) in all models 

# UVic
print('UVic')

maxSF_UV_ctrl_Pac = UVic_psi_Pac_ctrl.isel(time=0).sel(lat=slice(lat_above, 91)).max().values
maxSF_UV_vent_Pac = UVic_psi_Pac.isel(time=-1).sel(lat=slice(lat_above, 91)).max().values
print('\nPacific:')
print(f'UVic ctrl = {maxSF_UV_ctrl_Pac:.3f} Sv')
print(f'UVic Vent = {maxSF_UV_vent_Pac:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_UV_vent_Pac-maxSF_UV_ctrl_Pac):.3f}')


maxSF_UV_ctrl_Atl = UVic_psi_Atl_ctrl.isel(time=0).sel(lat=slice(lat_above, 91)).max().values
maxSF_UV_vent_Atl = UVic_psi_Atl.isel(time=-1).sel(lat=slice(lat_above, 91)).max().values
print('\nAtlantic:')
print(f'UVic ctrl = {maxSF_UV_ctrl_Atl:.3f} Sv')
print(f'UVic Vent = {maxSF_UV_vent_Atl:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_UV_vent_Atl-maxSF_UV_ctrl_Atl):.3f}')




# LC
print('\n\n\nLOVECLIM')

maxSF_LC_ctrl_Pac = LC_SFpac.isel(time=0).sel(lat=slice(lat_above, 91)).max().values
maxSF_LC_vent_Pac = LC_SFpac.isel(time=-1).sel(lat=slice(lat_above, 91)).max().values
print('\nPacific:')
print(f'LC ctrl = {maxSF_LC_ctrl_Pac:.3f} Sv')
print(f'LC Vent = {maxSF_LC_vent_Pac:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_LC_vent_Pac-maxSF_LC_ctrl_Pac):.3f}')

maxSF_LC_ctrl_Atl = LC_SFatl.isel(time=0).sel(lat=slice(lat_above, 91)).max().values
maxSF_LC_vent_Atl = LC_SFatl.isel(time=-1).sel(lat=slice(lat_above, 91)).max().values
print('\nAtlantic:')
print(f'LC ctrl = {maxSF_LC_ctrl_Atl:.3f} Sv')
print(f'LC Vent = {maxSF_LC_vent_Atl:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_LC_vent_Atl-maxSF_LC_ctrl_Atl):.3f}')





# LC-LGM
print('\n\n\nLOVECLIM-LGM')

maxSF_LCLGM_ctrl_Pac = LCLGM_SFpac_ctrl.sel(lat=slice(lat_above, 91)).max().values
maxSF_LCLGM_vent_Pac = LCLGM_SFpac.sel(lat=slice(lat_above, 91)).max().values
print('\nPacific:')
print(f'LC-LGM ctrl = {maxSF_LCLGM_ctrl_Pac:.3f} Sv')
print(f'LC-LGM Vent = {maxSF_LCLGM_vent_Pac:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_LCLGM_vent_Pac-maxSF_LCLGM_ctrl_Pac):.3f}')

maxSF_LCLGM_ctrl_Atl = LCLGM_SFatl_ctrl.sel(lat=slice(lat_above, 91)).max().values
maxSF_LCLGM_vent_Atl = LCLGM_SFatl.sel(lat=slice(lat_above, 91)).max().values
print('\nAtlantic:')
print(f'LC-LGM ctrl = {maxSF_LCLGM_ctrl_Atl:.3f} Sv')
print(f'LC-LGM Vent = {maxSF_LCLGM_vent_Atl:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_LCLGM_vent_Atl-maxSF_LCLGM_ctrl_Atl):.3f}')




# c-GENIE
print('\n\n\nc-GENIE')

maxSF_GENIE_ctrl_Pac = GENIE_SFpac_ctrl.sel(lat=slice(lat_above, 91)).max().values
maxSF_GENIE_vent_Pac = GENIE_SFpac_exp.sel(lat=slice(lat_above, 91)).max().values
print('\nPacific:')
print(f'GENIE ctrl = {maxSF_GENIE_ctrl_Pac:.3f} Sv')
print(f'GENIE Vent = {maxSF_GENIE_vent_Pac:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_GENIE_vent_Pac-maxSF_GENIE_ctrl_Pac):.3f}')

maxSF_GENIE_ctrl_Atl = GENIE_SFatl_ctrl.sel(lat=slice(lat_above, 91)).max().values
maxSF_GENIE_vent_Atl = GENIE_SFatl_exp.sel(lat=slice(lat_above, 91)).max().values
print('\nAtlantic:')
print(f'GENIE ctrl = {maxSF_GENIE_ctrl_Atl:.3f} Sv')
print(f'GENIE Vent = {maxSF_GENIE_vent_Atl:.3f} Sv')
print(f'Anom (vent-ctrl) = {(maxSF_GENIE_vent_Atl-maxSF_GENIE_ctrl_Atl):.3f}')


