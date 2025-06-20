#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 14:01:15 2025

@author: mgs23
"""



# %% (0) Set up


# Import modules
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt

xr.set_options(keep_attrs=True)



# Note, all vars in UVic will have same dims, (depth: 19, latitude: 100, longitude: 100):
#     * longitude  (longitude) float64 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
#     * latitude   (latitude) float64 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1
#     * depth      (depth) float64 17.5 82.5 177.5 ... 4.658e+03 5.202e+03 5.778e+03



# %% - - - - - - - - - - - - - - - - - - - - - - - - - -

# %% (1) Read in output 


# Load in PCO2, DIC, Alk output

pathname = '../results/';
vars_ctrl = xr.open_dataset(pathname+'UV_ctrl.nc', decode_times=False)[['alk', 'dic', 'PCO2']]
vars_pmoc = xr.open_dataset(pathname+'UV_pmoc.nc', decode_times=False)[['alk', 'dic', 'PCO2']]



# Load in flux

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg.02001.01.01.nc'
raw_flux_ctrl = xr.open_dataset(pathname+filename, decode_times=False)
raw_flux_ctrl = raw_flux_ctrl.F_dic.isel(time=0)

# Convert from upwards to downward carbon flux (*-1), and convert units from
#  mol m-2 s-1 to mol m-2 yr-1
flux_ctrl = raw_flux_ctrl*(-1*365*24*60*60)
flux_ctrl.attrs['long_name'] = 'upwards carbon flux'
flux_ctrl.attrs['units'] = 'mol m-2 yr-1'
# Add note
flux_ctrl.attrs['note'] = 'Ctrl'
# Rename coordinates
flux_ctrl = flux_ctrl.rename({'latitude':'lat', 'longitude':'lon'})



pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
raw_flux_pmoc = xr.open_dataset(pathname+filename, decode_times=False) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years
raw_flux_pmoc = raw_flux_pmoc.F_dic.isel(time=-1)

# Convert from upwards to downward carbon flux (*-1), and convert units from
#  mol m-2 s-1 to mol m-2 yr-1
flux_pmoc = raw_flux_pmoc*(-1*365*24*60*60)
flux_pmoc.attrs['long_name'] = 'upwards carbon flux'
flux_pmoc.attrs['units'] = 'mol m-2 yr-1'
# Add note
flux_pmoc.attrs['note'] = 'Ctrl'
# Rename coordinates
flux_pmoc = flux_pmoc.rename({'latitude':'lat', 'longitude':'lon'})




# %% - - - - - - - - - - - - - - - - - - - - - - - - - -

# %% (2) SUPP FIG: PCO2 and Flux, absolute and anoms

my_fig_width  = 11.6
my_fig_height = 10

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(my_fig_width, my_fig_height))

overall_title = r'      PCO$_2$                                                Air-Sea CO$_2$ Flux'

# Left col: PCO2
col_num = 0
var2plot = 'PCO2'

for ii in np.arange(3):
    ax = axs[ii,col_num]
    if ii == 0:
        var = vars_ctrl[var2plot].isel(depth=0)
        my_title = 'UVic-ctrl'; panel_label = 'a'
        my_vmin = 150; my_vmax=500; my_cmap='viridis'
    elif ii == 1:
        var = vars_pmoc[var2plot].isel(depth=0)
        my_title = 'UVic-NP'; panel_label = 'b'
        my_vmin = 150; my_vmax=500; my_cmap='viridis'
    elif ii == 2:
        var = vars_pmoc[var2plot].isel(depth=0) - vars_ctrl[var2plot].isel(depth=0)
        var = var.assign_attrs(long_name = 'anomaly')
        my_title = 'Anomaly: UVic-NP  – UVic-ctrl'; panel_label = 'c'
        my_vmin = -100; my_vmax=100; my_cmap='seismic'

    var.plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap) #,
                                         # cbar_kwargs={'label': 'years'})

    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_title(my_title, fontweight='bold')
    ax.set_title(panel_label, fontweight='bold', loc='left')
    ax.set_ylabel(r'Latitude [$\degree$N]')
    if ii == 2:
        ax.set_xlabel(r'Longitude [$\degree$E]')
    else:
        ax.set_xlabel('')


# Right col: Flux
col_num = 1

for ii in np.arange(3):
    ax = axs[ii,col_num]
    if ii == 0:
        var = flux_ctrl
        my_title = 'UVic-ctrl'; panel_label = 'd'
        my_vmin = -5; my_vmax=5; my_cmap='RdBu_r'
    elif ii == 1:
        var = flux_pmoc
        my_title = 'UVic-NP'; panel_label = 'e'
        my_vmin = -5; my_vmax=5; my_cmap='RdBu_r'
    elif ii == 2:
        var = flux_pmoc - flux_ctrl
        var = var.assign_attrs(long_name = 'anomaly')
        my_title = 'Anomaly: UVic-NP  – UVic-ctrl'; panel_label = 'f'
        my_vmin = -4; my_vmax=4; my_cmap='seismic'
        

    var.plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap) #,
                                         # cbar_kwargs={'label': 'years'})

    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_title(my_title, fontweight='bold')
    ax.set_title(panel_label, fontweight='bold', loc='left')
    ax.set_ylabel('')
    if ii == 2:
        ax.set_xlabel(r'Longitude [$\degree$E]')
    else:
        ax.set_xlabel('')

        
fig.suptitle(overall_title, fontweight='bold', fontsize=18, y=0.98)

plt.tight_layout()


# Save fig
plt.savefig('../results/SuppFig_GlobMap_PCO2FluxAnoms.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_GlobMap_PCO2FluxAnoms.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_GlobMap_PCO2FluxAnoms.pdf', dpi=600, bbox_inches='tight')

plt.show()



# %% - - - - - - - - - - - - - - - - - - - - - - - - - -

# %% (3) Zonal integrals of flux and PCO2

# %% Construct masks for Pac, IndoPac, and Atlantic

# Generate cmpi6 basin masks to get data from just one basin
# And just multiply variables by these masks of 0s and 1s

import cartopy.crs as ccrs
from cmip_basins.basins import generate_basin_codes
import cmip_basins.cmip6 as cmip6
# Note, also requires install of "regionmasks" module


grid_UVic = xr.Dataset()
grid_UVic["lon"] = xr.DataArray(flux_ctrl.lon.values, dims=("lon"))
grid_UVic["lat"] = xr.DataArray(flux_ctrl.lat.values, dims=("lat"))


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



# %% Also need UVic cell areas

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
dxdy_model = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['G_areaT']
# Rename coordinates
dxdy_model = dxdy_model.rename({'longitude':'lon', 'latitude':'lat'})

# <xarray.DataArray 'G_areaT' (lat: 100, lon: 100)> Size: 40kB
# [10000 values with dtype=float32]
# Coordinates:
#   * lon      (lon) float64 800B 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
#   * lat      (lat) float64 800B -89.1 -87.3 -85.5 -83.7 ... 83.7 85.5 87.3 89.1
# Attributes:
#     long_name:    tracer grid area
#     units:        m2
 

# %% GLOBAL zonal integral

#  Cut off Arctic w/ 65 degN
n_lat = 65

# flux * dxddy_model (m2) to go from mol m-2 yr-1 --> mol yr-1, then sum up zonally
int_glob_flux_pmoc = (flux_pmoc*dxdy_model).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_glob_flux_ctrl = (flux_ctrl*dxdy_model).sel(lat=slice(-90, n_lat)).sum(dim='lon')


int_glob_PCO2_pmoc = vars_pmoc.PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_glob_PCO2_ctrl = vars_ctrl.PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')

avg_glob_PCO2_pmoc = vars_pmoc.PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')
avg_glob_PCO2_ctrl = vars_ctrl.PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')


# # Un-Comment for rough figs

# int_glob_flux_pmoc.plot()
# int_glob_flux_ctrl.plot()
# (int_glob_flux_pmoc-int_glob_flux_ctrl).plot()

# int_glob_PCO2_pmoc.plot()
# int_glob_PCO2_ctrl.plot()
# (int_glob_PCO2_pmoc-int_glob_PCO2_ctrl).plot()

# avg_glob_PCO2_pmoc.plot()
# avg_glob_PCO2_ctrl.plot()
# (avg_glob_PCO2_pmoc-avg_glob_PCO2_ctrl).plot()


# %% PACIFIC zonal integral

#  Cut off Arctic w/ 65 degN
n_lat = 65

int_Pac_flux_pmoc = (flux_pmoc*dxdy_model*UVic_Pac_mask).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_Pac_flux_ctrl = (flux_ctrl*dxdy_model*UVic_Pac_mask).sel(lat=slice(-90, n_lat)).sum(dim='lon')


int_Pac_PCO2_pmoc = (vars_pmoc*UVic_Pac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_Pac_PCO2_ctrl = (vars_ctrl*UVic_Pac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')

avg_Pac_PCO2_pmoc = (vars_pmoc*UVic_Pac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')
avg_Pac_PCO2_ctrl = (vars_ctrl*UVic_Pac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')



# %% INDO-PACIFIC zonal integral

#  Cut off Arctic w/ 65 degN
n_lat = 65

int_IndoPac_flux_pmoc = (flux_pmoc*dxdy_model*UVic_IndoPac_mask).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_IndoPac_flux_ctrl = (flux_ctrl*dxdy_model*UVic_IndoPac_mask).sel(lat=slice(-90, n_lat)).sum(dim='lon')


int_IndoPac_PCO2_pmoc = (vars_pmoc*UVic_IndoPac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')
int_IndoPac_PCO2_ctrl = (vars_ctrl*UVic_IndoPac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).sum(dim='lon')

avg_IndoPac_PCO2_pmoc = (vars_pmoc*UVic_IndoPac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')
avg_IndoPac_PCO2_ctrl = (vars_ctrl*UVic_IndoPac_mask).PCO2.isel(depth=0).sel(lat=slice(-90, n_lat)).mean(dim='lon')



# %% SUPP FIG: Zon-Int Flux as func(lat), IndoPac and Glob

my_fig_width  = 12
my_fig_height = 4

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(my_fig_width, my_fig_height))

for ii in np.arange(2):
    
    # Left panel: Indo-Pacific integral
    if ii == 0:
        ctrl = int_IndoPac_flux_ctrl
        pmoc = int_IndoPac_flux_pmoc
        anom = pmoc - ctrl
        panel_title = 'Indian-Pacific Basin'
        panel_label = 'a'
        panel_ylabel = r'upwards carbon flux [mol C yr$^{-1}$]'
        
    else:
        ctrl = int_glob_flux_ctrl
        pmoc = int_glob_flux_pmoc
        anom = pmoc - ctrl
        panel_title = 'Global'
        panel_label = 'b'
        panel_ylabel = ''

    ax = axs[ii]
    
    # Plot a zero-line
    ax.axhline(y=0, color='k', linewidth=0.5)
    
    # Plot ctrl (orange), pmoc (blue), and anomaly btw the two (dashed black)
    ctrl.plot(ax=ax, label="UVic-ctrl", color='C0')
    pmoc.plot(ax=ax, label="UVic-NP", color='orange')
    anom.plot(ax=ax, label="Anomaly (UVic-NP - UVic-ctrl)", color='k', linestyle='--', linewidth=0.9)
    
    # Shade the area between 0 and the anomaly line
    ax.fill_between(anom.lat, 0, anom.where(anom<=0), color='k', alpha=0.2)  # alpha controls transparency
    
    # Adding labels and title
    ax.set_ylim(-0.75e13, 0.9e13)    
    ax.set_xlabel(r'Latitude [$\degree$N]')
    ax.set_ylabel(panel_ylabel)
    ax.set_title(panel_title, fontweight='bold')
    ax.set_title(panel_label, loc='left', fontweight='bold')
    
    # Shift the scientific notation "1e13"
    t = ax.yaxis.get_offset_text()
    t.set_x(-0.1)


fig.suptitle('Zonally-Integrated Air-Sea CO$_2$ Flux', fontweight='bold', fontsize=13, y=1.03)

# Add legend to figure
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncols=3, bbox_to_anchor=(0.5, -0.12), fontsize=12, frameon=False)


# Save fig
plt.savefig('../results/SuppFig_ZonIntFlux_byLat.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_ZonIntFlux_byLat.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_ZonIntFlux_byLat.pdf', dpi=600, bbox_inches='tight')

plt.show()



