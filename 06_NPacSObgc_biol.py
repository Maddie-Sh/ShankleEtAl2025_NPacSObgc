#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 13:30:52 2025

@author: mgs23
"""


# %% (0) Set up


# Import modules
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import PyCO2SYS as pyco2
import gsw

xr.set_options(keep_attrs=True)



# %% (1) Load in PO4 output 

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
UV_pmoc = xr.open_dataset(pathname+filename, decode_times=False)[['O_po4']].isel(time=-1) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years

filename = 'tavg.02001.01.01.nc'
UV_ctrl = xr.open_dataset(pathname+filename, decode_times=False)[['O_po4']].isel(time=0)
# time=1. array([2000.5]) (the end of the LGM simulation - our control)


# O_po4
    # valid_range:  [ -1. 100.]
    # FillValue:    9.96921e+36
    # long_name:    O_po4
    # units:        mol m-3



# Rename all dims to "depth", "lat", "lon"
UV_pmoc = UV_pmoc.rename({'longitude':'lon', 'latitude':'lat'})
UV_ctrl = UV_ctrl.rename({'longitude':'lon', 'latitude':'lat'})




# %% (2) Load in rho output

pathname = '../results/';
filename = 'UV_pmoc.nc'
UV_pmoc_rho = xr.open_dataset(pathname+filename, decode_times=False)['rho']
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years

filename = 'UV_ctrl.nc'
UV_ctrl_rho = xr.open_dataset(pathname+filename, decode_times=False)['rho']



# %% (3) Convert PO4 units mol m-3 --> assign_attrs(units = 'micro-moles kg-1')

UV_pmoc_po4 = UV_pmoc['O_po4']*(1/UV_pmoc_rho)*(10**6)
UV_ctrl_po4 = UV_ctrl['O_po4']*(1/UV_ctrl_rho)*(10**6)


# Re-assign attributes

UV_pmoc_po4 = UV_pmoc_po4.assign_attrs(units = 'μmol kg-1', long_name = 'phosphate', standard_name = 'phosphate')
UV_ctrl_po4 = UV_ctrl_po4.assign_attrs(units = 'μmol kg-1', long_name = 'phosphate', standard_name = 'phosphate')




# %% (4) Load in NPP 


# long_name:    ocean net primary production rate
# units:        mol N m-3 s-1

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
filename = 'tavg_tot.nc'
UV_pmoc_NPP = xr.open_dataset(pathname+filename, decode_times=False)[['O_phytnpp']].isel(time=-1) 
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
# Take last timestep at the end of the 1000 years

filename = 'tavg.02001.01.01.nc'
UV_ctrl_NPP = xr.open_dataset(pathname+filename, decode_times=False)[['O_phytnpp']].isel(time=0)
# time=1. array([2000.5]) (the end of the LGM simulation - our control)

# mol N m-3 s-1 --> yr-1
UV_pmoc_NPP = UV_pmoc_NPP.O_phytnpp *(365*24*60*60)
UV_ctrl_NPP = UV_ctrl_NPP.O_phytnpp *(365*24*60*60)

# Rename all dims to "depth", "lat", "lon"
UV_pmoc_NPP = UV_pmoc_NPP.rename({'npzdlev':'depth', 'longitude':'lon', 'latitude':'lat'})
UV_ctrl_NPP = UV_ctrl_NPP.rename({'npzdlev':'depth', 'longitude':'lon', 'latitude':'lat'})
# Re-assign attributes
UV_pmoc_NPP = UV_pmoc_NPP.assign_attrs(units = 'mol N m-3 yr-1')
UV_ctrl_NPP = UV_ctrl_NPP.assign_attrs(units = 'mol N m-3 yr-1')





# %% (5) SUPP FIG: NPP and PO43-, absolute and anoms

my_fig_width  = 11.6
my_fig_height = 10

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(my_fig_width, my_fig_height))

overall_title = r'PO$_4^{3-}$                                                       NPP         '

# Left col: PCO2
col_num = 0
var2plot = r'PO$_4^{3-}$'

for ii in np.arange(3):
    ax = axs[ii,col_num]
    if ii == 0:
        var = UV_ctrl_po4.isel(depth=0)
        my_title = 'UVic-ctrl'; panel_label = 'a'
# %%
        my_vmin = 0; my_vmax=2.5; my_cmap='viridis'

    elif ii == 1:
        var = UV_pmoc_po4.isel(depth=0)
        my_title = 'UVic-NP'; panel_label = 'b'
        my_vmin = 0; my_vmax=2.5; my_cmap='viridis'
    elif ii == 2:
        var = UV_pmoc_po4.isel(depth=0) - UV_ctrl_po4.isel(depth=0)
        var = var.assign_attrs(long_name = 'anomaly')
        my_title = 'Anomaly: UVic-NP  – UVic-ctrl'; panel_label = 'c'
        my_vmin = -1.5; my_vmax=1.5; my_cmap='seismic'

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
        var = UV_ctrl_NPP.isel(depth=0)
        my_title = 'UVic-ctrl'; panel_label = 'd'
        my_vmin = 0; my_vmax=0.15; my_cmap='viridis'
    elif ii == 1:
        var = UV_pmoc_NPP.isel(depth=0)
        my_title = 'UVic-NP'; panel_label = 'e'
        my_vmin = 0; my_vmax=0.15; my_cmap='viridis'
    elif ii == 2:
        var = UV_pmoc_NPP.isel(depth=0) - UV_ctrl_NPP.isel(depth=0)
        var = var.assign_attrs(long_name = 'anomaly')
        my_title = 'Anomaly: UVic-NP  – UVic-ctrl'; panel_label = 'f'
        my_vmin = -0.06; my_vmax=0.06; my_cmap='seismic'
        

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
plt.savefig('../results/SuppFig_GlobMap_NPPandPO43.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_GlobMap_NPPandPO43.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_GlobMap_NPPandPO43.pdf', dpi=600, bbox_inches='tight')

plt.show()



