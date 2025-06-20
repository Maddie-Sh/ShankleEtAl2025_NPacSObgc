# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 13:20:55 2023

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





# %% (1) Load in output

# Get data from:
pathname = '../results/';
UV_ctrl = xr.open_dataset(pathname+'UV_ctrl.nc', decode_times=False)
UV_pmoc = xr.open_dataset(pathname+'UV_pmoc.nc', decode_times=False)




# %% (2) Get weights from cos(lat), areas, vol

# https://docs.xarray.dev/en/latest/examples/area_weighted_temperature.html
# "For a rectangular grid the cosine of the latitude is proportional to the 
#   grid cell area."

# Maddie checked, weighting by cell area the same as weighting by cos(lat value)

weights = np.cos(np.deg2rad(UV_ctrl.lat))
weights.name = "weights"
del weights.attrs['units'] 
del weights.attrs['edges'] 
del weights.attrs['standard_name'] 
weights.attrs['long_name'] = 'weights'



# - - - - - AREA dx*dy - - - - - - #

# # Maddie checked:
#  - G_areaT is within < 0.0001% of (G_dxt*cos(deg2rad(lat)))*G_dyt
# It does appear that dxt is only an equatorial value of longitudinal width,
#  and therefore must be multiplied by cos(lat) to get lat-dep dx. See also dx
#  described in documentation as "longitudinal width of equatorial cell" 
#  (Pacanowski 1995, pg. 45)
pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
dxdy_model = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['G_areaT']
# Rename coordinates
dxdy_model = dxdy_model.rename({'longitude':'lon', 'latitude':'lat'})


# - - - - - THICKNESS dz - - - - - #

dz_model = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['G_dzt']



# - - - - VOLUME dz*dx*dy - - - - #

# Expand dxdy_model (latitude, longitude) down in depth dimension (n=19)

dV_model = dxdy_model.expand_dims(dim={"depth": 19}).assign_coords({"depth": dz_model.depth.values})
dV_model = dV_model*dz_model
dV_model.attrs['long_name'] = 'tracer grid volume'
dV_model.attrs['units'] = 'm3'



# - - - - - AREA dx*dz - - - - - - #

dxdz_model = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['G_dxt']
# (longitude: 100) width t grid [m] ... but only at eqtr

# Expand out and assign dims to lat and depth
dxdz_model = dxdz_model.expand_dims(dim={"latitude": 100})
dxdz_model = dxdz_model.expand_dims(dim={"depth": 19})

dxdz_model = dxdz_model.assign_coords({"depth": dV_model.depth.values})
dxdz_model = dxdz_model.assign_coords({"latitude": dV_model.lat.values})

dxdz_model = dxdz_model*np.cos(np.deg2rad(dxdz_model.latitude))*dz_model
# dx * cos(lat) modulates width dx as a function of latitude
# then * dz

# Rename coordinates
dxdz_model = dxdz_model.rename({'longitude':'lon', 'latitude':'lat'})




# %% (3) Func for taking NPac avg PCO2

def PacBasin_fig(var, var_descr, fig_width, fig_height, lon_of_interest, my_vmin=None, my_vmax=None, my_cmap='seismic', depth_top=None, depth_btm=None, lat_S=None, lat_N=None, lon_E=140, lon_W=240, box_color='yellow', box_lw=2, var_depth_top=None, var_depth_btm=None, var_lat_S=None, var_lat_N=None, var_lon_E=140, var_lon_W=240, var_box_color='lime', weights=None):
    """
    This func will make one plot.
    - 'var' must have dims named 'depth', 'lat', and 'lon'.
    - 'var_descr' is a string simply describing the variable 'var'.
    - 'weights' must be 3D of same dims as var. Will usually be volume of cells,
        i.e. dVol_model.
    -  Note, matplotlib will automatically choose vmin/vmax if none are specified
        by user in calling the function.
    - It will also draw boxes, if depth_top/btm & lat_S/N are specified
    """
    fig, ax = plt.subplots(figsize=(fig_width, fig_height)) # 6.4, 4.8 default
    var.interp(lon=lon_of_interest).plot(vmin=my_vmin, vmax=my_vmax, cmap='seismic', cbar_kwargs={'label': var.units})

    # Add yellow Delta_PCO2 box (box corresponding to max PCO2 anom, same box in all plots)
    if depth_top is not None:
        ax.plot([lat_S, lat_N, lat_N, lat_S, lat_S], \
                [depth_btm, depth_btm, depth_top, depth_top, depth_btm], \
                 color=box_color, linewidth=box_lw*2) 
        # Also print out averages
        var_subset     = var.sel(depth=slice(depth_top, depth_btm), lat=slice(lat_S, lat_N), lon=slice(lon_E, lon_W))
        weights_subset = weights.sel(depth=slice(depth_top, depth_btm), lat=slice(lat_S, lat_N), lon=slice(lon_E, lon_W))
        
        print('\n'+var_descr+' '+var.units)
        print('Average over PCO2-anom box (yellow): ')
        print('lats '+str([lat_S, lat_N])+' and lons '+str([lon_E, lon_W])+' and depths '+str([depth_top, depth_btm]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)

    # Add green var-specific box (corresponding to this var's max anom)
    if var_depth_top is not None:
        ax.plot([var_lat_S, var_lat_N, var_lat_N, var_lat_S, var_lat_S], \
                [var_depth_btm, var_depth_btm, var_depth_top, var_depth_top, var_depth_btm], \
                 color=var_box_color, linewidth=box_lw) 
        # Also print out averages
        var_subset     = var.sel(depth=slice(var_depth_top, var_depth_btm), lat=slice(var_lat_S, var_lat_N), lon=slice(var_lon_E, var_lon_W))
        weights_subset = weights.sel(depth=slice(var_depth_top, var_depth_btm), lat=slice(var_lat_S, var_lat_N), lon=slice(var_lon_E, var_lon_W))
        
        print('Average over var-specific box (green): ')
        print('lats '+str([var_lat_S, var_lat_N])+' and lons '+str([var_lon_E, var_lon_W])+' and depths '+str([var_depth_top, var_depth_btm]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)
        # print('\n') 
           
    # Aesthetics
    ax.set_ylim(0, 5000); ax.invert_yaxis(); ax.set_facecolor('darkgrey')
    ax.set_ylabel('Depth [m]'); ax.set_xlabel('Latitude [$^\circ$N]')
    ax.set_title(var_descr+' at '+str(lon_of_interest)+'$^\circ$E', fontsize=13)



# %% Set lon_of_interest, fig width and height

lon_of_interest = 200   # degE
my_fig_width = 8; my_fig_height=4



# %% Take avg: PCO2 anom, depthxlat

NPac_lat_n = 60; NPac_lat_s = 30; NPac_depth_top = 500; NPac_depth_btm = 2000

# Print before averages printed
print('% - - - - - North Pacific Regional Averages - - - - - %')

PacBasin_fig(UV_pmoc.PCO2 - UV_ctrl.PCO2, 'PCO$_2$ Anomaly', my_fig_width, my_fig_height, lon_of_interest, -300, 300, \
        depth_top=NPac_depth_top, depth_btm=NPac_depth_btm, lat_S=NPac_lat_s, lat_N=NPac_lat_n, weights=dV_model)    
    
    
# %% Take avg: PO4 anom, depthxlat


# Print before averages printed
print('% - - - - - North Pacific Regional Averages - - - - - %')

PacBasin_fig(UV_pmoc.po4 - UV_ctrl.po4, 'PO$_4$ Anomaly', my_fig_width, my_fig_height, lon_of_interest, -1.75, 1.75, \
        depth_top=NPac_depth_top, depth_btm=NPac_depth_btm, lat_S=NPac_lat_s, lat_N=NPac_lat_n, weights=dV_model)    
    
    
    
# %% (4) Func for taking avgs across basin: 45degN, 30degS

def DepthLon_fig(var, var_descr, fig_width, fig_height, lat_of_interest, my_vmin=None, my_vmax=None, my_cmap='seismic', depth_top=None, depth_btm=None, lon_E=140, lon_W=240, box_color='yellow', box_lw=2, var_depth_top=None, var_depth_btm=None, var_lon_E=140, var_lon_W=240, var_box_color='lime', weights=None):
    """
    ! ! ! This is copied from section (4). Changes = 
     - var.interp(lon_of_interest)... --> var.interp(lat_of_interest)...
     - added: ax.set_xlim(lon_east, lon_west)
     - ax.set_xlabel('Latitude [$^\circ$N]') --> ax.set_xlabel('Longitude [$^\circ$E')
     - ax.set_title(var_descr+' at '+str(lon_of_interest)+'$^\circ$E', fontsize=13) -->
       ax.set_title(var_descr+' at '+str(lat_of_interest)+'$^\circ$N', fontsize=13)
     - drawing boxes: instead of lats & depths, --> lons & depths
     - Changed printing avgs: now avg at lat_of_interest, !!! weight by dxdz_model !!!
     - Removed inputs: lat_S, lat_N, var_lat_S, var_lat_N (don't need anymore; avg'ing AT a not OVER lat)
     
    This func will make one plot.
    - 'var' must have dims named 'latitude', 'longitude', and'depth'.
    - 'var_descr' is a string simply describing the variable 'var'.
    - 'weights' must be 3D of same dims as var. Will usually be volume of cells,
        i.e. dVol_model. !!! No, now dxdz_model !!!
    -  Note, matplotlib will automatically choose vmin/vmax if none are specified
        by user in calling the function.
    - It will also draw boxes, if depth_top/btm & lat_S/N are specified
    """
    fig, ax = plt.subplots(figsize=(fig_width, fig_height)) # 6.4, 4.8 default
    var.interp(lat=lat_of_interest).plot(vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, cbar_kwargs={'label': var.units})

    # Add yellow Delta_PCO2 box (box corresponding to max PCO2 anom, same box in all plots)
    if depth_top is not None:
        ax.plot([lon_W, lon_E, lon_E, lon_W, lon_W], \
                [depth_btm, depth_btm, depth_top, depth_top, depth_btm], \
                 color=box_color, linewidth=box_lw*2) 
        # Also print out averages
        var_subset     = var.interp(lat=lat_of_interest).sel(depth=slice(depth_top, depth_btm), lon=slice(lon_E, lon_W))
        weights_subset = weights.interp(lat=lat_of_interest).sel(depth=slice(depth_top, depth_btm), lon=slice(lon_E, lon_W))
        
        print('\n'+var_descr+' '+var.units)
        print('Average over PCO2-anom box (yellow) at '+str(lat_of_interest)+' degN: ')
        print('lons '+str([lon_E, lon_W])+' and depths '+str([depth_top, depth_btm]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)

    # Add green var-specific box (corresponding to this var's max anom)
    if var_depth_top is not None:
        ax.plot([var_lon_W, var_lon_E, var_lon_E, var_lon_W, var_lon_W], \
                [var_depth_btm, var_depth_btm, var_depth_top, var_depth_top, var_depth_btm], \
                 color=var_box_color, linewidth=box_lw) 
        # Also print out averages
        var_subset     = var.interp(lat=lat_of_interest).sel(depth=slice(var_depth_top, var_depth_btm), lon=slice(var_lon_E, var_lon_W))
        weights_subset = weights.interp(lat=lat_of_interest).sel(depth=slice(var_depth_top, var_depth_btm), lon=slice(var_lon_E, var_lon_W))
        
        print('Average over var-specific box (green) at '+str(lat_of_interest)+' degN: ')
        print('lons '+str([var_lon_E, var_lon_W])+' and depths '+str([var_depth_top, var_depth_btm]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)
        # print('\n') 
           
    # Aesthetics
    ax.set_xlim(lon_east, lon_west); ax.set_ylim(0, 5000); ax.invert_yaxis(); ax.set_facecolor('darkgrey')
    ax.set_ylabel('Depth [m]'); ax.set_xlabel('Longitude [$^\circ$E]')
    ax.set_title(var_descr+' at '+str(lat_of_interest)+'$^\circ$N', fontsize=13)



# %% FIG. 4: PCO2, PO4, and sigma2

my_fig_width  = 10.6   # 13
my_fig_height = 5.6    #  7

PCO2_anom = UV_pmoc.PCO2 - UV_ctrl.PCO2
PO4_anom  = UV_pmoc.po4  - UV_ctrl.po4

list_lats = [45, -30]
panel_labels = ['a', 'b', 'c', 'd']

cont_LW = 2.5  # contour line width and font size
cont_FS = 11
cont_accent = 'goldenrod'

# East-west cut_off of Pac basin [degE]:
lon_east = 140; lon_west = 280



fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(my_fig_width, my_fig_height))

for (row, col), ax in np.ndenumerate(axs):  
    # Note:
    # print(row*2 + col) # goes 0, 1, 2, 3; top left right, bottom left right
    
    if row == 0: 
        lat_to_use = list_lats[0]
    else: 
        lat_to_use = list_lats[1]
    
    if col == 0:
        var_to_use = PCO2_anom
        cbar_str = '[µatm]'
        my_cmap = 'seismic'
        if row == 0:
            my_vmin = -300
        else:
            my_vmin = -110
    else:
        var_to_use = PO4_anom
        cbar_str = '[µmol kg-1]'
        my_cmap = 'seismic'
        if row == 0:
            my_vmin = -2.2
        else:
            my_vmin = -0.8

    a = var_to_use.interp(lat=lat_to_use).plot(ax=ax, vmin=my_vmin, vmax=-my_vmin, cmap=my_cmap, cbar_kwargs={'label':cbar_str}) #, add_colorbar=False)


    red_levels = [36.0, 36.8]; yellow_levels = [36.5]
    
    # Add sigma2
    var1 = UV_pmoc.sigma2.interp(lat=lat_to_use)
    CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS, inline_spacing=-15)        
    CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS, inline_spacing=-15)         
                


    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_xlim((lon_east, lon_west)); ax.set_ylim((0,5500))
    ax.invert_yaxis()
    if row == 0:
        ax.set_xlabel('')
        if col == 0:
            ax.set_title('PCO$_2$ anomaly', fontweight='bold', fontsize = 18, y = 1.04)
        else:
            ax.set_title('PO$_4$$^{3-}$ anomaly', fontweight='bold', fontsize = 18, y = 1.03)
    else:
        ax.set_xlabel(r'Longitude [$\degree$E]')
        ax.set_title('')
    if col == 0:
        ax.set_ylabel('Depth [m]')
    else:
        ax.set_ylabel('')

    # Add 'a', 'b', 'c', 'd' annotations
    ax.text(-0.09, 1.05, panel_labels[row*2 + col], transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    # Add '45degN' or '30degS' labels
    if row == 0:
        ax.text(0.95, 0.13, r'45$\degree$N', transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', \
                     bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
    else:
        ax.text(0.95, 0.13, r'30$\degree$S', transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', \
                     bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))   
    

# Save fig
plt.savefig('../results/Fig4_PCO2po4Sigma2_45N30S.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig4_PCO2po4Sigma2_45N30S.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig4_PCO2po4Sigma2_45N30S.pdf', dpi=600, bbox_inches='tight')

plt.show()





# %% Take avg: PCO2 anom 30degSN

lat_of_interest = -30   # degN
S30_depth_top = 1000; S30_depth_btm = 2500
my_fig_width = 7
my_fig_height = 4

# Print before averages printed
print('% - - - - - Averages at '+str(lat_of_interest)+' latitude - - - - - %')


DepthLon_fig(UV_pmoc.PCO2 - UV_ctrl.PCO2, 'PCO$_2$ Anomaly', my_fig_width, my_fig_height, lat_of_interest, -100, 100, \
        depth_top=S30_depth_top, depth_btm=S30_depth_btm, lon_W=280, \
        var_depth_top=S30_depth_top, var_depth_btm=S30_depth_btm, var_lon_W=200, weights=dxdz_model)    

    
    
# %% Take avg: PO4 anom 30degSN

lat_of_interest = -30   # degN
S30_depth_top = 1000; S30_depth_btm = 2500
my_fig_width = 7
my_fig_height = 4

# Print before averages printed
print('% - - - - - Averages at '+str(lat_of_interest)+' latitude - - - - - %')


DepthLon_fig(UV_pmoc.po4 - UV_ctrl.po4, 'PO$_4$ Anomaly', my_fig_width, my_fig_height, lat_of_interest, -1, 1, \
        depth_top=S30_depth_top, depth_btm=S30_depth_btm, lon_W=280, \
        var_depth_top=S30_depth_top, var_depth_btm=S30_depth_btm, var_lon_W=200, weights=dxdz_model)    

    
    
    
# %% (5) Func for taking avgs across S.O. Surface 

def SOsurf_fig(var, var_descr, fig_width, fig_height, lat_cutoff, my_vmin=None, my_vmax=None, my_cmap='seismic', lat_N=None, lat_S=None, lon_W=0, lon_E=360, box_color='yellow', box_lw=2, var_lat_N=None, var_lat_S=None, var_lon_W=0, var_lon_E=360, var_box_color='lime', weights=None):
    """
    ! ! ! This is copied from section (5a). Changes = 
     - var.interp(lat_of_interest)... --> var.sel(depth=0, lat=slice(-90, lat_cutoff+1))...
     - changed title, ylabel. removed xlims and ax.invert_yaxis()
     - ylim(0,5000) --> ax.set_ylim(-80, lat_cutoff);
     - drawing boxes: instead of depths & lons --> lats and lons
     - printing avgs: var.interp(lat_of_interest).sel(depth=slice(), longitude=slice())
         --> var.sel(depth=0, method='nearest').sel(lat=slice(lat_S, lat_N), lon=slice(lon_W, lon_E))
     - ! ! ! weights will need to be a 3D array of dxdy areas ! ! ! at surf
     - default lons read in = 140:240/280 --> 0:360 (all the way around by default)
     
    This is copied from section (4). Changes = 
     - var.interp(lon_of_interest)... --> var.interp(lat_of_interest)...
     - added: ax.set_xlim(lon_east, lon_west)
     - ax.set_xlabel('Latitude [$^\circ$N]') --> ax.set_xlabel('Longitude [$^\circ$E')
     - ax.set_title(var_descr+' at '+str(lon_of_interest)+'$^\circ$E', fontsize=13) -->
       ax.set_title(var_descr+' at '+str(lat_of_interest)+'$^\circ$N', fontsize=13)
     - drawing boxes: instead of lats & depths, --> lons & depths
     - Changed printing avgs: now avg at lat_of_interest, !!! weight by dxdz_model !!!
     - Removed inputs: lat_S, lat_N, var_lat_S, var_lat_N (don't need anymore; avg'ing AT a not OVER lat)
     
    This func will make one plot.
    - 'var' must have dims named 'depth', 'lat', and'lon'.
    - 'var_descr' is a string simply describing the variable 'var'.
    - 'weights' must be 3D of same dims as var. Will usually be volume of cells,
        i.e. dVol_model. !!! No, now dxdz_model !!!
    -  Note, matplotlib will automatically choose vmin/vmax if none are specified
        by user in calling the function.
    - It will also draw boxes, if depth_top/btm & lat_S/N are specified
    """
    fig, ax = plt.subplots(figsize=(fig_width, fig_height)) # 6.4, 4.8 default
    var.sel(depth=0, method='nearest').sel(lat=slice(-90, lat_cutoff+1)).plot(vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, cbar_kwargs={'label': var.units})
    
    # Add yellow Delta_PCO2 box (box corresponding to max PCO2 anom, same box in all plots)
    if lat_N is not None:
        ax.plot([lon_W, lon_E, lon_E, lon_W, lon_W], \
                [lat_S, lat_S, lat_N, lat_N, lat_S], \
                 color=box_color, linewidth=box_lw*2) 
        # Also print out averages
        var_subset     = var.sel(depth=0, method='nearest').sel(lat=slice(lat_S, lat_N), lon=slice(lon_W, lon_E))
        weights_subset = weights.sel(lat=slice(lat_S, lat_N), lon=slice(lon_W, lon_E))
        
        print('\n'+var_descr+' '+var.units)
        print('Average over PCO2-anom box (yellow) at S.O. surf: ')
        print('lats '+str([lat_S, lat_N])+' and lons '+str([lon_W, lon_E]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)

    # Add green var-specific box (corresponding to this var's max anom)
    if var_lat_N is not None:
        ax.plot([var_lon_W, var_lon_E, var_lon_E, var_lon_W, var_lon_W], \
                [var_lat_S, var_lat_S, var_lat_N, var_lat_N, var_lat_S], \
                 color=var_box_color, linewidth=box_lw) 
        # Also print out averages
        var_subset     = var.sel(depth=0, method='nearest').sel(lat=slice(var_lat_S, var_lat_N), lon=slice(var_lon_W, var_lon_E))
        weights_subset = weights.sel(lat=slice(var_lat_S, var_lat_N), lon=slice(var_lon_W, var_lon_E))
        
        print('Average over var-specific box (green) at '+str(lat_of_interest)+' degN: ')
        print('lats '+str([var_lat_S, var_lat_N])+' and lons '+str([var_lon_W, var_lon_E]))
        # print(var_subset.mean().values)
        # print('This is the cell-vol-weighted mean:')
        print(var_subset.weighted(weights_subset).mean().values)
        # print('\n') 
           
    # Aesthetics
    ax.set_ylim(-80, lat_cutoff); ax.set_facecolor('darkgrey')
    ax.set_ylabel('Latitude [$^\circ$N]'); ax.set_xlabel('Longitude [$^\circ$E]')
    ax.set_title('Southern Ocean Surface '+var_descr, fontsize=13)

   

# %% Take avg: S.O. surf PCO2 anom

lat_cutoff = -30   # degN
my_fig_width = 10; my_fig_height=3.5

# Print before averages printed
print('% - - - - - Averages at S.O. Surface - - - - - %')

SOsurf_fig(UV_pmoc.PCO2 - UV_ctrl.PCO2, 'PCO$_2$ Anomaly', my_fig_width, my_fig_height, lat_cutoff, -30, 30, \
        lat_N=-40, lat_S=-70, var_lat_N = -47, var_lat_S=-70, var_lon_W=20, var_lon_E=290, weights=dxdy_model)    
      
  

# %% Take avg & fig: GLOBAL surf PCO2 anom
# New 2024-06-03

lat_cutoff = 60   # degN
my_fig_width = 10; my_fig_height=5 # 3.5

# Print before averages printed
print('% - - - - - Averages at S.O. Surface - - - - - %')

SOsurf_fig(UV_pmoc.PCO2 - UV_ctrl.PCO2, 'PCO$_2$ Anomaly', my_fig_width, my_fig_height, lat_cutoff, -30, 30, \
        lat_N=lat_cutoff, lat_S=-70, var_lat_N = lat_cutoff, var_lat_S=-70, var_lon_W=150, var_lon_E=280, weights=dxdy_model)    
          
    
anomX=( (UV_pmoc.PCO2 - UV_ctrl.PCO2).isel(depth=0) /  UV_ctrl.PCO2.isel(depth=0) ).mean()
print('\nanom % = '+str((anomX*100).values))  
anomXPac=((UV_pmoc.PCO2 - UV_ctrl.PCO2).sel(lon=slice(150,280)).isel(depth=0) / UV_ctrl.PCO2.sel(lon=slice(150,280)).isel(depth=0)).mean()
print('\nanom Pac % = '+str((anomXPac*100).values))       
          
    



# %% Take avg: S.O. surf PO4 anom

lat_cutoff = -30   # degN
my_fig_width = 10; my_fig_height=3.5

# Print before averages printed
print('% - - - - - Averages at S.O. Surface - - - - - %')

SOsurf_fig(UV_pmoc.po4 - UV_ctrl.po4, 'PO$_4$ Anomaly', my_fig_width, my_fig_height, lat_cutoff, -1, 1, \
        lat_N=-45, lat_S=-70, var_lat_N = -47, var_lat_S=-70, var_lon_W=20, var_lon_E=290, weights=dxdy_model)    
      


# %% Take avg & fig: GLOBAL surf PO4 anom
# New 2024-06-03

lat_cutoff = 60   # degN
my_fig_width = 10; my_fig_height=5 # 3.5

# Print before averages printed
print('% - - - - - Averages at S.O. Surface - - - - - %')

SOsurf_fig(UV_pmoc.po4 - UV_ctrl.po4, 'PO$_4$ Anomaly', my_fig_width, my_fig_height, lat_cutoff, -1, 1, \
        lat_N=lat_cutoff, lat_S=-70, var_lat_N = lat_cutoff, var_lat_S=-70, var_lon_W=150, var_lon_E=280, weights=dxdy_model)    
          
    

anomX=( (UV_pmoc.po4 - UV_ctrl.po4).isel(depth=0) /  UV_ctrl.po4.isel(depth=0) ).mean()
print('\nanom % = '+str((anomX*100).values))  
anomXPac=((UV_pmoc.po4 - UV_ctrl.po4).sel(lon=slice(150,280)).isel(depth=0) / UV_ctrl.po4.sel(lon=slice(150,280)).isel(depth=0)).mean()
print('\nanom Pac % = '+str((anomXPac*100).values))             

anomXSOPac=((UV_pmoc.po4 - UV_ctrl.po4).sel(lon=slice(150,280), lat=slice(-70,-45)).isel(depth=0) / UV_ctrl.po4.sel(lon=slice(150,280),lat=slice(-70,-45)).isel(depth=0)).mean()
print('\n\nanom SO Pac % = '+str((anomXSOPac*100).values))         
anomSO=((UV_pmoc.po4 - UV_ctrl.po4).sel(lat=slice(-70,-45)).isel(depth=0) / UV_ctrl.po4.sel(lat=slice(-70,-45)).isel(depth=0)).mean()
print('\nanom SO % = '+str((anomSO*100).values))         

# %% Largest PCO2 anomaly S.O.
(UV_pmoc.PCO2 - UV_ctrl.PCO2).sel(lat=slice(-90,-40)).sel(depth=0, method='nearest').min()
    
 
# %% Largest PO4 anomaly S.O.
(UV_pmoc.po4 - UV_ctrl.po4).sel(lat=slice(-90,-45)).sel(depth=0, method='nearest').min()
    


# %% (6) CALC CO2 FLUX by hand (to decompose) 

# %% Load in potT, Sal - for use later

ctrl_potT = UV_ctrl.potT
exp_potT  = UV_pmoc.potT

ctrl_sal = UV_ctrl.sal
exp_sal  = UV_pmoc.sal




# %% (6a) Sort winds - apply wind anomaly output to wind climatology following UVic code's procedure

# Model output
pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';

# - - - - - Load in windspeed anomalies - - - - - 

# eastward wind anomaly   # m s-1
ctrl_awindX = xr.open_dataset(pathname+'tavg.02001.01.01.nc', decode_times=False)['A_awindX'].isel(time=0)
exp_awindX  = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['A_awindX'].isel(time=-1)

# northward wind anomaly   # m s-1
ctrl_awindY = xr.open_dataset(pathname+'tavg.02001.01.01.nc', decode_times=False)['A_awindY'].isel(time=0)
exp_awindY  = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['A_awindY'].isel(time=-1)

# Rename coordinates
ctrl_awindX = ctrl_awindX.rename({'longitude_V':'longitude', 'latitude_V':'latitude'})
ctrl_awindY = ctrl_awindY.rename({'longitude_V':'longitude', 'latitude_V':'latitude'})
exp_awindX  = exp_awindX.rename({'longitude_V':'longitude', 'latitude_V':'latitude'})
exp_awindY  = exp_awindY.rename({'longitude_V':'longitude', 'latitude_V':'latitude'})
  # * longitude  (longitude) float64 3.6 7.2 10.8 14.4 ... 352.8 356.4 360.0
  # * latitude   (latitude) float64 -88.2 -86.4 -84.6 -82.8 ... 86.4 88.2 90.0


# - - - - - Load in windspeed [m/s] - - - - - - - 

pathname_ws = '../data/Menviel Wind Calc/'
# A_wind = xr.open_dataset(pathname_ws+'A_wind.nc')
A_windsur = xr.open_dataset(pathname_ws+'A_windsur.nc')
# Note different lat/lon values, and weird time values? [year]
  # * time       (time) float64 0.04247 0.1233 0.2041 ... 0.7904 0.874 0.9575
  # * longitude  (longitude) float32 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
  # * latitude   (latitude) float32 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1

model_windspeed = A_windsur.A_windspd.isel(time=-1) # .interp(longitude = ctrl_awindY.longitude, latitude = ctrl_awindY.latitude)
model_windangle = A_windsur.A_windang.isel(time=-1) # .interp(longitude = ctrl_awindY.longitude, latitude = ctrl_awindY.latitude)
  # * longitude  (longitude) float32 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
  # * latitude   (latitude) float32 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1
# [m s-1]
# [radians]


# - - - Interp anoms onto windspeed (match sea ice) - - - - - - 
ctrl_awindX = ctrl_awindX.interp(longitude = model_windspeed.longitude, latitude = model_windspeed.latitude)
ctrl_awindY = ctrl_awindY.interp(longitude = model_windspeed.longitude, latitude = model_windspeed.latitude)
exp_awindX  =  exp_awindX.interp(longitude = model_windspeed.longitude, latitude = model_windspeed.latitude)
exp_awindY  =  exp_awindY.interp(longitude = model_windspeed.longitude, latitude = model_windspeed.latitude)
  # * longitude  (longitude) float32 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
  # * latitude   (latitude) float32 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1


# - - - Add surface anom to windspeed - - - - - - 

p25 = 0.25  # doing some kind of average of 4 points
contr = 0.8
angle = 20.0 # degrees
cosa = np.cos(20.0*(np.pi/180))
sina = np.sin(20.0*(np.pi/180))

# arrays to be filled
final_windspeed_ctrl = np.empty(np.shape(ctrl_awindX))
final_windspeed_ctrl[:] = np.nan 
final_windspeed_exp = np.empty(np.shape(ctrl_awindX))
final_windspeed_exp[:] = np.nan 


# CTRL

# Loop thru 1:11
for ii in np.arange(1, np.size(ctrl_awindX.longitude)):
    for jj in np.arange(1, np.size(ctrl_awindX.latitude)):
        # print('CTRL ii = '+str(ii)+', jj = '+str(jj))
     
        component1 = ctrl_awindX.isel(longitude = ii,   latitude = jj)
        component2 = ctrl_awindX.isel(longitude = ii-1, latitude = jj)
        component3 = ctrl_awindX.isel(longitude = ii,   latitude = jj-1)
        component4 = ctrl_awindX.isel(longitude = ii-1, latitude = jj-1)
        ax = p25*(component1 + component2 + component3 + component4)

        component1y = ctrl_awindY.isel(longitude = ii,   latitude = jj)
        component2y = ctrl_awindY.isel(longitude = ii-1, latitude = jj)
        component3y = ctrl_awindY.isel(longitude = ii,   latitude = jj-1)
        component4y = ctrl_awindY.isel(longitude = ii-1, latitude = jj-1)
        ay = p25*(component1y + component2y + component3y + component4y)
        
        x = contr*(ax*cosa - ay*sina) + np.cos(model_windangle.isel(longitude=ii, latitude=jj))*model_windspeed.isel(longitude=ii, latitude=jj)
        y = contr*(ax*sina + ay*cosa) + np.sin(model_windangle.isel(longitude=ii, latitude=jj))*model_windspeed.isel(longitude=ii, latitude=jj)

        final_windspeed_ctrl[ii, jj] = np.sqrt(x**2 + y**2)
        
       
# EXP

# Loop thru 1:11
for ii in np.arange(1, np.size(exp_awindX.longitude)):
    for jj in np.arange(1, np.size(exp_awindX.latitude)):
        # print('EXP ii = '+str(ii)+', jj = '+str(jj))
        
        component1 = exp_awindX.isel(longitude = ii,   latitude = jj)
        component2 = exp_awindX.isel(longitude = ii-1, latitude = jj)
        component3 = exp_awindX.isel(longitude = ii,   latitude = jj-1)
        component4 = exp_awindX.isel(longitude = ii-1, latitude = jj-1)
        ax = p25*(component1 + component2 + component3 + component4)

        component1y = exp_awindY.isel(longitude = ii,   latitude = jj)
        component2y = exp_awindY.isel(longitude = ii-1, latitude = jj)
        component3y = exp_awindY.isel(longitude = ii,   latitude = jj-1)
        component4y = exp_awindY.isel(longitude = ii-1, latitude = jj-1)
        ay = p25*(component1y + component2y + component3y + component4y)
        
        x = contr*(ax*cosa - ay*sina) + np.cos(model_windangle.isel(longitude=ii, latitude=jj))*model_windspeed.isel(longitude=ii, latitude=jj)
        y = contr*(ax*sina + ay*cosa) + np.sin(model_windangle.isel(longitude=ii, latitude=jj))*model_windspeed.isel(longitude=ii, latitude=jj)

        final_windspeed_exp[ii, jj] = np.sqrt(x**2 + y**2)




# This will leave two empty rows at the top of each array of windspeed, just
# fill in with the values next to it so that grey bar doesn't appear on left/W 
# side of Fig. 5
final_windspeed_ctrl[0,:] = final_windspeed_ctrl[2,:]
final_windspeed_ctrl[1,:] = final_windspeed_ctrl[2,:]
final_windspeed_exp[0,:] = final_windspeed_exp[2,:]
final_windspeed_exp[1,:] = final_windspeed_exp[2,:]


# Rename coordinates
ctrl_awindY  =  ctrl_awindY.rename({'longitude':'lon', 'latitude':'lat'})

# Construct xr.DataArrays

final_windspeed_ctrl_da = xr.DataArray(
    data=np.transpose(final_windspeed_ctrl),
    dims=['lat', 'lon'],
    attrs=dict(
        units="m s1",
    ),
)    
final_windspeed_ctrl_da = final_windspeed_ctrl_da.assign_coords(lat=ctrl_awindY.lat, lon=ctrl_awindY.lon)


final_windspeed_exp_da = xr.DataArray(
    data=np.transpose(final_windspeed_exp),
    dims=['lat', 'lon'],
    attrs=dict(
        units="m s1",
    ),
)    
final_windspeed_exp_da = final_windspeed_exp_da.assign_coords(lat=ctrl_awindY.lat, lon=ctrl_awindY.lon)


  # * lat      (lat) float32 -89.1 -87.3 -85.5 -83.7 -81.9 ... 83.7 85.5 87.3 89.1
  # * lon      (lon) float32 1.8 5.4 9.0 12.6 16.2 ... 347.4 351.0 354.6 358.2


# %% (6b) Calc Schmidt number (Sc) and piston_vel (kw)
 
# - - - - - Load in exp sea ice fraction - - - - - 

exp_seaice  = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['O_icefra'].isel(time=-1)
ctrl_seaice  = xr.open_dataset(pathname+'tavg.02001.01.01.nc', decode_times=False)['O_icefra'].isel(time=0)
  # * longitude  (longitude) float64 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
  # * latitude   (latitude) float64 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1
# Rename coordinates
exp_seaice  =  exp_seaice.rename({'longitude':'lon', 'latitude':'lat'})
ctrl_seaice = ctrl_seaice.rename({'longitude':'lon', 'latitude':'lat'})


# - - - - - Sc number calculation - - - - - 

A = 2073.1; B = 125.62; C = 3.6276; D = 0.043219
temp_use = exp_potT.sel(depth=0, method='nearest') # degC is correct for eqn
Sc_exp = A - (B*temp_use) + (C*(temp_use**2)) - (D*(temp_use**3))
temp_use = ctrl_potT.sel(depth=0, method='nearest') # degC is correct for eqn
Sc_ctrl = A - (B*temp_use) + (C*(temp_use**2)) - (D*(temp_use**3))


# - - - - - kw calculations - - - - -  # I believe is cm/hr so convert --> m/s

xconv = (0.337/(60*60)) # this takes care of /hr --> /sec; cm to m taken care of below
# 0.337 a proportionality factor btw gas transfer velocity (cm/hr) and long-term
#  avg 10m windspeed squared; gives kw in cm/hr. See Wannikhof 1998 JGR-O.
# The /(60*60) is to further convert kw into sec-1 (so cm/sec.)
# Final answer will be m/s (see *(1/100) below)

# Piston velocity [m/s] in control and perturbed simulations
kw_exp =   (1-exp_seaice)*xconv* (final_windspeed_exp_da**2)*( (Sc_exp/660)**-0.5) * (1/100)
kw_ctrl = (1-ctrl_seaice)*xconv*(final_windspeed_ctrl_da**2)*((Sc_ctrl/660)**-0.5) * (1/100)

# For quantification purposes, hold all values at perturbation-values except 1,
#  to e.g. quantify the influence of sea ice alone on outgassing, or windspeed, 
#  or Schmidt number, etc
# All else == exp
kw_ctrlSI = (1-ctrl_seaice)*xconv*(final_windspeed_exp_da**2) *((Sc_exp/660)**-0.5) * (1/100)
kw_ctrlVs = (1-exp_seaice)* xconv*(final_windspeed_ctrl_da**2)*((Sc_exp/660)**-0.5) * (1/100)
kw_ctrlSc = (1-exp_seaice)* xconv*(final_windspeed_exp_da**2) *((Sc_ctrl/660)**-0.5) * (1/100)



# %% (6c) Calc Catm, from solubility (ff) from UVic code (roughly follows Wannikhof 1992)


# Calc dissolved CO2 (CO2*) in equilibrium with the atmosphere
#  NOTE: ctrl = 191 ppm, exp = 185 ppm.
#  NOTE: lower atmospheric ppm would actually act to !enhance! outgassing, all 
#         else held equal, and yet we still see reduced outgassing signal.
#  Ctrl pCO2 = 191e-6 atm
#  Exp pCO2 = 185e-6 atm

A1 = -162.8301; A2 = 218.2968; A3 = 90.9241; A4 = -1.47696; 
B1 = 0.025695; B2 = -0.025225; B3 = 0.0049867

tk100_exp  =  (exp_potT.sel(depth=0, method='nearest')+273.15)/100
tk100_ctrl = (ctrl_potT.sel(depth=0, method='nearest')+273.15)/100
# * longitude  (longitude) float64 1.8 5.4 9.0 12.6 ... 347.4 351.0 354.6 358.2
# * latitude   (latitude) float64 -89.1 -87.3 -85.5 -83.7 ... 85.5 87.3 89.1
# (and sal is the same. Units 1e-3 == permille which ff equation needs, good)

# According to Wannikhof 1992, gives ff in units of mol kg-1 atm-1. Do rough conversion to m-3 below (*1027 kg m-3)
ff_exp = np.exp( A1 + A2/tk100_exp  + A3*np.log(tk100_exp) + A4*(tk100_exp**2) \
                + exp_sal.sel(depth=0, method='nearest')*(B1 + B2*tk100_exp + B3*(tk100_exp**2) ) )

ff_ctrl = np.exp( A1 + A2/tk100_ctrl  + A3*np.log(tk100_ctrl) + A4*(tk100_ctrl**2) \
                + ctrl_sal.sel(depth=0, method='nearest')*(B1 + B2*tk100_ctrl + B3*(tk100_ctrl**2) ) )

    
#  gives 
# [mol m-3]  =  [mol kg-1 atm-1]  *  [atm]  * [approx density in kg m-3]
exp_Catm =       ff_exp * (185 * (10**-6)) * 1027

#  gives 
# [mol m-3]  =  [mol L-1 atm-1]  *  [atm]  * [approx density in kg m-3]
ctrl_Catm =    ff_ctrl * (191 * (10**-6)) * 1027


# Change in Catm is in part due to ppm change, and in part due to solubility (ff)
#  change. Decompose the total change into its two constituent parts (new Catm
#  keeping ppm constant and just changing ff, vs new Catm keeping ff/solubility
#  constant and just changing ppm)

# Testing impact of just ppm going 191>185 (use ctrl ppm and exp ff)
ctrl_Catm191 = ff_exp * (191 * (10**-6)) * 1027

# Testing impact of just solubility changing (use exp ppm and ctrl ff)
ctrl_CatmFF = ff_ctrl * (185 * (10**-6)) * 1027



# %% (6d) Calc Csw from dic and alk

input_P = 0 # dbar

# Set PyCO2SYS parameters (same as in 01_NPacSObgc_calcPCO2_wPO4_andSigma.py)
my_opt_pH_scale     = 2  # 2 = sw
my_opt_k_carbonic   = 5  # 5 = H73a, H73b and MCHP73 refit by DM87
my_opt_k_bisulfate  = 1  # 1 = D90a
my_opt_k_fluoride   = 1  # 1 = DR79
my_opt_total_borate = 2  # 2 = LKB10
my_opt_buffers_mode = 2 # how to calculate the various buffer factors (or not)
# 1: using automatic differentiation, which accounts for the effects of all equilibrating solutes (default).
# 2: using explicit equations reported in the literature, which only account for carbonate, borate and water alkalinity.
# 0: not at all.


# - - - - - exp_Csw - - - - -
input_alk  = UV_pmoc.alk.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 1)
input_dic  = UV_pmoc.dic.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 2)
input_sal  = UV_pmoc.sal.sel(depth=0, method='nearest')  # salinity 1e-3 (needs to be prac sal? ~35, ~good)
input_temp = UV_pmoc.potT.sel(depth=0, method='nearest')  # degC, good
input_po4  = UV_pmoc.po4.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys)

# μmol kg-1
exp_Csw = pyco2.sys(par1 = input_alk, par2 = input_dic, par1_type=1, par2_type=2, \
                          salinity = input_sal, temperature = input_temp, pressure = input_P, \
                              opt_pH_scale = my_opt_pH_scale, opt_k_carbonic = my_opt_k_carbonic, \
                                  opt_k_bisulfate = my_opt_k_bisulfate, opt_k_fluoride = my_opt_k_fluoride, \
                                      opt_total_borate = my_opt_total_borate, opt_buffers_mode=my_opt_buffers_mode, \
                                            total_phosphate = input_po4)['aqueous_CO2'];    
# Convert from to [μmol kg-1] to [mol m-3]   (Note, negligible diff (<0.5%) btw converting w/ ctrl-vs-exp rho)  
exp_Csw = exp_Csw * (10**-6) * UV_pmoc.rho.sel(depth=0, method='nearest')
exp_Csw = exp_Csw.assign_attrs(units="mol m-3", long_name='aqueous CO2', note='all exp values into pyco2sys')
del exp_Csw.attrs['standard_name']



# - - - - - ctrl_Csw - - - - -
input_alk  = UV_ctrl.alk.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 1)
input_dic  = UV_ctrl.dic.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 2)
input_sal  = UV_ctrl.sal.sel(depth=0, method='nearest')  # salinity 1e-3 (needs to be prac sal? ~35, ~good)
input_temp = UV_ctrl.potT.sel(depth=0, method='nearest') # degC, good
input_po4  = UV_ctrl.po4.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys)

# μmol kg-1
ctrl_Csw = pyco2.sys(par1 = input_alk, par2 = input_dic, par1_type=1, par2_type=2, \
                          salinity = input_sal, temperature = input_temp, pressure = input_P, \
                              opt_pH_scale = my_opt_pH_scale, opt_k_carbonic = my_opt_k_carbonic, \
                                  opt_k_bisulfate = my_opt_k_bisulfate, opt_k_fluoride = my_opt_k_fluoride, \
                                      opt_total_borate = my_opt_total_borate, opt_buffers_mode=my_opt_buffers_mode, \
                                            total_phosphate = input_po4)['aqueous_CO2'];    
# Convert from to [μmol kg-1] to [mol m-3]   (Note, negligible diff (<0.5%) btw converting w/ ctrl-vs-exp rho)  
ctrl_Csw = ctrl_Csw * (10**-6) * UV_ctrl.rho.sel(depth=0, method='nearest')
ctrl_Csw = ctrl_Csw.assign_attrs(units="mol m-3", long_name='aqueous CO2', note='all ctrl values into pyco2sys')
del ctrl_Csw.attrs['standard_name']



# To quantify the impact of T&S changes on outgassing, calc Csw with all exp values
# except ctrl_T&S, and compare the two calculations
# - - - - - ctrl_CswTS - - - - -
input_alk  =      UV_pmoc.alk.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 1)
input_dic  =      UV_pmoc.dic.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 2)
input_sal  = UV_ctrl.sal.sel(depth=0, method='nearest')  # salinity 1e-3 (needs to be prac sal? ~35, ~good)
input_temp = UV_ctrl.potT.sel(depth=0, method='nearest') # degC, good
input_po4  =      UV_pmoc.po4.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys)

# μmol kg-1
ctrl_CswTS = pyco2.sys(par1 = input_alk, par2 = input_dic, par1_type=1, par2_type=2, \
                          salinity = input_sal, temperature = input_temp, pressure = input_P, \
                              opt_pH_scale = my_opt_pH_scale, opt_k_carbonic = my_opt_k_carbonic, \
                                  opt_k_bisulfate = my_opt_k_bisulfate, opt_k_fluoride = my_opt_k_fluoride, \
                                      opt_total_borate = my_opt_total_borate, opt_buffers_mode=my_opt_buffers_mode, \
                                            total_phosphate = input_po4)['aqueous_CO2'];    
# Convert from to [μmol kg-1] to [mol m-3]   (Note, negligible diff (<0.5%) btw converting w/ ctrl-vs-exp rho)   
ctrl_CswTS = ctrl_CswTS * (10**-6) * UV_pmoc.rho.sel(depth=0, method='nearest')
ctrl_CswTS = ctrl_CswTS.assign_attrs(units="mol m-3", long_name='aqueous CO2', note='all exp values into pyco2sys except T,S')
del ctrl_CswTS.attrs['standard_name']



# To quantify the impact of Alk&DIC changes on outgassing, calc Csw with all exp values
# except ctrl_Alk&DIC, and compare the two calculations
# - - - - - ctrl_CswAlkDIC - - - - -
input_alk  = UV_ctrl.alk.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 1)
input_dic  = UV_ctrl.dic.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys) (par type 2)
input_sal  =      UV_pmoc.sal.sel(depth=0, method='nearest')  # salinity 1e-3 (needs to be prac sal? ~35, ~good)
input_temp =      UV_pmoc.potT.sel(depth=0, method='nearest')  # degC, good
input_po4  =      UV_pmoc.po4.sel(depth=0, method='nearest')  # umol/kg (correct for pyCO2sys)

# μmol kg-1
ctrl_CswAlkDIC = pyco2.sys(par1 = input_alk, par2 = input_dic, par1_type=1, par2_type=2, \
                          salinity = input_sal, temperature = input_temp, pressure = input_P, \
                              opt_pH_scale = my_opt_pH_scale, opt_k_carbonic = my_opt_k_carbonic, \
                                  opt_k_bisulfate = my_opt_k_bisulfate, opt_k_fluoride = my_opt_k_fluoride, \
                                      opt_total_borate = my_opt_total_borate, opt_buffers_mode=my_opt_buffers_mode, \
                                            total_phosphate = input_po4)['aqueous_CO2'];    
# Convert from to [μmol kg-1] to [mol m-3]   (Note, negligible diff (<0.5%) btw converting w/ ctrl-vs-exp rho)   
ctrl_CswAlkDIC = ctrl_CswAlkDIC * (10**-6) * UV_pmoc.rho.sel(depth=0, method='nearest')
ctrl_CswAlkDIC = ctrl_CswAlkDIC.assign_attrs(units="mol m-3", long_name='aqueous CO2', note='all exp values into pyco2sys except alk,DIC,po4')
del ctrl_CswAlkDIC.attrs['standard_name']



# %% (6e) Calc fluxes and flux anoms

# f = -1*k*[Catm-Csw]       # -1 so that (+)f = outgassing (Csw > Catm == outgassing)
#  units: m/s * mol/m3 == mol m-2 s-1. Good.

f_exp = -1*kw_exp*(exp_Catm - exp_Csw)

f_ctrl = -1*kw_ctrl*(ctrl_Catm - ctrl_Csw)

f_ctrlIce = -1*kw_ctrlSI*(exp_Catm - exp_Csw)
f_ctrlVs  = -1*kw_ctrlVs*(exp_Catm - exp_Csw) 
f_ctrlSc  = -1*kw_ctrlSc*(exp_Catm - exp_Csw) 

f_ctrlCatm    = -1*kw_exp*(ctrl_Catm - exp_Csw)
f_ctrlCatm191 = -1*kw_exp*(ctrl_Catm191 - exp_Csw)
f_ctrlCatmff  = -1*kw_exp*(ctrl_CatmFF - exp_Csw)

f_ctrlSW     = -1*kw_exp*(exp_Catm - ctrl_Csw)
f_ctrlTS     = -1*kw_exp*(exp_Catm - ctrl_CswTS)
f_ctrlAlkDIC = -1*kw_exp*(exp_Catm - ctrl_CswAlkDIC)

# Later convert from mol m-2 s-1 to mol m-2 yr-1
# re-do long_name and units, del standard_name, add note, 

# Give proper units and description

f_exp.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_exp.attrs['units'] = 'mol m-2 s-1'
del f_exp.attrs['standard_name']

f_ctrl.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrl.attrs['units'] = 'mol m-2 s-1'
del f_ctrl.attrs['standard_name']

f_ctrlIce.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlIce.attrs['units'] = 'mol m-2 s-1'
del f_ctrlIce.attrs['standard_name']

f_ctrlVs.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlVs.attrs['units'] = 'mol m-2 s-1'
del f_ctrlVs.attrs['standard_name']

f_ctrlSc.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlSc.attrs['units'] = 'mol m-2 s-1'
del f_ctrlSc.attrs['standard_name']

f_ctrlCatm.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlCatm.attrs['units'] = 'mol m-2 s-1'
del f_ctrlCatm.attrs['standard_name']

f_ctrlCatm191.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlCatm191.attrs['units'] = 'mol m-2 s-1'
del f_ctrlCatm191.attrs['standard_name']

f_ctrlCatmff.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlCatmff.attrs['units'] = 'mol m-2 s-1'
del f_ctrlCatmff.attrs['standard_name']

f_ctrlSW.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlSW.attrs['units'] = 'mol m-2 s-1'
del f_ctrlSW.attrs['standard_name']

f_ctrlTS.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlTS.attrs['units'] = 'mol m-2 s-1'
del f_ctrlTS.attrs['standard_name']

f_ctrlAlkDIC.attrs['long_name'] = 'CO2 flux from ocean to atms'
f_ctrlAlkDIC.attrs['units'] = 'mol m-2 s-1'
del f_ctrlAlkDIC.attrs['standard_name']



f_exp.attrs['description'] = 'all exp values'

f_ctrl.attrs['description'] = 'all ctrl values'

f_ctrlIce.attrs['description'] = 'all exp values except ctrl sea ice fraction'
f_ctrlVs.attrs['description'] = 'all exp values except ctrl windspeeds'
f_ctrlSc.attrs['description'] = 'all exp values except ctrl Sc number'

f_ctrlCatm.attrs['description'] = 'all exp values except ctrl Catm (f(ppm, solubility))'
f_ctrlCatm191.attrs['description'] = 'all exp values except ctrl ppm (and exp ff/solubility)'
f_ctrlCatmff.attrs['description'] = 'all exp values  except ctrl ff/solubility (and exp ppm)'

f_ctrlSW.attrs['description'] = 'all exp values except ctrl Csw (f(T,S,Alk,DIC,po4))'
f_ctrlTS.attrs['description'] = 'all exp values except ctrl T&S into Csw calc'
f_ctrlAlkDIC.attrs['description'] = 'all exp values except ctrl Alk&DIC into Csw calc'



# %% (6f) Calc flux yr-1 anom signals due to various factors


f_anom = (f_exp - f_ctrl)*(365*24*60*60)

# BGC
f_anom_dueToAlkDIC = (f_exp - f_ctrlAlkDIC)*(365*24*60*60)
f_anom_dueToTS     = (f_exp - f_ctrlTS)*(365*24*60*60)
# *

# Atms/Kh
f_anom_dueToCatm = (f_exp - f_ctrlCatm)*(365*24*60*60)
f_anom_dueToCatm191 = (f_exp - f_ctrlCatm191)*(365*24*60*60)
f_anom_duetoCatmKh = (f_exp - f_ctrlCatmff)*(365*24*60*60)


# Sea ice, winds
f_anom_dueToSeaIce = (f_exp - f_ctrlIce)*(365*24*60*60) # v. small
f_anom_dueToWinds  = (f_exp - f_ctrlVs)*(365*24*60*60)  # extremely small, order of magnitude
f_anom_dueToSc     = (f_exp - f_ctrlSc)*(365*24*60*60)  # v. small


# * 
f_anom_dueToAllBGC = (f_exp - f_ctrlSW)*(365*24*60*60)



f_exp = f_exp*(365*24*60*60)
f_ctrl = f_ctrl*(365*24*60*60)

f_exp = f_exp.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux', note='Ventilated')
f_ctrl = f_ctrl.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux', note='Ctrl')


f_anom = f_anom.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')

f_anom_dueToAlkDIC = f_anom_dueToAlkDIC.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')
f_anom_dueToTS = f_anom_dueToTS.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')

f_anom_dueToCatm = f_anom_dueToCatm.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')
f_anom_dueToCatm191 = f_anom_dueToCatm191.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')
f_anom_duetoCatmKh = f_anom_duetoCatmKh.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')

f_anom_dueToSeaIce = f_anom_dueToSeaIce.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')
f_anom_dueToWinds = f_anom_dueToWinds.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')
f_anom_dueToSc = f_anom_dueToSc.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')

f_anom_dueToAllBGC = f_anom_dueToAllBGC.assign_attrs(units="mol m-2 yr-1", long_name='upwards carbon flux anomaly', note='Exp - Ctrl')



del f_anom.attrs['description']

del f_anom_dueToAlkDIC.attrs['description']
del f_anom_dueToTS.attrs['description']

del f_anom_dueToCatm.attrs['description']
del f_anom_dueToCatm191.attrs['description']
del f_anom_duetoCatmKh.attrs['description']

del f_anom_dueToSeaIce.attrs['description']
del f_anom_dueToWinds.attrs['description']
del f_anom_dueToSc.attrs['description']

del f_anom_dueToAllBGC.attrs['description']




# Save fluxes for use later

# Don't need these dims
f_exp = f_exp.drop(['time', 'depth'])
f_ctrl = f_ctrl.drop(['time', 'depth'])

f_exp_ds  = f_exp.to_dataset(name='flux')
f_ctrl_ds = f_ctrl.to_dataset(name='flux')

f_exp_ds.to_netcdf('../results/UV_pmoc_flux.nc')
f_ctrl_ds.to_netcdf('../results/UV_ctrl_flux.nc')





# %% (6g) Func to plot flux anoms

def plot_flux_anoms(data, title_str, my_vmin=None, my_vmax=None):
    
    fig, ax = plt.subplots(figsize=(10,3.5)) # S.O. plot
    
    data.sel(lat=slice(-90, -29.700001)).plot(ax=ax, cmap='seismic', vmin=my_vmin, vmax=my_vmax, cbar_kwargs={'label': data.units})
    
    # Aesthetics
    ax.set_facecolor('darkgrey')
    ax.set_ylabel('Latitude ($^\circ$N)'); ax.set_xlabel('Longitude ($^\circ$E)')
    ax.set_title(title_str)
    ax.set_ylim(-80, -30)




# %% (6h) Visuals - Outgas anoms mol m-2 yr-1

vmin_forplot = -10/3; vmax_forplot = -1*vmin_forplot # mol m-2 yr-1


plot_flux_anoms(f_anom, 'Outgassing Anomaly (Exp-Ctrl)', vmin_forplot, vmax_forplot) # -6, 6 also works
plot_flux_anoms(f_anom_dueToAllBGC, 'Outgassing Anomaly due to ALL Ocn BGC Change (T,S,Alk,DIC,PO4)', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_dueToAlkDIC, 'Outgassing Anomaly due to Alk+DIC Changes ONLY', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_dueToTS, 'Outgassing Anomaly due to T+S Changes ONLY', vmin_forplot, vmax_forplot)

plot_flux_anoms(f_anom_dueToCatm, 'Outgassing Anomaly due to Change in $C_{atm}$ (191-->185ppm)', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_dueToCatm191, 'Outgassing Anomaly due to only Atms Atms 191>185ppm', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_duetoCatmKh, 'Outgassing Anomaly due to Kh exp vs ctrl', vmin_forplot, vmax_forplot)

plot_flux_anoms(f_anom_dueToSeaIce, 'Outgassing Anomaly due to Sea Ice Changes', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_dueToWinds, 'Outgassing Anomaly due to Wind Changes', vmin_forplot, vmax_forplot)
plot_flux_anoms(f_anom_dueToSc, 'Outgassing Anomaly due to Sc Changes', vmin_forplot, vmax_forplot)







# %% (7) FIG. 5 w sigma2: 4-panel S.O. (w PO4 now): po4, PCO2, outgas, and outgas-due-to-DIC&Alk

lat_cutoff = -30

my_fig_width  = 9
my_fig_height = 14  # 11   # 10
my_cmap = 'seismic'
panel_labels = ['a', 'b', 'c', 'd']

cont_LW = 1.5
cont_FS = 11
cont_accent = 'goldenrod'


# Vars: po4_anom, PCO2 anom, outgassing anom, and outgassing anom due only to DIC&Alk changes
PCO2_anom = (UV_pmoc.PCO2 - UV_ctrl.PCO2).sel(depth=0, method='nearest')
po4_anom =  (UV_pmoc.po4  - UV_ctrl.po4).sel(depth=0, method='nearest')

list_vars = [po4_anom, PCO2_anom, f_anom, f_anom_dueToAlkDIC]



fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(my_fig_width, my_fig_height))

for ii in np.arange(np.size(axs)):
    
    # Panel-specific details
    if ii == 0:     # PCO2 anom
        my_vmin = -0.7; my_vmax = -my_vmin
        title_str = r'PO$_4^{3-}$ Anomaly'; cbar_str = '[µmol kg-1]'
    elif ii == 1:   # PCO2 anomaly
        my_vmin = -30; my_vmax = -my_vmin
        title_str = r'PCO$_2$ Anomaly'; cbar_str = '[µatm]'    
    elif ii == 2:   # outgassing anomaly
        my_vmin = -3; my_vmax = -my_vmin
        title_str = r'Air-Sea CO$_2$ Flux Anomaly'; cbar_str = '[mol m-2 yr-1]'
    elif ii == 3:   # outgassing anomaly due only to DIC & Alk changes
        my_vmin = -3; my_max = -my_vmin
        title_str = r'CO$_2$ Flux Anomaly due ONLY to changes in DIC & Alk.'
        cbar_str = '[mol m-2 yr-1]'

    
    var = list_vars[ii]
    a = var.sel(lat=slice(-90, lat_cutoff+1)).plot(ax=axs[ii], vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, cbar_kwargs={'label':cbar_str}) # , add_colorbar=False)
    
    
    # # Customise colorbar
    # if ii == 0:
    #     # Create a new axis for colorbar
    #     cbar_ax = fig.add_axes([0.93, 0.663, 0.025, 0.215])
    # elif ii == 1:
    #     cbar_ax = fig.add_axes([0.93, 0.394, 0.025, 0.215])
    # else:
    #     cbar_ax = fig.add_axes([0.93, 0.124, 0.025, 0.215]) 
    # cbar = plt.colorbar(a, cax=cbar_ax)
    # cbar.set_label(label=cbar_str, fontsize=12, rotation=-90, va='bottom', labelpad=1)
    # cbar.ax.tick_params(labelsize=12)   
    
    # Add sigma2
    red_levels = [36.0, 36.8]; yellow_levels = [36.5]

    var1 = UV_pmoc.sigma2.interp(depth=500)
    CS = var1.plot.contour(ax=axs[ii], colors='red', linewidths=cont_LW, levels=red_levels); axs[ii].clabel(CS, fontsize=cont_FS)        
    CS = var1.plot.contour(ax=axs[ii], colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); axs[ii].clabel(CS, fontsize=cont_FS)         
                        
    
    
    # Aesthetics
    axs[ii].set_ylim(-80, lat_cutoff); axs[ii].set_facecolor('darkgrey')
    axs[ii].set_ylabel('Latitude [$^\circ$N]')
    if ii == 3:
        axs[ii].set_xlabel('Longitude [$^\circ$E]')
    else:
        axs[ii].set_xlabel('')
    axs[ii].set_title(title_str, fontsize=13, fontweight='bold')
    
    # Add a, b,annotations
    axs[ii].text(-0.05, 1.13, panel_labels[ii], transform=axs[ii].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    
# Adjust panels
fig.subplots_adjust(hspace=0.25)
    
# Save fig
plt.savefig('../results/Fig5_po4_PCO2_outgas_anoms_sigma2.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig5_po4_PCO2_outgas_anoms_sigma2.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/Fig5_po4_PCO2_outgas_anoms_sigma2.pdf', dpi=600, bbox_inches='tight')

plt.show()


    
    

# %% (8) SUPP. FIG. S.O. Flux Decomp w/ sigma2


my_fig_width  = 14.4   # 16
my_fig_height =  9.9   # 11   # 10
my_cmap = 'seismic'
my_vmin = -3; my_vmax = -my_vmin

cont_LW = 1.5
cont_FS = 11
cont_accent = 'goldenrod'


# The loop below will go left,right (top row), then left,right (mid), left,right (btm)
# So be careful about list of variables and a,b,c,d,e,f
list_vars = [f_anom, f_anom_dueToCatm, f_anom_dueToAlkDIC, f_anom_dueToSeaIce, f_anom_dueToTS, f_anom_dueToWinds]
panel_labels = ['a', 'd', 'b', 'e', 'c', 'f']
panel_titles = ['Air-Sea CO$_2$ Flux Anomaly', 'Flux Anomaly due to Change in C$_{atm}$', \
                'Flux Anomaly due to Changes in DIC & Alkalinity', 'Flux Anomaly due to Changes in Sea Ice Fraction', \
                'Flux Anomaly due to Changes in Temperature & Salinity', 'Flux Anomaly due to Changes in Windspeed']

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(my_fig_width, my_fig_height))

# Loop through models (Note, (row*2 + col) goes [0, 1, 2, 3, 4, 5, 6])
for (row, col), ax in np.ndenumerate(axs):  
    
    # Plot PCO2
    var = list_vars[row*2 + col]    

    a = var.sel(lat=slice(-90, lat_cutoff+1)).plot(ax=ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, add_colorbar=False)
    
    
    # Add sigma2
    red_levels = [36.0, 36.8]; yellow_levels = [36.5]

    var1 = UV_pmoc.sigma2.interp(depth=500)
    CS = var1.plot.contour(ax=ax, colors='red', linewidths=cont_LW, levels=red_levels); ax.clabel(CS, fontsize=cont_FS)        
    CS = var1.plot.contour(ax=ax, colors=cont_accent, linewidths=cont_LW, levels=yellow_levels); ax.clabel(CS, fontsize=cont_FS)         
                     
    
    # Customise colorbar
    
    # Aesthetics
    ax.set_ylim(-80, lat_cutoff); ax.set_facecolor('darkgrey')
    ax.set_title(panel_titles[row*2+col], fontweight='bold')

    if (row*2+col == 4) or (row*2+col == 5):
        ax.set_xlabel('Longitude [$^\circ$E]')
    else:
        ax.set_xlabel('')

    if (row*2+col == 0) or (row*2+col == 2) or (row*2+col == 4):
        ax.set_ylabel('Latitude [$^\circ$N]')
    else:
        ax.set_ylabel('')
    
    
    # Add a, b,annotations
    ax.text(-0.05, 1.11, panel_labels[row*2 + col], transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# Create a new axis for colorbar
cbar_ax = fig.add_axes([0.92, 0.124, 0.015, 0.756])  # Adjust the position and size as needed

# Add common colorbar
cbar = plt.colorbar(a, cax=cbar_ax)
cbar.set_label(label='[mol m$^{-2}$ yr$^{-1}$]', fontsize=14, weight='bold', rotation=-90, va='bottom', labelpad=1)
cbar.ax.tick_params(labelsize=14)

# Adjust panels
fig.subplots_adjust(hspace=0.25, wspace=0.12)
    
# Save fig
plt.savefig('../results/SuppFig_OutgasAnom_Decomposed_sigma2.eps', format='eps', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_OutgasAnom_Decomposed_sigma2.png', dpi=600, bbox_inches='tight')
plt.savefig('../results/SuppFig_OutgasAnom_Decomposed_sigma2.pdf', dpi=600, bbox_inches='tight')

plt.show()





# %% - - - - - - - - - - - - - - - - - - - - - - - - 

# %% (9) Putting some numbers on outgassing reductions


# %% - - - - - - -


# %% Write func for surface area-wt'd averages

def print_SO_areawtd_avgs(data_descr, ctrl_data, exp_data, anom_data, weights, cutoff_lat, east_lon, west_lon, S_cutoff_lat=None):
    """
    Prints out area-weighted averages across whole S.O. (90S up to cutoff_lat)
    and just Pac sector S.O. (90S to cutoff_lat, across west_lon to east_lon).
    Fed in data must already be at depth (surface) desired.
    
    data_descr : string describing the data
    ctrl_data, exp_data, anom_data : must have dims "latitude" and "longitude",
      must also have an attribute called "units"
    weights: must have dim "latitude"
    cutoff_lat : latitude below which S.O. starts (degN)
    east_lon, west_lon: edges of Pac sector of S.O.
    
    """

    # Calc area-wt'd averages
    
    # For S.O. winds:
    if S_cutoff_lat is not None:
        ctrl_val = ctrl_data.sel(lat=slice(S_cutoff_lat, cutoff_lat)).weighted(weights).mean().values
        exp_val  = exp_data.sel(lat=slice(S_cutoff_lat, cutoff_lat)).weighted(weights).mean().values
        anom_val = anom_data.sel(lat=slice(S_cutoff_lat, cutoff_lat)).weighted(weights).mean().values
        
        Pac_ctrl_val = ctrl_data.sel(lat=slice(S_cutoff_lat, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values
        Pac_exp_val  = exp_data.sel(lat=slice(S_cutoff_lat, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values
        Pac_anom_val = anom_data.sel(lat=slice(S_cutoff_lat, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values    
    else:    
        ctrl_val = ctrl_data.sel(lat=slice(-90, cutoff_lat)).weighted(weights).mean().values
        exp_val  = exp_data.sel(lat=slice(-90, cutoff_lat)).weighted(weights).mean().values
        anom_val = anom_data.sel(lat=slice(-90, cutoff_lat)).weighted(weights).mean().values
        
        Pac_ctrl_val = ctrl_data.sel(lat=slice(-90, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values
        Pac_exp_val  = exp_data.sel(lat=slice(-90, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values
        Pac_anom_val = anom_data.sel(lat=slice(-90, cutoff_lat), lon=slice(east_lon, west_lon)).weighted(weights).mean().values
    
    # Print them out
    
    print('\n\n\nUVic mean '+data_descr+' ['+anom_data.units+']')
    
    print('South of '+str(-1*cutoff_lat)+'S:')
    print(' CTRL  = '+str(ctrl_val))
    print(' EXP   = '+str(exp_val))
    print(' ANOM  = '+str(anom_val)+' (avg of anomaly field)')
    print(' ANOM val = '+'{:.1f}'.format(anom_val/ctrl_val*100)+' % of CTRL val')
    
    print('\nS.O. Pac sector ('+str(east_lon)+' degE - '+str(360-west_lon)+' degW):')
    print(' CTRL  = '+str(Pac_ctrl_val))
    print(' EXP   = '+str(Pac_exp_val))
    print(' ANOM  = '+str(Pac_anom_val)+' (avg of anomaly field)')
    print(' ANOM val = '+'{:.1f}'.format(Pac_anom_val/Pac_ctrl_val*100)+' % of CTRL val')




# %% Get cell areas, and weights from cos(lat)

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
all_data = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)

# G_areaT = tracer grid cell area; same coords as flux.
grid_areas_T = all_data.G_areaT
# Rename coordinates
grid_areas_T = grid_areas_T.rename({'longitude':'lon', 'latitude':'lat'})


weights = np.cos(np.deg2rad(all_data.G_areaT.latitude))
weights.name = "weights"
del weights.attrs['units'] 
del weights.attrs['edges'] 
del weights.attrs['standard_name'] 
weights.attrs['long_name'] = 'weights'




# %% Print flux/Outgassing anoms

# East-west cut_off of INDO-Pac basin [degE]:
east_lon = 25; west_lon = 280

cutoff_lat = -40

# Print out area-wt'd averages (careful, plain old f_ctrl and f_exp not converted to yr-1 yet)
print_SO_areawtd_avgs('fDIC into atms', f_ctrl, f_exp, f_anom, weights, cutoff_lat, east_lon, west_lon)


