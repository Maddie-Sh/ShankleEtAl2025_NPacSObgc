# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 18:06:25 2023

@author: mgs23
"""

# %% (0) Set up

# Import modules

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


# Set parameters

depth_for_strat = 500   # (500m - surf: greater number = more stratified)





# %% (1) Read in output: rho, wind stress


# - - - - Wind Stress - - - - - - #

# Perturbed simulation ("UVic-NP" in text)

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
x = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)
# time=10. array([100.5, 196., 296., 396., 496., 596., 696., 796., 896., 996.])
tauX_UV_pmoc  = x.O_tauX.isel(time=-1)
tauY_UV_pmoc  = x.O_tauY.isel(time=-1)

# Control simulation ("UVic-ctrl" in text)

x = xr.open_dataset(pathname+'tavg.02001.01.01.nc', decode_times=False)
# time=1. array([2000.5]) (the end of the LGM experiment)
tauX_UV_ctrl  = x.O_tauX.isel(time=0)
tauY_UV_ctrl  = x.O_tauY.isel(time=0)


# - - - - rho - - - - - - - - - - #

pathname = '../results/';
rho_UV_pmoc = xr.open_dataset(pathname+'UV_pmoc.nc', decode_times=False)['rho']
rho_UV_ctrl = xr.open_dataset(pathname+'UV_ctrl.nc', decode_times=False)['rho']





# %% Wind stress: rename dims to depth, lat, lon

tauX_UV_ctrl = tauX_UV_ctrl.rename({'longitude_V':'lon', 'latitude_V':'lat'})
tauX_UV_pmoc = tauX_UV_pmoc.rename({'longitude_V':'lon', 'latitude_V':'lat'})
tauY_UV_ctrl = tauY_UV_ctrl.rename({'longitude_V':'lon', 'latitude_V':'lat'})
tauY_UV_pmoc = tauY_UV_pmoc.rename({'longitude_V':'lon', 'latitude_V':'lat'})




# %% Load dA (dx*dy) for each model

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
dA_forTau = xr.open_dataset(pathname+'tavg_tot.nc',decode_times=False)['G_areaU'] # velocity grid area [m2]
dA_forRho = xr.open_dataset(pathname+'tavg_tot.nc',decode_times=False)['G_areaT'] # velocity grid area [m2]

# Rename dims to lat, lon
dA_forTau = dA_forTau.rename({'longitude_V':'lon', 'latitude_V':'lat'})
dA_forRho = dA_forRho.rename({'longitude':'lon', 'latitude':'lat'})

# The same for perturbed & ctrl ('tavg_tot.nc' & 'tavg.02001.01.01.nc')





# %% --------- (2) STRATIFICATION -------------


# %% Calc strat (500m-surf)

# depth_for_strat = 500     # set above
                                                                                       # = 17.5 m 
strat_UV_ctrl = rho_UV_ctrl.interp(depth=depth_for_strat) - rho_UV_ctrl.sel(depth=0, method='nearest')
strat_UV_pmoc = rho_UV_pmoc.interp(depth=depth_for_strat) - rho_UV_pmoc.sel(depth=0, method='nearest')
strat_UV_anom = strat_UV_pmoc - strat_UV_ctrl

                                                                                       

# %% Calc area-wt'd mean strat & anomaly south of 40S

lower_lat=-90; upper_lat=-40
Pac_west_lim = 150; Pac_east_lim=280


# First apply weights

wt_strat_ctrl = (strat_UV_ctrl).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)
wt_strat_pmoc = (strat_UV_pmoc).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)
wt_strat_anom = (strat_UV_anom).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)

wt_strat_ctrl_Pac = (strat_UV_ctrl).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_strat_pmoc_Pac = (strat_UV_pmoc).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_strat_anom_Pac = (strat_UV_anom).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)


# Then take avg

meanStratSO_UV14_ctrl_AreaWtd = wt_strat_ctrl.mean()
meanStratSO_UV14_pmoc_AreaWtd = wt_strat_pmoc.mean()
meanStratSO_UV14_anom_AreaWtd = wt_strat_anom.mean()

meanStratPacSO_UV14_ctrl_AreaWtd = wt_strat_ctrl_Pac.mean()
meanStratPacSO_UV14_pmoc_AreaWtd = wt_strat_pmoc_Pac.mean()
meanStratPacSO_UV14_anom_AreaWtd = wt_strat_anom_Pac.mean()



# Print out all the numbers
print('Mean Stratification (Grad-rho, 500m - surf), Full S.O.:')
print("UVIC '14")
print('  meanStratSO_UV14_ctrl_AreaWtd = '+str(meanStratSO_UV14_ctrl_AreaWtd.values))
print('  meanStratSO_UV14_pmoc_AreaWtd = '+str(meanStratSO_UV14_pmoc_AreaWtd.values))
print('  meanStratSO_UV14_anom_AreaWtd = '+str(meanStratSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanStratSO_UV14_anom_AreaWtd/meanStratSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)


print('\nMean Stratification (Grad-rho, 500m - surf), Pac-Sector S.O.:')
print("UVIC '14")
print('  meanStratPacSO_UV14_ctrl_AreaWtd = '+str(meanStratPacSO_UV14_ctrl_AreaWtd.values))
print('  meanStratPacSO_UV14_pmoc_AreaWtd = '+str(meanStratPacSO_UV14_pmoc_AreaWtd.values))
print('  meanStratPacSO_UV14_anom_AreaWtd = '+str(meanStratPacSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanStratPacSO_UV14_anom_AreaWtd/meanStratPacSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)


# %% Calc area-wt'd mean strat & anomaly - GLOBALLY 

lower_lat=-90; upper_lat=90
Pac_west_lim = 150; Pac_east_lim=280


# First apply weights

wt_strat_ctrl = (strat_UV_ctrl).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)
wt_strat_pmoc = (strat_UV_pmoc).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)
wt_strat_anom = (strat_UV_anom).sel(lat=slice(lower_lat, upper_lat)).weighted(dA_forRho)

wt_strat_ctrl_Pac = (strat_UV_ctrl).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_strat_pmoc_Pac = (strat_UV_pmoc).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_strat_anom_Pac = (strat_UV_anom).sel(lat=slice(lower_lat,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)


# Then take avg

meanStratSO_UV14_ctrl_AreaWtd = wt_strat_ctrl.mean()
meanStratSO_UV14_pmoc_AreaWtd = wt_strat_pmoc.mean()
meanStratSO_UV14_anom_AreaWtd = wt_strat_anom.mean()

meanStratPacSO_UV14_ctrl_AreaWtd = wt_strat_ctrl_Pac.mean()
meanStratPacSO_UV14_pmoc_AreaWtd = wt_strat_pmoc_Pac.mean()
meanStratPacSO_UV14_anom_AreaWtd = wt_strat_anom_Pac.mean()



# Print out all the numbers
print('Mean GLOBAL Stratification (Grad-rho, 500m - surf), Full GLOBAL:')
print("UVIC '14")
print('  meanStratSO_UV14_ctrl_AreaWtd = '+str(meanStratSO_UV14_ctrl_AreaWtd.values))
print('  meanStratSO_UV14_pmoc_AreaWtd = '+str(meanStratSO_UV14_pmoc_AreaWtd.values))
print('  meanStratSO_UV14_anom_AreaWtd = '+str(meanStratSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanStratSO_UV14_anom_AreaWtd/meanStratSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)


print('\nMean GLOBAL Stratification (Grad-rho, 500m - surf), Pac-Sector GLOBAL:')
print("UVIC '14")
print('  meanStratPacSO_UV14_ctrl_AreaWtd = '+str(meanStratPacSO_UV14_ctrl_AreaWtd.values))
print('  meanStratPacSO_UV14_pmoc_AreaWtd = '+str(meanStratPacSO_UV14_pmoc_AreaWtd.values))
print('  meanStratPacSO_UV14_anom_AreaWtd = '+str(meanStratPacSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanStratPacSO_UV14_anom_AreaWtd/meanStratPacSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)



# %% --------- (3) WIND STRESS CURL ---------


# %% Define func to calc wind stress curl
# checked, good
def calc_WScurl(WSU_array, WSV_array, DPhi, model_name_str, make_plots=True):
    '''
    Parameters
    ----------
    WSU_array : 2D xr.DataArray of zonal wind stress, must have dims 'lat', 
                'lon'. Assumes lat/y dim is 1st (0) dimension.
    WSV_array : As above but for meridional wind stress, must be same as above
                 in dims, coords, and units
    DPhi : latitude spacing for the model, in [deg]. A single constant number.
    model_name_str : string describing model, used in plot titles.
    make_plots : logical to specify whether to make plots or not (Fig 1: 2x2,
                  WSU, WSV, dTaux/dy, dTauy/dx. Fig 2: WS curl w/ WS vectors.
                  Fig 3: 1x2, zonal mean WSU & zonal mean dTaux/dy)

    Returns
    -------
    curl_Tau : An wind stress curl xr.DataArray of same dims & coords as 
               WSU/WSV_array. Units = units of WSU/WSV /m
    fig1, fig2, fig3 : if make_plots=True, be sure to use this function as: 
                       curl, f1, f2, f3 = calc_WScurl(...)
    '''
    import numpy as np
    radius_Earth_m = 6371000 # NASA Earth Fact Sheet (https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)

    # For plotting
    lims_Taux = [-0.3, 0.3]
    lims_Tauy = [-0.2, 0.2]
    lims_dTxdy = [-3e-7, 3e-7]
    lims_dTydx = [-1.5e-7, 1.5e-7]
    lims_curlT = [-4e-7, 4e-7]
    my_cmap = 'RdBu_r'

    # Calc dy value (dPhi/2pi radians = dy meters/2piRadiusEarth)
    dy = np.deg2rad(DPhi)*radius_Earth_m  # just a constant scalar value
    # Will calc 'dx' as a func of lat w/in zonal grad loop

    # Meridional gradient (dTaux/dy) (https://github.com/fspaolo/altimpy/blob/master/build/lib/altimpy/util.py)    
    ydim=0
    dwsu_dy = np.gradient(WSU_array, dy)[ydim] # (N/m^3)
     
    # Zonal gradient (dTauy/dx)
    dwsv_dx = np.zeros(dwsu_dy.shape)
    num_lons = WSU_array.sizes['lon']
    num_lats = WSU_array.sizes['lat']
    lats = WSU_array.lat.values
    
    for lat_index in np.arange(np.size(lats)):
        dx_temp = (2*np.pi*radius_Earth_m*np.cos(np.deg2rad(lats[lat_index])))/num_lons
        # ^Circumf of Earth at that lat, divided by # longitude cells
        # Calc gradient (a 1D array), horz around the globe across a band of lat
        zon_grad = np.gradient(WSV_array.isel(lat=lat_index), dx_temp)
        # Build up the gradient array one row at a time
        dwsv_dx[lat_index,:] = zon_grad
        
    curl_Tau_data = dwsv_dx - dwsu_dy # (N/m^3 for UV, m/s /m for LC(?)

    # Build xr.DataArray
    curl_Tau = xr.DataArray(data=curl_Tau_data,
        dims=["lat", "lon"], coords=dict(
            lon=(["lon"], WSV_array.lon.values),
            lat=(["lat"], WSV_array.lat.values), ),
        attrs=dict(
            long_name="Wind Stress Curl",
            ), ) # units=WSV_array.units+' /m',
    
    if make_plots == True:
        # Fig 1
        fig1, axs = plt.subplots(2,2,figsize=(9,5))
        # Row 1: wind stress U and V
        cur_ax = axs[0,0]; cur_ax.set_facecolor('grey')
        WSU_array.plot(ax=cur_ax, vmin=lims_Taux[0], vmax=lims_Taux[1], cmap=my_cmap)
        cur_ax.set_title(WSU_array.long_name); cur_ax.set_xticklabels(''); cur_ax.set_xlabel('')
        
        cur_ax = axs[0,1]; cur_ax.set_facecolor('grey')
        WSV_array.plot(ax=cur_ax, vmin=lims_Tauy[0], vmax=lims_Tauy[1], cmap=my_cmap)
        cur_ax.set_title(WSV_array.long_name); cur_ax.set_ylabel(''); cur_ax.set_xticklabels(''); cur_ax.set_xlabel('')
        
        # Row 2: dTaux/dy, dTauy/dx
        cur_ax = axs[1,0]; cur_ax.set_facecolor('grey')
        CS = cur_ax.contourf(dwsu_dy, levels=np.arange(lims_dTxdy[0],lims_dTxdy[1],(lims_dTxdy[1]*2/100)), cmap=my_cmap)
        cbar = fig1.colorbar(CS); cur_ax.set_title('dTau_x/dy')
        
        cur_ax = axs[1,1]; cur_ax.set_facecolor('grey')
        CS = cur_ax.contourf(dwsv_dx, levels=np.arange(lims_dTydx[0],lims_dTydx[1],(lims_dTxdy[1]*2/100)), cmap=my_cmap)
        cbar = fig1.colorbar(CS); cur_ax.set_title('dTau_y/dx')
   
        fig1.suptitle(model_name_str)
 
   
        # Fig 2
        # Must build DataSet U first Then V
        ws_DataSet = xr.merge([WSU_array, WSV_array])
        fig2, ax = plt.subplots()
        curl_Tau.plot(ax=ax, vmin=lims_curlT[0], vmax=lims_curlT[1], cmap=my_cmap)
        # In quiver plot, take every 5th point in lon and every 3rd point in lat
        my_quiv = ws_DataSet.isel(lon=np.arange(0,num_lons,5), lat=np.arange(0,num_lats,3)).plot.quiver(x='lon', y='lat', u=list(ws_DataSet.keys())[0], v=list(ws_DataSet.keys())[1], ax=ax)
        veclength = 0.2 # m/s (legend arrow)
        maxstr = '%3.1f ' % veclength + WSU_array.units
        plt.quiverkey(my_quiv,0.1,0.1,veclength,maxstr,labelpos='S', coordinates='axes').set_zorder(11)
        plt.title(model_name_str+': WS Curl with WS Vectors')
        
        # Fig 3
        fig3, axs = plt.subplots()
        line1, = WSU_array.mean(dim='lon').plot()
        line2, = axs.plot(WSU_array.lat.values, np.nanmean(dwsu_dy,1)*1e6)
        axs.axhline(y=0.0, color='k', linestyle='-')
        axs.legend([line1, line2], ['Zonal wind stress (Tau_x)', 'Merid. Derivative (dTau_x/dy) *1e6'])
        plt.title(model_name_str)

        
    if make_plots == True:
        return curl_Tau, fig1, fig2, fig3
    else:
        return curl_Tau
    


# %% Calc wind stress curls
# checked, good.
# Remember to run as: 'curlTau, fig1, fig2, fig3 = calc_WScurl(...)' if 
#  make_plots=True (and to calc anom as '_pmoc[0] - _ctrl[0]' b/c will be tuples,
#  4 assignments made to it: curlTau, fig1, fig2, and fig3)

# Tested this func on LC14, looks good
# Checked all plots, look good (all same patterns and similar magnitudes). Now
#  turn off 'make_plots'. 

#                 calc_WScurl(WSU_array,    WSV_array,   DPhi, model_name_str, make_plots=True)
curlTau_UV_ctrl = calc_WScurl(tauX_UV_ctrl, tauY_UV_ctrl, 1.8, 'UV14 Ctrl', make_plots=False)
curlTau_UV_pmoc = calc_WScurl(tauX_UV_pmoc, tauY_UV_pmoc, 1.8, 'UV14 PMOC', make_plots=False)
curlTau_UV_anom = curlTau_UV_pmoc - curlTau_UV_ctrl



# %% Fig: S.O. WS Curl

lower_lat = -80
upper_lat = -30

my_facecolor = [0.4, 0.4, 0.4]

my_vmin = -3e-7
my_vmax = -my_vmin
my_cmap = 'PuOr'
my_vmin_anom = -1e-7
my_vmax_anom = -my_vmin_anom
my_cmap_anom = 'seismic'

fig, axs = plt.subplots(3,1, figsize=(11,9))


# ----- TOP ROW: ctrl ----- #

# ax 0: UVic, Ctrl
cur_ax = axs[0]; cur_ax.set_facecolor(my_facecolor)
var2plot = curlTau_UV_ctrl.sel(lat=slice(lower_lat, upper_lat))
my_plt = var2plot.plot(ax=cur_ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, add_colorbar=False)
cur_ax.axhline(y=-50, color='grey')
cur_ax.set_xlabel('')
cur_ax.set_ylabel('Control', fontweight='bold')
cur_ax.set_title('UVic', fontweight = 'bold')


# ----- MID ROW: pmoc ----- #

# ax 1: UVic, PMOC
cur_ax = axs[1]; cur_ax.set_facecolor(my_facecolor)
var2plot = curlTau_UV_pmoc.sel(lat=slice(lower_lat, upper_lat))
my_plt = var2plot.plot(ax=cur_ax, vmin=my_vmin, vmax=my_vmax, cmap=my_cmap, add_colorbar=False)
cur_ax.axhline(y=-50, color='grey')
cur_ax.set_xlabel('')
cur_ax.set_ylabel('Experiment', fontweight='bold')
cur_ax.set_title('')


# ----- BTM ROW: Anom ----- #
# On Anomaly plots, add contour lines of where (0) WScurl is in ctrl (--) & exper (-)
# And plot grey horizontal line at -50 (below which take average in next section)
# ax 2: UVic, Anom
cur_ax = axs[2]; cur_ax.set_facecolor(my_facecolor)
var2plot = curlTau_UV_anom.sel(lat=slice(lower_lat, upper_lat))
my_plt_2 = var2plot.plot(ax=cur_ax, vmin=my_vmin_anom, vmax=my_vmax_anom, cmap=my_cmap_anom, add_colorbar=False)
curlTau_UV_ctrl.sel(lat=slice(-65,-30)).plot.contour(ax=cur_ax, levels=[0], \
     colors='k', linewidths=1, linestyles='--')
curlTau_UV_pmoc.sel(lat=slice(-65,-30)).plot.contour(ax=cur_ax, levels=[0], \
     colors='k', linewidths=1)
cur_ax.axhline(y=-50, color='grey')
cur_ax.set_xlabel('Longitude [$^\circ$E]')
cur_ax.set_ylabel('Anomaly', fontweight='bold')
cur_ax.set_title('')
cur_ax.text(10, -75, 'Black lines denote region of upwelling, ctrl (dashed) & exp (solid)', bbox=dict(fc="lightgray"), fontsize=9)
 

# ----- Add colorbar ----- #
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.815, 0.39, 0.03, 0.49])
fig.colorbar(my_plt, cax=cbar_ax).set_label('Wind Stress Curl [Pa/m]', fontweight='bold')

# ----- Add 2nd colorbar ----- #
# fig.subplots_adjust(right=0.8)
cbar_ax_2 = fig.add_axes([0.815, 0.125, 0.03, 0.22])
fig.colorbar(my_plt_2, cax=cbar_ax_2).set_label('Curl Tau anomaly [Pa/m]', fontweight='bold')

# ----- Add overall title ----- #
fig.suptitle('Wind Stress Curl', fontweight='bold', fontsize=16, y=0.94)




# %% Calc area-wtd mean w.s. curl south of 50S

# For example see 'Area weighted average': 
# https://xgcm.readthedocs.io/en/latest/xgcm-examples/05_autogenerate.html?highlight=area-weighted%20average#Area-weighted-average


SO_upper_lat = -50
# Pac_west_lim = 150; Pac_east_lim=280 # 280 already set above


# First apply weights

wt_curlTau_ctrl = curlTau_UV_ctrl.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forTau)
wt_curlTau_pmoc = curlTau_UV_pmoc.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forTau)
wt_curlTau_anom = curlTau_UV_anom.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forTau)


wt_curlTau_ctrl_Pac = curlTau_UV_ctrl.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forTau)
wt_curlTau_pmoc_Pac = curlTau_UV_pmoc.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forTau)
wt_curlTau_anom_Pac = curlTau_UV_anom.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forTau)


# Then take avg

meanCurlTauSO_UV14_ctrl_AreaWtd = wt_curlTau_ctrl.mean()
meanCurlTauSO_UV14_pmoc_AreaWtd = wt_curlTau_pmoc.mean()
meanCurlTauSO_UV14_anom_AreaWtd = wt_curlTau_anom.mean()

meanCurlTauPacSO_UV14_ctrl_AreaWtd = wt_curlTau_ctrl_Pac.mean()
meanCurlTauPacSO_UV14_pmoc_AreaWtd = wt_curlTau_pmoc_Pac.mean()
meanCurlTauPacSO_UV14_anom_AreaWtd = wt_curlTau_anom_Pac.mean()



# Print out results

# Full numbers (for copying into Excel):
print('\nArea-Wtd Avg Wind Stress Curl, All S.O. (<50degS)')
print("UVIC '14")
print('  ctrl       = '+str(meanCurlTauSO_UV14_ctrl_AreaWtd.values))
print('  pmoc       = '+str(meanCurlTauSO_UV14_pmoc_AreaWtd.values))
print('  anom field = '+str(meanCurlTauSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanCurlTauSO_UV14_anom_AreaWtd/meanCurlTauSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)

# Print out results - Pac

print('\nArea-Wtd Avg Wind Stress Curl, Pac-Sector S.O. (<50degS)')
print("UVIC '14")
print('  ctrl       = '+str(meanCurlTauPacSO_UV14_ctrl_AreaWtd.values))
print('  pmoc       = '+str(meanCurlTauPacSO_UV14_pmoc_AreaWtd.values))
print('  anom field = '+str(meanCurlTauPacSO_UV14_anom_AreaWtd.values))
anom_asPerc = (meanCurlTauPacSO_UV14_anom_AreaWtd/meanCurlTauPacSO_UV14_ctrl_AreaWtd).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)




# %% --------- (4) BIOLOGY -------------------
# %% Read in output: npp

pathname = '../data/Output_2014_UVic_U-fwf_U-fNA/';
raw_npp_ctrl = xr.open_dataset(pathname+'tavg.02001.01.01.nc', decode_times=False)['O_phytnpp']
raw_npp = xr.open_dataset(pathname+'tavg_tot.nc', decode_times=False)['O_phytnpp']

# ocean net primary production rate    # mol N m-3 s-1
ctrl_npp = raw_npp_ctrl.isel(time=0); exp_npp = raw_npp.isel(time=-1)
anom_npp = exp_npp - ctrl_npp



# Rename dims from 'depth' to 'npzdlev'
ctrl_npp = ctrl_npp.rename({'latitude':'lat', 'longitude':'lon'})
exp_npp  =  exp_npp.rename({'latitude':'lat', 'longitude':'lon'})
anom_npp = anom_npp.rename({'latitude':'lat', 'longitude':'lon'})


# Convert units --> yr-1
ctrl_npp = ctrl_npp*365*24*60*60
exp_npp  = exp_npp*365*24*60*60
anom_npp = anom_npp*365*24*60*60
ctrl_npp = ctrl_npp.assign_attrs(units="mol N m-3 yr-1")
exp_npp  = exp_npp.assign_attrs(units="mol N m-3 yr-1")
anom_npp = anom_npp.assign_attrs(units="mol N m-3 yr-1")

# Also calc percentage-anom at each grid cell
percanom_npp = (anom_npp/ctrl_npp)*100
percanom_npp = percanom_npp.assign_attrs(units="% change of control")



# %% Fig: npp anomaly

# (npzdlev=0) == upper ~35m

# Anomaly
fig, ax = plt.subplots(figsize=(12,4))
anom_npp.isel(npzdlev=0).plot(ax=ax, cmap='seismic')
ax.set_facecolor('darkgrey')
ax.set_ylim(-80, -30)
plt.suptitle('NPP Rate Anomaly')

# Anomaly - as percentage of ctrl
fig, ax = plt.subplots(figsize=(12,4))
npp_vmin = 0.05; npp_vmax = -npp_vmin
percanom_npp.isel(npzdlev=0).plot(ax=ax, vmin=-50, vmax=50, cmap='seismic')
ax.set_facecolor('darkgrey')
ax.set_ylim(-80, -30)
plt.suptitle('NPP Rate Anomaly (% of ctrl)')


# %% Calc area-wt'd mean NPP & anomaly south of 40S

# Note, npp on same tracer grid as rho, therefore use dA_forRho

SO_upper_lat = -50

# First apply weights

wt_npp_ctrl = ctrl_npp.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forRho)
wt_npp_pmoc =  exp_npp.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forRho)
wt_npp_anom = anom_npp.sel(lat=slice(-90,SO_upper_lat)).weighted(dA_forRho)


wt_npp_ctrl_Pac = ctrl_npp.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_npp_pmoc_Pac = exp_npp.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)
wt_npp_anom_Pac = anom_npp.sel(lat=slice(-90,upper_lat),lon=slice(Pac_west_lim,Pac_east_lim)).weighted(dA_forRho)


# Then take avg

mean_npp_wt_ctrl = wt_npp_ctrl.mean()
mean_npp_wt_exp  = wt_npp_pmoc.mean()
mean_npp_wt_anom = wt_npp_anom.mean()

mean_npp_wt_ctrl_Pac = wt_npp_ctrl_Pac.mean()
mean_npp_wt_exp_Pac  = wt_npp_pmoc_Pac.mean()
mean_npp_wt_anom_Pac = wt_npp_anom_Pac.mean()



# Print out results

print('\nArea-Wtd Avg NPP, All S.O. (<50degS)')
print("UVIC '14")
print('  ctrl       = '+str(mean_npp_wt_ctrl.values))
print('  pmoc       = '+str(mean_npp_wt_exp.values))
print('  anom field = '+str(mean_npp_wt_anom.values))
anom_asPerc = (mean_npp_wt_anom/mean_npp_wt_ctrl).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)

# Print out results - Pac

print('\nArea-Wtd Avg NPP, Pac-Sector S.O. (<50degS)')
print("UVIC '14")
print('  ctrl       = '+str(mean_npp_wt_ctrl_Pac.values))
print('  pmoc       = '+str(mean_npp_wt_exp_Pac.values))
print('  anom field = '+str(mean_npp_wt_anom_Pac.values))
anom_asPerc = (mean_npp_wt_anom_Pac/mean_npp_wt_ctrl_Pac).values*100
print('  (anom/ctrl)*100 = %.4f %%' % anom_asPerc)


