# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:31:26 2023

@author: Donald Nkabinde
"""
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
# import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['text.usetex'] = True
import matplotlib
# matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'
#-------------------------------------------------------------------------------------------------------
# Set up figure titles
#-------------------------------------------------------------------------------------------------------
timesteps = ['t = -24 hours', 't = -12 hours', 't = 0 hours', 't = +12 hours', 't = +24 hours']
#-------------------------------------------------------------------------------------------------------
# Colour map stuff
#-------------------------------------------------------------------------------------------------------
jet = matplotlib.colormaps['jet']
ccm = jet(range(256))
new_ccm = ccm[128:-1, :]
zonal_wind_cmap = ListedColormap(new_ccm)

jet = matplotlib.colormaps['jet_r']
ccm = jet(range(256))
new_ccm = ccm[128:-1,:]
jet_blue = ListedColormap(new_ccm)

cmap = matplotlib.colormaps['bwr']
ccm = cmap(range(256))
white = np.array([1, 1, 1, 1])
ccm[120:137, :] = white
new_bwr = ListedColormap(ccm)

cmap = matplotlib.colormaps['coolwarm']
ccm = cmap(range(256))
white = np.array([1, 1, 1, 1])
ccm[120:137, :] = white
new_coolwarm = ListedColormap(ccm)
#-------------------------------------------------------------------------------------------------------
def plot_zonal_wind_composites(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot zonal wind evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the basic state PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/bs_pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                BSPV = fh.variables[f"bs_pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                U = fh.variables[f"u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the perturbation zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eddy_u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                XF = fh.variables[f"eddy_u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the perturbation meridional wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eddy_v_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                YF = fh.variables[f"eddy_v_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Momentum flux
                #-------------------------------------------------------------------------------------------
                MF = XF * YF
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with zonal wind. Label: {comp_label}. Max: {np.max(np.max(U))}. Min: {np.min(np.min(U))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                
                levels = np.arange(28, 43) 
                
                chandles0 = plt.contourf(trimmed_lons, lats, U, levels = levels, cmap = zonal_wind_cmap)
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2, -1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.contour(trimmed_lons, lats, BSPV, levels = np.arange(-2, -1), colors = 'blue', linewidths = 2, linestyles = 'solid')
                if (name_part == "P2" or name_part == "LC1") and (len(np.arange(round(np.min(np.min(MF))), -20, 100)) > 0):
                    chandles3 = plt.contour(trimmed_lons, lats, MF, levels = np.arange(round(np.min(np.min(MF))), -20, 100), colors = 'blue', linewidths = 0.5, linestyles = 'dashed')
                else:
                    chandles3 = plt.contour(trimmed_lons, lats, MF, levels = np.arange(20, round(np.max(np.max(MF))) + 100, 100), colors = 'red', linewidths = 0.5, linestyles = 'solid')
                
                chandles4 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3, ::3], YF[::3, ::3], color = 'black', scale = 250, width = 0.0028, pivot = 'mid', headwidth = 6, headlength = 6)
                ax.quiverkey(chandles4, X = 0.13, Y = 0.95, U = 10, label = r"10 m s$^{-1}$", labelpos='E')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':',linewidth = 0.5)
                ax.clabel(chandles3, chandles3.levels[::2], fontsize = 8)
                ax.clabel(chandles1, chandles1.levels, fontsize = 8)
                ax.clabel(chandles2, chandles2.levels, fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_zonal_wind/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_eke_composites(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot eddy kinetic energy evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the EKE file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eke_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                EKE = fh.variables[f"eke_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                U = fh.variables[f"u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the perturbation zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eddy_u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                XF = fh.variables[f"eddy_u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the perturbation meridional wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eddy_v_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                YF = fh.variables[f"eddy_v_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with EKE. Label: {comp_label}. Max: {np.max(np.max(U))}. Min: {np.min(np.min(U))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                
                if name_part == "LC1":
                    if q == 310:
                        levels = np.arange(120, 275 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(120, 275 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(120, 200 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(120, 155 + 5, 5) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(120, 145 + 5, 5) # Levels we want to show on the colorbar
                elif name_part == "P2":                    
                    if q == 310:
                        levels = np.arange(120, 235 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(120, 250 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(120, 245 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(120, 240 + 5, 5) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(120, 210 + 5, 5) # Levels we want to show on the colorbar
                elif name_part == "P1":
                    if q == 310:
                        levels = np.arange(120, 190 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(120, 270 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(120, 260 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(120, 190 + 5, 5) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(120, 165 + 5, 5) # Levels we want to show on the colorbar
                else:
                    if q == 310:
                        levels = np.arange(120, 235 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(120, 245 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(120, 180 + 5, 5) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(120, 190 + 5, 5) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(120, 200 + 5, 5) # Levels we want to show on the colorbar
                            
                chandles0 = plt.contourf(trimmed_lons, lats, EKE, levels = levels, cmap = 'cubehelix_r')
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2, -1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3,::3], YF[::3,::3], color = 'black', scale = 250, width = 0.0028, pivot = 'mid', headwidth = 6, headlength = 6)
                ax.quiverkey(chandles2, X = 0.13, Y = 0.95, U = 10, label = r'10 m s$^{-1}$', labelpos='E')
                chandles3 = plt.contour(trimmed_lons, lats, U, levels = np.arange(28, 43 + 6, 6), colors = 'black', linewidths = 1.5, linestyles = 'dashed')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':', linewidth = 0.5)
                ax.clabel(chandles1, chandles1.levels[::2], fontsize = 8)
                ax.clabel(chandles3, chandles3.levels, fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_eke/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_ageo_flux_div_composites(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot ageostrophic flux divergence evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the EKE file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eke_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                EKE = fh.variables[f"eke_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the ageostrophic flux divergence file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/ageo_horizontal_flux_div_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                F = fh.variables[f"ageo_horizontal_flux_div_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                F = -F * 1e4
                #-------------------------------------------------------------------------------------------
                # Open the zonal flux file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/ageo_flux_x_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                XF = fh.variables[f"ageo_flux_x_{qlv}"][:, min_lon_idx:max_lon_idx + 1]/4
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the meridional flux file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/ageo_flux_y_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                YF = fh.variables[f"ageo_flux_y_{qlv}"][:, min_lon_idx:max_lon_idx + 1]/4
                fh.close()
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with ageo flux div. Label: {comp_label}. Max: {np.max(np.max(EKE))}. Min: {np.min(np.min(EKE))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                
                if name_part == "LC1":
                    if q == 310:
                        levels = np.arange(-50, 50 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 320:
                        levels = np.arange(-60, 60 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 330:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 340:
                        levels = np.arange(-18, 18 + 3, 3) # Levels we want to show on the colorbar
                        eke_levels = np.arange(100, np.max(np.max(EKE)) + 30, 30)
                    else:
                        levels = np.arange(-12, 12 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(100, np.max(np.max(EKE)) + 30, 30)
                elif name_part == "P2":                    
                    if q == 310:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 320:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 330:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 340:
                        levels = np.arange(-35 , 35 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    else:
                        levels = np.arange(-50, 50 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                elif name_part == "P1":
                    if q == 310:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(125, np.max(np.max(EKE)) + 25, 25) # EKE levels we want to show
                    elif q == 320:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 330:
                        levels = np.arange(-50, 50 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(180, np.max(np.max(EKE)) + 25, 25)
                    elif q == 340:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    else:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(130, np.max(np.max(EKE)) + 25, 25)
                else:
                    if q == 310:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 320:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 330:
                        levels = np.arange(-28, 28 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30, 30)
                    elif q == 340:
                        levels = np.arange(-10, 10 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30 ,30)
                    else:
                        levels = np.arange(-15, 15 + 3, 3) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30, 30)
                            
                chandles0 = plt.contourf(trimmed_lons, lats, F, levels = levels, cmap = new_bwr)
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2,-1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3, ::3], YF[::3, ::3], color = 'black', scale = 11000, width = 0.0028, pivot = 'mid', headwidth = 6, headlength = 6)
                ax.quiverkey(chandles2, X = 0.13, Y = 0.95, U = 500, label = r'500 m$^{2}$ s$^{-2}$ day$^{-1}$', labelpos='E')
                chandles3 = plt.contour(trimmed_lons, lats, EKE, levels = eke_levels, colors = 'black', linewidths = 0.8, linestyles = 'solid')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':', linewidth = 0.5)
                ax.clabel(chandles1, chandles1.levels[::2], fontsize = 8)
                ax.clabel(chandles3, chandles3.levels[::3], fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_ageo_flux_div/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_eke_tendency_composites(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot ageostrophic flux divergence evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the EKE file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eke_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                EKE = fh.variables[f"eke_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the EKE tendency file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/eke_tendency_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                F = fh.variables[f"eke_tendency_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with EKE tendency. Label: {comp_label}. Max: {np.max(np.max(EKE))}. Min: {np.min(np.min(EKE))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                
                if name_part == "LC1":
                    if q == 310:
                        levels = np.arange(-75, 75 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 320:
                        levels = np.arange(-75, 75 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 330:
                        levels = np.arange(-50, 50 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                    elif q == 340:
                        levels = np.arange(-20, 20 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(100, np.max(np.max(EKE)) + 30, 30)
                    else:
                        levels = np.arange(-14, 14 + 1) # Levels we want to show on the colorbar
                        eke_levels = np.arange(100, np.max(np.max(EKE)) + 30, 30)
                elif name_part == "P2":                    
                    if q == 310:
                        levels = np.arange(-45, 45 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 320:
                        levels = np.arange(-60, 60 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 330:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    elif q == 340:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30) # EKE levels we want to show
                    else:
                        levels = np.arange(-36, 36 + 3, 3) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 30, 30)
                elif name_part == "P1":
                    if q == 310:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(125, np.max(np.max(EKE)) + 25, 25) # EKE levels we want to show
                    elif q == 320:
                        levels = np.arange(-60, 60 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 330:
                        levels = np.arange(-60, 60 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(180, np.max(np.max(EKE)) + 25, 25)
                    elif q == 340:
                        levels = np.arange(-26, 26 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    else:
                        levels = np.arange(-16, 16 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(130, np.max(np.max(EKE)) + 25, 25)
                else:
                    if q == 310:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 320:
                        levels = np.arange(-40, 40 + 5, 5) # Levels we want to show on the colorbar
                        eke_levels = np.arange(150, np.max(np.max(EKE)) + 25, 25)
                    elif q == 330:
                        levels = np.arange(-30, 30 + 2, 2) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30, 30)
                    elif q == 340:
                        levels = np.arange(-36, 36 + 3, 3) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30 ,30)
                    else:
                        levels = np.arange(-36, 36 + 3, 3) # Levels we want to show on the colorbar
                        eke_levels = np.arange(120, np.max(np.max(EKE)) + 30, 30)
                            
                chandles0 = plt.contourf(trimmed_lons, lats, F, levels = levels, cmap = new_bwr)
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2,-1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.contour(trimmed_lons, lats, EKE, levels = eke_levels, colors = 'black', linewidths = 0.8, linestyles = 'solid')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':', linewidth = 0.5)
                ax.clabel(chandles1, chandles1.levels[::2], fontsize = 8)
                ax.clabel(chandles2, chandles2.levels[::3], fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_eke_tendency/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_downstream_development(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot phi prime and agestrophic wind evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the perturbation geopotential file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/z_{qlv}_{comp_label}_3D.nc"
                fh = nc.Dataset(fn,'r')
                Z = fh.variables[f"z_{qlv}"][9, :, min_lon_idx:max_lon_idx + 1] # phi prime at 250hPa
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                U = fh.variables[f"u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the zonal ageostrophic wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/ua_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                XF = fh.variables[f"ua_{qlv}"][:, min_lon_idx:max_lon_idx + 1] * 10
                fh.close()
                XF[XF > 40] = np.nan
                XF[XF < -40] = np.nan
                #-------------------------------------------------------------------------------------------
                # Open the meridional ageostrophic wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/va_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                YF = fh.variables[f"va_{qlv}"][:, min_lon_idx:max_lon_idx + 1] * 10
                fh.close()
                YF[YF > 40] = np.nan
                YF[YF < -40] = np.nan
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with downstream development. Label: {comp_label}. Max: {np.max(np.max(U))}. Min: {np.min(np.min(U))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                
                if name_part == "LC1":
                    if q == 310:
                        levels = np.arange(-1600, 1600 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(-1600, 1600 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(-1300, 1300 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(-1000, 1000 + 100, 100) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(-800, 800 + 100, 100) # Levels we want to show on the colorbar
                elif name_part == "P2":                    
                    if q == 310:
                        levels = np.arange(-1700, 1700 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(-2100, 2100 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(-2100, 2100 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(-1800, 1800 + 100, 100) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(-1500, 1500 + 100, 100) # Levels we want to show on the colorbar
                elif name_part == "P1":
                    if q == 310:
                        levels = np.arange(-1200, 1200 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(-1900, 1900 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(-2100, 2100 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(-1600, 1600 + 100, 100) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(-1200, 1200 + 100, 100) # Levels we want to show on the colorbar
                else:
                    if q == 310:
                        levels = np.arange(-1800, 1800 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 320:
                        levels = np.arange(-1500, 1500 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 330:
                        levels = np.arange(-1100, 1100 + 100, 100) # Levels we want to show on the colorbar
                    elif q == 340:
                        levels = np.arange(-800, 800 + 100, 100) # Levels we want to show on the colorbar
                    else:
                        levels = np.arange(-600, 600 + 50, 50) # Levels we want to show on the colorbar
                            
                chandles0 = plt.contourf(trimmed_lons, lats, Z, levels = levels, cmap = new_coolwarm)
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2,-1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles5 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3, ::3]/2, YF[::3, ::3]/2, color = 'black', scale = 250, width = 0.0028, pivot = 'mid', headwidth = 6, headlength = 6)
                ax.quiverkey(chandles5, X = 0.13, Y = 0.95, U = 10, label = r'10 m s$^{-1}$', labelpos='E')
                chandles3 = plt.contour(trimmed_lons, lats, U, levels = np.arange(28, 43 + 6, 6), colors = 'black', linewidths = 1.5, linestyles = 'dashed')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':', linewidth = 0.5)
                ax.clabel(chandles1, chandles1.levels[::2], fontsize = 8)
                ax.clabel(chandles3, chandles3.levels, fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_phi_prime/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_strain_rate_composites(comp_dir, rwb_types, lons, lats, min_lon, max_lon, qlv_name, labels, save_plots = False):
    """
    Plot rate of strain evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                PV = fh.variables[f"pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the basic state PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/bs_pv_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                BSPV = fh.variables[f"bs_pv_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the zonal wind file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/u_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                U = fh.variables[f"u_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                #-------------------------------------------------------------------------------------------
                # Open the magnitude of strain file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/strain_{qlv}_{comp_label}.nc"
                fh = nc.Dataset(fn,'r')
                SR = fh.variables[f"strain_{qlv}"][:, min_lon_idx:max_lon_idx + 1]
                fh.close()
                SR = SR * 1e6
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} events on the {qlv}K surface with strain rate. Label: {comp_label}. Max: {np.max(np.max(U))}. Min: {np.min(np.min(U))}')
                #-------------------------------------------------------------------------------------------
                ax = plt.axes(projection = ccrs.PlateCarree())
                                
                levels = np.arange(3, 8.5, .5) 
            
                chandles0 = plt.contourf(trimmed_lons, lats, SR, levels = levels, cmap = jet_blue)
                chandles1 = plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2,-1), colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.contour(trimmed_lons, lats, BSPV, levels = np.arange(-2,-1), colors = 'blue', linewidths = 2, linestyles = 'solid')
                chandles3 = plt.contour(trimmed_lons, lats, U, levels = np.arange(28, 50, 5), colors = 'black', linewidths = 2, linestyles = 'dashed')               
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude', fontsize = 12)
                ax.set_ylabel('relative latitude', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = np.arange(min_lon, max_lon + 30, 30)
                ax.set_xticks(x_tick_labels)
                y_tick_labels = np.arange(-20, 30, 10)
                ax.set_yticks(y_tick_labels)
                ax.set_ylim([min(lats), max(lats)])
                ax.set_xlim([min(trimmed_lons), max(trimmed_lons)])
                lon_formatter = cticker.LongitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                lat_formatter = cticker.LatitudeFormatter()
                ax.yaxis.set_major_formatter(lat_formatter)
                plt.grid(linestyle = ':',linewidth = 0.5)
                ax.clabel(chandles1, chandles1.levels[::2], fontsize = 8)
                ax.clabel(chandles2, chandles2.levels[::2], fontsize = 8)
                ax.clabel(chandles3, chandles3.levels, fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_strain_rate/RWB_{qlv}_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

def plot_meridionally_averaged_composites(comp_dir, rwb_types, lons, min_lon, max_lon, qlv_name, pressure_levels, labels, save_plots = False):
    """
    Plot meridionally averaged variables evolution composite means.

    Args:
        comp_dir: 
            The directory where the netCDF files are stored containing the data for the composites.
        rwb_types: 
            A list of strings representing the different types of RWB that will be plotted.
        lons and lats: 
            Arrays containing the longitudes and latitudes of the data grid.
        min_lon and max_lon: 
            The minimum and maximum longitudes to be included in the plot.
        qlv_name: 
            A list of integers representing the isentropic levels to be plotted.
        labels: 
            A list of strings representing the labels for each composite plot.
        save_plots: 
            A Boolean variable that determines whether or not the plot should be saved to a file.
    """
    min_lon_idx = np.where(lons == min_lon)[0][0]
    max_lon_idx = np.where(lons == max_lon)[0][0]
    trimmed_lons = lons[min_lon_idx:max_lon_idx + 1]
    for q in qlv_name:
        qlv = str(q)
        char_iterator = 0
        for name_part in rwb_types:
            for ic in np.arange(len(labels)):
                letter = chr(97 + char_iterator)
                comp_label = labels[ic]
                comp_time = timesteps[ic]
                title_labels = f"({letter}) {comp_time}"
                #-------------------------------------------------------------------------------------------
                # Open the PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/meridional_average_pv_{qlv}_{comp_label}.txt"
                PV = np.loadtxt(fn, delimiter=',')
                PV = PV * 1e6
                PV = PV[min_lon_idx:max_lon_idx + 1, :]
                #-------------------------------------------------------------------------------------------
                # Open the basic state PV file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/meridional_average_bs_pv_{qlv}_{comp_label}.txt"
                BPV = np.loadtxt(fn, delimiter=',')
                BPV = BPV * 1e6
                BPV = BPV[min_lon_idx:max_lon_idx + 1, :]
                #-------------------------------------------------------------------------------------------
                # Open the eddy geopotential file
                #-------------------------------------------------------------------------------------------
                fn = f"{comp_dir}/{name_part}/meridional_average_z_{qlv}_{comp_label}.txt"
                Z = np.loadtxt(fn, delimiter=',')
                Z = Z[min_lon_idx:max_lon_idx + 1, :]
                #-------------------------------------------------------------------------------------------
                # PV anomaly
                #-------------------------------------------------------------------------------------------
                APV = PV - BPV
                #-------------------------------------------------------------------------------------------
                print(f'Plotting {name_part} cross sections of events on the {qlv}K surface. Label: {comp_label}. Max: {np.max(np.max(APV))}. Min: {np.min(np.min(APV))}')
                #-------------------------------------------------------------------------------------------
                fig, ax = plt.subplots(figsize=(8, 8))
            
                # levels = np.arange(-84,84+4,4)
                levels = np.arange(-1, 1.1, 0.1)
                
                chandles0 = plt.contourf(trimmed_lons, np.arange(1,18), np.transpose(APV), levels = levels, cmap = new_bwr)
                chandles1 = plt.contour(trimmed_lons, np.arange(1,18), np.transpose(PV), levels = np.arange(-2,-1),colors = 'black', linewidths = 2, linestyles = 'solid')
                chandles2 = plt.contour(trimmed_lons, np.arange(1,18), np.transpose(Z), levels = np.arange(-1000,-50,50), colors = 'red', linewidths = 1, linestyles = 'dashed')
                chandles3 = plt.contour(trimmed_lons, np.arange(1,18), np.transpose(Z), levels = np.arange(100,1050,50), colors = 'blue', linewidths = 1, linestyles = 'solid')               
                #-------------------------------------------------------------------------------------------
                plt.yscale('log')
                ax.set_xlabel('latitude', fontsize = 12)
                ax.set_ylabel('pressure level (hPa)', fontsize = 12)
                ax.tick_params(direction = 'in')    # Tick direction... Similar to m_grid('tickdir','in') in MATLAB
                x_tick_labels = trimmed_lons[::4]
                ax.set_xticks(x_tick_labels)
                ax.set_xticklabels(['60°W','50°W','40°W','30°W','20°W','10°W','0°','10°E','20°E','30°E','40°E','50°E','60°E'])
                y_tick_labels = np.arange(1,18)
                ax.set_yticks(y_tick_labels)
                ax.set_yticklabels(["".join(item) for item in pressure_levels.astype(str)])
                plt.grid(linestyle = ':',linewidth = 0.5)
                ax.clabel(chandles1,chandles1.levels, fontsize = 6)
                ax.clabel(chandles2,chandles2.levels[::3], fontsize = 6)
                ax.clabel(chandles3,chandles3.levels[::3], fontsize = 6)
                
                divider = make_axes_locatable(ax)
                cax = divider.new_vertical(size = "5%",
                                    pad = 0.4,
                                    pack_start = True)
                fig.add_axes(cax)
                cbar = plt.colorbar(chandles0,cax = cax, orientation = 'horizontal', ticks = levels[::5])
                cbar.ax.set_xticklabels(["".join(item) for item in np.arange(-1,1.5,0.5).astype(str)])
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_meridional_averages_composites/RWB_{qlv}_{char_iterator}.pdf"
                ax.set_title(figure_title, fontsize = 20)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return
