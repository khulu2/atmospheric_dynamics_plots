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
plt.rcParams['text.usetex'] = True
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#-------------------------------------------------------------------------------------------------------
# Colour map stuff
#-------------------------------------------------------------------------------------------------------
jet = matplotlib.colormaps['jet']
ccm = jet(range(256))
new_ccm = ccm[128:-1, :]
zonal_wind_cmap = ListedColormap(new_ccm)

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
    #-------------------------------------------------------------------------------------------------------
    # Set up figure titles
    #-------------------------------------------------------------------------------------------------------
    timesteps = ['t = -24 hours', 't = -12 hours', 't = 0 hours', 't = +12 hours', 't = +24 hours']
    #-------------------------------------------------------------------------------------------------------
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
                plt.contour(trimmed_lons, lats, PV, levels = np.arange(-2, -1), colors = 'black', linewidths = 2, linestyles = 'solid')
                plt.contour(trimmed_lons, lats, BSPV, levels = np.arange(-2, -1), colors = 'blue', linewidths = 2, linestyles = 'solid')
                if (name_part == "P2" or name_part == "LC1") and (len(np.arange(round(np.min(np.min(MF))), -20, 100)) > 0):
                    chandles1 = plt.contour(trimmed_lons, lats, MF, levels = np.arange(round(np.min(np.min(MF))), -20, 100), colors = 'blue', linewidths = 0.5, linestyles = 'dashed')
                else:
                    chandles1 = plt.contour(trimmed_lons, lats, MF, levels = np.arange(20, round(np.max(np.max(MF))) + 100, 100), colors = 'red', linewidths = 0.5, linestyles = 'solid')
                
                chandles2 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3,::3], YF[::3,::3], color = 'black', scale = 250, width = 0.0009, pivot = 'mid', headwidth = 12, headlength = 12)
                ax.quiverkey(chandles2, X = 0.13, Y = 0.95, U = 10, label = r"10 m s$^{-1}$", labelpos='E')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude')
                ax.set_ylabel('relative latitude')
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
                # ax.clabel(chandles4,chandles4.levels[::3], fontsize = 8)
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_zonal_wind/RWB_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 15)
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
    #-------------------------------------------------------------------------------------------------------
    # Set up figure titles
    #-------------------------------------------------------------------------------------------------------
    timesteps = ['t = -24 hours', 't = -12 hours', 't = 0 hours', 't = +12 hours', 't = +24 hours']
    #-------------------------------------------------------------------------------------------------------
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
                # chandles2 = plt.contour(lons,lats,BSPV,levels = np.arange(-2,-1,1),colors = 'blue', linewidths = 2, linestyles = 'solid')
                chandles5 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3,::3], YF[::3,::3], color = 'black', scale = 250, width = 0.0009, pivot = 'mid', headwidth = 12, headlength = 12)
                ax.quiverkey(chandles5, X = 0.13, Y = 0.95, U = 10, label = r'10 m s$^{-1}$', labelpos='E')
                plt.contour(trimmed_lons, lats, U, levels = np.arange(28, 43 + 6, 6), colors = 'black', linewidths = 1.5, linestyles = 'dashed')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude')
                ax.set_ylabel('relative latitude')
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
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_eke/RWB_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 15)
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
    #-------------------------------------------------------------------------------------------------------
    # Set up figure titles
    #-------------------------------------------------------------------------------------------------------
    timesteps = ['t = -24 hours', 't = -12 hours', 't = 0 hours', 't = +12 hours', 't = +24 hours']
    #-------------------------------------------------------------------------------------------------------
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
                # chandles2 = plt.contour(lons,lats,BSPV,levels = np.arange(-2,-1,1),colors = 'blue', linewidths = 2, linestyles = 'solid')
                chandles5 = plt.quiver(trimmed_lons[::3], lats[::3], XF[::3,::3], YF[::3,::3], color = 'black', scale = 11000, width = 0.0009, pivot = 'mid', headwidth = 12, headlength = 12)
                ax.quiverkey(chandles5, X = 0.13, Y = 0.95, U = 500, label = r'500 m$^{2}$ s$^{-2}$ day$^{-1}$', labelpos='E')
                plt.contour(trimmed_lons, lats, EKE, levels = eke_levels, colors = 'black', linewidths = 0.8, linestyles = 'solid')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude')
                ax.set_ylabel('relative latitude')
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
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_ageo_flux_div/RWB_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 15)
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
    #-------------------------------------------------------------------------------------------------------
    # Set up figure titles
    #-------------------------------------------------------------------------------------------------------
    timesteps = ['t = -24 hours', 't = -12 hours', 't = 0 hours', 't = +12 hours', 't = +24 hours']
    #-------------------------------------------------------------------------------------------------------
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
                print(f'Plotting {name_part} events on the {qlv}K surface with ageo flux div. Label: {comp_label}. Max: {np.max(np.max(EKE))}. Min: {np.min(np.min(EKE))}')
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
                plt.contour(trimmed_lons, lats, EKE, levels = eke_levels, colors = 'black', linewidths = 0.8, linestyles = 'solid')
                #-------------------------------------------------------------------------------------------
                ax.set_xlabel('relative longitude')
                ax.set_ylabel('relative latitude')
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
                plt.colorbar(chandles0, orientation = 'horizontal')
                figure_title = f"{title_labels} ({name_part} at {qlv}K)" 
                save_title = f"e:/dissertation_figures/evolution_ageo_flux_div/RWB_{char_iterator}.pdf"
                plt.title(figure_title, fontsize = 15)
                char_iterator += 1
                if save_plots:
                    plt.savefig(save_title, dpi = 1000, bbox_inches = 'tight', pad_inches = 0)
                    plt.close()
                else:
                    plt.show()
                    plt.close()
    return

