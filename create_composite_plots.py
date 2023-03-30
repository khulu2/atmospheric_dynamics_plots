# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:55:47 2023

@author: Donald Nkabinde
"""
import plotting_module as pm
import numpy as np
#-------------------------------------------------------------------------------------------------------
comp_dir = 'f:/composite_means'
qlv_name = list(np.arange(310, 360, 10))
rwb_type = ['LC1', 'P2', 'LC2', 'P1']
labels = ['m04', 'm02', 'm00', 'p02', 'p04']
lons = np.arange(-180, 182.5, 2.5)
lats = np.arange(27.5, -30, -2.5)
pressure_levels = [1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10]
pressure_levels = np.array(pressure_levels)
#-------------------------------------------------------------------------------------------------------
pm.plot_zonal_wind_composites(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_eke_composites(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_ageo_flux_div_composites(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_eke_tendency_composites(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_downstream_development(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_strain_rate_composites(comp_dir, rwb_type, lons, lats, -60, 60, qlv_name, labels, True)
pm.plot_meridionally_averaged_composites(comp_dir, rwb_type, lons, -60, 60, qlv_name, pressure_levels, labels, True)