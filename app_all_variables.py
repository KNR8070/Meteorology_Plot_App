#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 20:54:39 2024

@author: knreddy
"""
#%% [markdown]
## load modules

import streamlit as st
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature
import calendar
import cmcrameri.cm as cmc
from streamlit_option_menu import option_menu
import geopandas as gpd

st.set_page_config(layout="wide")
#%% [markdown] 
## Load wind data from NetCDF
@st.cache_data
def load_shapefile():
    India = gpd.read_file('data/India_Shapefile_/india_administrative_outline_boundary.shp', crs="epsg:4326")
    # If CRS is missing, set it manually
    if India.crs is None:
        India.set_crs("EPSG:4326", inplace=True)
    
    return India
@st.cache_data
def load_temp_data():
    ds_temp = xr.open_dataset("data/air.2m.mon.ltm.1991-2020.nc")
    return ds_temp
@st.cache_data
def load_uwind_data():
    ds_u = xr.open_dataset("data/uwnd.mon.ltm.1991-2020.nc")
    return ds_u
@st.cache_data
def load_vwind_data():
    ds_v = xr.open_dataset("data/vwnd.mon.ltm.1991-2020.nc")
    return ds_v
@st.cache_data
def load_pr_data():
    ds_pr = xr.open_dataset("data/precip.mon.ltm.1991-2020.nc")
    return ds_pr
@st.cache_data
def load_rh_data():   
    ds_rh = xr.open_dataset("data/rhum.mon.ltm.1991-2020.nc")
    return ds_rh

# %% Anomaly plotting
def anomaly_plotting(clim_var,mon, region):
    [lat_min, lat_max, lon_min, lon_max] = user_input_region(region)
    # Load the data
    var_data = xr.open_dataset('data/era5_2024_2025_'+ clim_var.name+'.nc')
    var_data = convert_180_180(var_data).sel(lat=slice(85,-85),lon=slice(-176,176))
    var_subset = var_data[clim_var.name].sel(lat=slice(lat_max, lat_min), 
                                        lon=slice(lon_min, lon_max),
                                        level=1000)
    clim_var_subset = clim_var.sel(lat=slice(lat_max, lat_min),
                                    lon=slice(lon_min, lon_max))
    
    if len(var_subset.time[var_subset.time.dt.month == mon]) == 2:
        year = np.unique(var_subset.time.dt.year.values)[1]
        sub_data = var_subset[var_subset.time.dt.year == year]
        var_plot_data = np.squeeze(sub_data.isel(time=mon-1))-273.15
    else:
        year = np.unique(var_subset.time.dt.year.values)[0]
        sub_data = var_subset[var_subset.time.dt.year == year]
        var_plot_data = np.squeeze(sub_data.isel(time=mon-1))-273.15
###### Plotting
    x_size, y_size = calculate_x_y_size(lat_min, lat_max, lon_min, lon_max)
    fig, ax3 = plt.subplots(figsize=(x_size,y_size),
                            subplot_kw={'projection': ccrs.PlateCarree()})#,
    ax3.set_extent([var_subset.lon.values.min(),#lon_min,
                    var_subset.lon.values.max(),#lon_max, 
                    var_subset.lat.values.min(),#lat_min,
                    var_subset.lat.values.max()],#lat_max,], 
                    crs=ccrs.PlateCarree()) 
    India = load_shapefile()
    shape_feature = ShapelyFeature(India.geometry,
                                ccrs.PlateCarree(), edgecolor='black', 
                                facecolor='none')
    if (region != 'Global') and (region != 'India') and (region != 'China'):  
        ax3.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax3.add_feature(cfeature.COASTLINE)
    elif (region == 'India'):
        ax3.add_feature(shape_feature, linewidth=1.0, facecolor='none')
        #India.plot(facecolor='none',edgecolor='black',ax=ax3)
        ax3.add_feature(cfeature.COASTLINE)
    elif (region == 'China'):
        ax3.add_feature(shape_feature, linewidth=1.0, facecolor='none')
        ax3.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax3.add_feature(cfeature.COASTLINE)
    else:
        ax3.add_feature(cfeature.COASTLINE)
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(clim_var.lon, clim_var.lat)
    level_in_feet = {"1000.0": 'Surface',
                     "925.0": '3000 ft',
                     "850.0": '5000 ft',
                     "700.0": '10000 ft',
                     "500.0": '18000 ft'}

    if clim_var.var_desc=='Air temperature':
        clim_plot_data = np.squeeze(clim_var_subset.isel(time=mon-1))-273.15
        plot_data = var_plot_data - clim_plot_data
        s_plot = ax3.contourf(lons,lats,plot_data,
                              cmap=cmc.vik, 
                              vmin=-20, vamx=20,
                              levels=np.linspace(-20,20,41),
                              extend='both')
        if x_size<y_size:
            cbar = fig.colorbar(s_plot, ax=ax3,shrink=0.3)# label="2m Temperature (degC)",                  
        else:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)# label="2m Temperature (degC)", 
        cbar.set_label('2m Temperature (degC)',size='xx-small')
        ax3.set_title('Anomaly in '+calendar.month_name[mon]+'  '+ str(year)+ 
                      ' compared to 1991-2021 clim.', size='x-small')            
        
    else:
        plot_data = np.squeeze(var_subset.isel(time=mon-1))
        s_plot = ax3.contourf(lons,lats,plot_data,
                              cmap=cmc.batlowW, 
                              vmin=0,vmax=100,
                              levels=np.linspace(0,100,26))
        if x_size<y_size:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.3)# label="Relative humidity (%)",                     
        else:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)# label="Relative humidity (%)",                  
        cbar.set_label('Relative humidity (%)',size='xx-small')
        if str(var_subset.level.values) in level_in_feet:
            ax3.set_title('Month:'+calendar.month_name[mon]+
                         '  Level:'+str(var_subset.level.values)+
                         ' '+var_subset.level.GRIB_name+' ('+level_in_feet[str(var_subset.level.values)]+')', 
                         size='medium')
        else:
            ax3.set_title('Month:'+calendar.month_name[mon]+
                         '  Level:'+str(var_subset.level.values)+
                         ' '+var_subset.level.GRIB_name,
                         size='medium')  
    cbar.ax.tick_params(labelsize='xx-small')
    
    #ax3.set_title(calendar.month_name[mon][:3],size='small')
    #ax3.set_title('Month:'+calendar.month_name[mon][:3]+
    #             '  Level:'+str(var_subset.level.values)+
    #             ' '+var_subset.level.GRIB_name, size='x-small')
    
    

    ax3.set_xlabel('Longitude',size='x-small')
    ax3.set_ylabel('Latitude',size='x-small')

    ax3.set_xticks(np.linspace(np.floor(clim_var.lon.values.min()),#lon_min,
                               np.floor(clim_var.lon.values.max()),#lon_max,
                               num=5,endpoint=True))#format='%2.2f'))
    ax3.set_yticks(np.linspace(np.floor(clim_var.lat.values.min()),#lat_min,
                               np.floor(clim_var.lat.values.max()),#lat_max,
                               num=5,endpoint=True))#,format='%2.2f'))
    ax3.set_xticklabels(np.linspace(np.floor(clim_var.lon.values.min()),#lon_min,
                                    np.floor(clim_var.lon.values.max()),#lon_max,
                                    num=5,endpoint=True),#format='%2.2f'),
                                    size='xx-small')
    ax3.set_yticklabels(np.linspace(np.floor(clim_var.lat.values.min()),#lat_min,
                                    np.floor(clim_var.lat.values.max()),#lat_max,
                                    num=5,endpoint=True),#format='%2.2f'),
                                    size='xx-small')
    st.pyplot(fig)
#%% [markdown]
# Calculate wind speed and direction from u10 and v10
def calculate_wind(uwnd,vwnd):
    speed = np.sqrt(uwnd**2 + vwnd**2)
    direction = (np.arctan2(vwnd, uwnd) * 180 / np.pi + 180) % 360
    return speed, direction
#%% [markdown] 
# Function to select the spatial plot box
def select_box(lat,lon):
    # Latitude and Longitude Box Selection
    if ((lat_min > lat_max) or lon_min > lon_max):
        st.text("ERROR: Minimum value greater than Maximum")    
    return lat_min, lat_max, lon_min, lon_max
#%% [markdown] 
# Function to plot windrose Plot
def plot_wind_rose(speed_pwr, direction_pwr,lat_l,lon_l):
    fig_ws = plt.figure(figsize=(4, 4))
    ax = WindroseAxes.from_ax()
    ax.bar(direction_pwr, speed_pwr, normed=True, opening=0.8, edgecolor='white')
    ax.set_title('Latitude = '+str(lat_l)+' and Longitude = '+str(lon_l))
    ax.set_legend(title="Wind Speed (m/s)",loc='best')
    ax.text(0.7,-0.0,'Data Source: '+ds_temp.attrs['source'],
            fontsize=6,transform=ax.transAxes)
    st.pyplot()  
#%% [markdown]
# calculating figsize for spatial plots
def calculate_x_y_size(lat_min, lat_max, lon_min, lon_max):
    ### Dynamic plot size
    lat_range = lat_max-lat_min
    lon_range = lon_max-lon_min
    
    x_ratio = lon_range/(lon_range+lat_range)
    y_ratio = lat_range/(lon_range+lat_range)
    
    if x_ratio>y_ratio:
        ratio_xy = np.round(x_ratio/y_ratio)
        x_size = np.round(5)
        y_size = np.round(7)
    elif x_ratio==y_ratio:
        x_size = np.round(7)
        y_size = np.round(5)
    else:
        ratio_yx = np.round(y_ratio/x_ratio)
        x_size = np.round(6)
        y_size = np.round(5)
    return x_size, y_size
#%% [markdown] 
# Function to plot wind vector plot
def plot_wind_vectors(ds_u,ds_v, lat_min, lat_max, lon_min, lon_max, time_s, region):
    # Select data within specified lat/lon box   
    speed_mean = np.sqrt(ds_u**2 + ds_v**2)
    x_size, y_size = calculate_x_y_size(lat_min, lat_max, lon_min, lon_max)
    fig, ax = plt.subplots(figsize=(x_size,y_size),
                           subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([speed_mean.lon.values.min(),#lon_min, 
                   speed_mean.lon.values.max(),#lon_max, 
                   speed_mean.lat.values.min(),#lat_min, 
                   speed_mean.lat.values.max(),#lat_max
                   ], 
                   crs=ccrs.PlateCarree())
    India = load_shapefile()
    shape_feature = ShapelyFeature(India.geometry,
                                ccrs.PlateCarree(), edgecolor='black', 
                                facecolor='none')
    if (region != 'Global') and (region != 'India'):  
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE)
    elif region == 'India':
        #India.plot(facecolor='none',edgecolor='black',ax=ax)
        ax.add_feature(shape_feature, linewidth=1.0)
        ax.add_feature(cfeature.COASTLINE)
    else:
        ax.add_feature(cfeature.COASTLINE)
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(ds_u.lon, ds_u.lat)
    speed_plot = ax.contourf(lons, lats, speed_mean, cmap=cmc.batlowW_r, extend='both')
    if x_size<y_size:
        cbar = fig.colorbar(speed_plot, ax=ax,shrink=0.3)# label="Wind Speed (m/s)",
    else:
        cbar = fig.colorbar(speed_plot, ax=ax, shrink=0.5)#label="Wind Speed (m/s)",
                      
    cbar.ax.tick_params(labelsize='xx-small')
    cbar.set_label('Wind Speed (m/s)',size='xx-small')

    if (lon_max-lon_min)<60:
        alt_num = 1
    elif (lon_max-lon_min)>60 and (lon_max-lon_min)<200:
        alt_num = 2
    elif (lon_max-lon_min)>200:
        alt_num = 4
    Q = ax.quiver(lons[::alt_num,::alt_num], 
              lats[::alt_num,::alt_num], 
              ds_u[::alt_num,::alt_num], 
              ds_v[::alt_num,::alt_num],
              angles='xy', scale_units='xy',alpha=0.6)
    if ds_u.level.values<800:
        qk = ax.quiverkey(Q, 0.8, 0.9, 15, r'$15 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    else:
        qk = ax.quiverkey(Q, 0.8, 0.9, 5, r'$5 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    ax.set_xticks(np.linspace(speed_mean.lon.values.min(),#lon_min,
                              speed_mean.lon.values.max(),#lon_max,
                              num=5,endpoint=True))
    ax.set_yticks(np.linspace(speed_mean.lat.values.min(),#lon_min,
                              speed_mean.lat.values.max(),#lon_max,
                              num=5,endpoint=True))
    ax.set_xticklabels(np.linspace(speed_mean.lon.values.min(),#lon_min,
                                   speed_mean.lon.values.max(),#lon_max,
                                   num=5,endpoint=True),
                                   size='xx-small')
    ax.set_yticklabels(np.linspace(speed_mean.lat.values.min(),#lon_min,
                                   speed_mean.lat.values.max(),#lon_max,
                                   num=5,endpoint=True),
                                   size='xx-small')
    ax.set_xlabel('Longitude',size='x-small')
    ax.set_ylabel('Latitude',size='x-small')
    level_in_feet = {"1000.0": 'Surface',
                     "925.0": '3000 ft',
                     "850.0": '5000 ft',
                     "700.0": '10000 ft',
                     "500.0": '18000 ft'}
    ax.set_title('Month:'+calendar.month_name[time_s]+
                 '  Level:'+str(ds_u.level.values)+
                 ' '+ds_u.level.GRIB_name+' ('+level_in_feet[str(ds_u.level.values)]+')', 
                 size='x-small')
    ax.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax.transAxes)
    st.pyplot(fig)
#%% [markdown]
#Function to plot spatial variation in variables
def plot_spatial2(var_subset,lat_min, lat_max, lon_min, lon_max,time_s, region):
    #plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    x_size, y_size = calculate_x_y_size(lat_min, lat_max, lon_min, lon_max)
    fig, ax3 = plt.subplots(figsize=(x_size,y_size),
                            subplot_kw={'projection': ccrs.PlateCarree()})#,
    ax3.set_extent([var_subset.lon.values.min(),#lon_min,
                    var_subset.lon.values.max(),#lon_max, 
                    var_subset.lat.values.min(),#lat_min,
                    var_subset.lat.values.max()],#lat_max,], 
                    crs=ccrs.PlateCarree()) 
    India = load_shapefile()
    shape_feature = ShapelyFeature(India.geometry,
                                ccrs.PlateCarree(), edgecolor='black', 
                                facecolor='none')
    if (region != 'Global') and (region != 'India') and (region != 'China'):  
        ax3.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax3.add_feature(cfeature.COASTLINE)
    elif (region == 'India'):
        ax3.add_feature(shape_feature, linewidth=1.0, facecolor='none')
        #India.plot(facecolor='none',edgecolor='black',ax=ax3)
        ax3.add_feature(cfeature.COASTLINE)
    elif (region == 'China'):
        ax3.add_feature(shape_feature, linewidth=1.0, facecolor='none')
        ax3.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax3.add_feature(cfeature.COASTLINE)
    else:
        ax3.add_feature(cfeature.COASTLINE)
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(var_subset.lon, var_subset.lat)
    level_in_feet = {"1000.0": 'Surface',
                     "925.0": '3000 ft',
                     "850.0": '5000 ft',
                     "700.0": '10000 ft',
                     "500.0": '18000 ft'}

    if var_subset.var_desc=='Air temperature':
        plot_data = np.squeeze(var_subset.isel(time=time_s-1))-273.15
        s_plot = ax3.contourf(lons,lats,plot_data,
                              cmap=cmc.broc, 
                              vmin=-40, vamx=40,
                              levels=np.linspace(-40,40,81),
                              extend='both')
        if x_size<y_size:
            cbar = fig.colorbar(s_plot, ax=ax3,shrink=0.3)# label="2m Temperature (degC)",                  
        else:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)# label="2m Temperature (degC)", 
        cbar.set_label('2m Temperature (degC)',size='xx-small')
        ax3.set_title('Climatology in '+calendar.month_name[time_s]+
                 ' (1991-2021)', size='x-small')            
        
    elif var_subset.var_desc=='Precipitation':
        plot_data = np.squeeze(var_subset.isel(time=time_s-1))
        s_plot = ax3.contourf(lons,lats,plot_data,
                              cmap=cmc.batlowW_r, 
                              vmin=0, vamx=30,
                              levels=np.linspace(0,30,31),
                              extend='max')
        if x_size<y_size:
            cbar = fig.colorbar(s_plot, ax=ax3,shrink=0.3)# label="Mean Precipitation (mm/day)",                      
        else:
            cbar = fig.colorbar(s_plot, ax=ax3,shrink=0.5)# label="Mean Precipitation (mm/day)",                  
        cbar.set_label('Mean Precipitation (mm/day)',size='xx-small')
        ax3.set_title('Month:'+calendar.month_name[time_s]+
                 '  Level: Surface', size='x-small')
    else:
        plot_data = np.squeeze(var_subset.isel(time=time_s-1))
        s_plot = ax3.contourf(lons,lats,plot_data,
                              cmap=cmc.batlowW, 
                              vmin=0,vmax=100,
                              levels=np.linspace(0,100,26))
        if x_size<y_size:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.3)# label="Relative humidity (%)",                     
        else:
            cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)# label="Relative humidity (%)",                  
        cbar.set_label('Relative humidity (%)',size='xx-small')
        if str(var_subset.level.values) in level_in_feet:
            ax3.set_title('Month:'+calendar.month_name[time_s]+
                         '  Level:'+str(var_subset.level.values)+
                         ' '+var_subset.level.GRIB_name+' ('+level_in_feet[str(var_subset.level.values)]+')', 
                         size='medium')
        else:
            ax3.set_title('Month:'+calendar.month_name[time_s]+
                         '  Level:'+str(var_subset.level.values)+
                         ' '+var_subset.level.GRIB_name,
                         size='medium')  
    cbar.ax.tick_params(labelsize='xx-small')
    
    #ax3.set_title(calendar.month_name[time_s][:3],size='small')
    #ax3.set_title('Month:'+calendar.month_name[time_s][:3]+
    #             '  Level:'+str(var_subset.level.values)+
    #             ' '+var_subset.level.GRIB_name, size='x-small')
    
    

    ax3.set_xlabel('Longitude',size='x-small')
    ax3.set_ylabel('Latitude',size='x-small')

    ax3.set_xticks(np.linspace(np.floor(var_subset.lon.values.min()),#lon_min,
                               np.floor(var_subset.lon.values.max()),#lon_max,
                               num=5,endpoint=True))#format='%2.2f'))
    ax3.set_yticks(np.linspace(np.floor(var_subset.lat.values.min()),#lat_min,
                               np.floor(var_subset.lat.values.max()),#lat_max,
                               num=5,endpoint=True))#,format='%2.2f'))
    ax3.set_xticklabels(np.linspace(np.floor(var_subset.lon.values.min()),#lon_min,
                                    np.floor(var_subset.lon.values.max()),#lon_max,
                                    num=5,endpoint=True),#format='%2.2f'),
                                    size='xx-small')
    ax3.set_yticklabels(np.linspace(np.floor(var_subset.lat.values.min()),#lat_min,
                                    np.floor(var_subset.lat.values.max()),#lat_max,
                                    num=5,endpoint=True),#format='%2.2f'),
                                    size='xx-small')
    #if x_size<y_size:
    #    ax3.text(0.7,-0.4,'Data Source: '+ds_temp.attrs['source'],
    #             fontsize=4,transform=ax3.transAxes)
    #else:
    #    ax3.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],
    #             fontsize=4,transform=ax3.transAxes)
    st.pyplot(fig)   
#%% [markdown] 
# Function to plot wind speed and direction time series
def plot_time_series(speed,direction): 
    # plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    fig, ax1 = plt.subplots(figsize=(12,6))
    color = 'tab:blue'
    ax1.plot(np.arange(1,13), speed, label='Wind Speed',color=color)
    ax1.set_xlabel("Month")
    ax1.set_xticks(np.arange(1,13))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                        'Sep','Oct','Nov', 'Dec'])
    ax1.set_ylabel("Wind Speed (m/s)",color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    color = 'tab:red'
    ax2 = ax1.twinx()
    ax2.set_ylabel('Direction',color=color)
    ax2.scatter(np.arange(1,13), direction, label='Wind Direction', marker='o', color=color)
    ax2.set_ylim([0,360])
    ax2.set_yticks(np.arange(0,361,45))
    ax2.set_yticklabels(['N','NE','E','SE','S','SW','W','NW','N'],color=color)
    ax1.set_title('Latitude = '+str(speed_loc.lat.values)+
                  ' and Longitude = '+str(speed_loc.lon.values)+
                  '  Level:'+str(speed.level.values)+
                  ' '+speed.level.GRIB_name)
    ax1.text(0.7,-0.1,'Data Source: '+ds_temp.attrs['source'],
             fontsize=6,transform=ax1.transAxes)
    st.pyplot(fig)
#%% [markdown]
# Function to plot time series of variables
def plot_time_series2(var):
    # plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    fig, ax1 = plt.subplots(figsize=(12,6))

    color = 'tab:blue'
    if var.var_desc == 'Air temperature':        
        ax1.plot(np.arange(1,13), var-273.15,color=color)
        ax1.set_ylabel("2m Temperature (degC)",color=color, fontsize=14)
        ax1.set_yticks(np.linspace(np.floor(min(var.values - 273.15)),np.floor(max(var.values - 273.15)+1),num=5))
    elif var.var_desc == 'Precipitation':
        ax1.plot(np.arange(1,13), var,color=color)
        ax1.set_ylabel("Mean Precipitation (mm/day)",color=color, fontsize=14)
        ax1.set_yticks(np.linspace(np.floor(min(var.values)),np.floor(max(var.values)+1),num=5))
    else:
        ax1.plot(np.arange(1,13), var.values,color=color)
        ax1.set_xlabel("Month")
        ax1.set_ylabel("Relative humidity (%)",color=color, fontsize=14)
        ax1.set_yticks(np.linspace(np.floor(min(var.values)),
                                   np.floor(max(var.values)+1),num=5))
    
    ax1.set_xlabel("Month")
    ax1.set_xticks(np.arange(1,13))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                        'Sep','Oct','Nov', 'Dec'])
    ax1.text(0.7,-0.1,'Data Source: '+ds_temp.attrs['source'],
             fontsize=6,transform=ax1.transAxes)
    st.pyplot(fig)
#%% [markdown]
# Plotting vertical wind profile
def plot_vertical_wind(u_var, v_var, mon): #var has dim's level and month
    speed_loc,_ = calculate_wind(u_var,v_var)
    #speed_loc.reindex(level=list(reversed(speed_loc.level)))
    fig, ax1 = plt.subplots(figsize=(3,6))
    ax1.plot(speed_loc.values[mon,:],speed_loc.level.values)
    ax1.set_xlabel('Wind Speed (m/s)',size='x-small')
    ax1.set_ylabel('Vertical level (hPa)', size='x-small')
    ax1.set_xticks(np.linspace(np.floor(speed_loc.values.min()),
                               np.floor(speed_loc.values.max()),5,endpoint=True))
    ax1.set_xticklabels(np.linspace(np.floor(speed_loc.values.min()),
                               np.floor(speed_loc.values.max()),5,endpoint=True),size='x-small')
    ax1.set_yticks(speed_loc.level.values[::2])
    ax1.set_yticklabels(speed_loc.level.values[::2],size='x-small')
    ax1.invert_yaxis()
    ax1.set_title('Month: '+calendar.month_name[mon+1])
    #ax1.quiver(u_var.level.values,np.zeros(len(u_var.level.values)),
    #           u_var.values[mon,:],v_var.values[mon,:])
    st.pyplot(fig)
#%% [markdown]
# Plotting vertical wind profile
def plot_vertical_rh(var, mon): #var has dim's level and month
    fig, ax1 = plt.subplots(figsize=(3,6))
    ax1.plot(var.values[mon,:],var.level.values)
    ax1.set_xlabel('Relative Humidity (%)',size='x-small')
    ax1.set_ylabel('Vertical level (hPa)', size='x-small')
    ax1.set_xticks(np.linspace(np.floor(var.values.min()),
                               np.floor(var.values.max()),5,endpoint=True))
    ax1.set_xticklabels(np.linspace(np.floor(var.values.min()),
                               np.floor(var.values.max()),5,endpoint=True),
                               size='x-small')
    ax1.set_yticks(var.level.values[::2])
    ax1.set_yticklabels(var.level.values[::2],size='x-small')
    ax1.invert_yaxis()
    ax1.set_title('Month: '+calendar.month_name[mon+1])
    #ax1.quiver(u_var.level.values,np.zeros(len(u_var.level.values)),
    #           u_var.values[mon,:],v_var.values[mon,:])
    st.pyplot(fig)
#%% [markdown]
# Function to covert 0 360 to -180 to 180
def convert_180_180(ds_var):
    ds_var.coords['lon'] = (ds_var.coords['lon'] + 179.0625) % 358.125 - 179.0625
    ds_var2 = ds_var.sortby(ds_var.lon)
    return (ds_var2)
#%% [markdown]
# function to get the user input spatial box for plotting
def user_input_box(lat,lon):
    lat_min = st.sidebar.number_input("Enter Lat. min.", min_value=float(str(lat.values.min())), 
                                    max_value=float(str(lat.values.max())), 
                                    value=-84.00,step=0.01, format='%2.2f')
    lat_max = st.sidebar.number_input("Enter Lat. max.", float(str(lat.values.min())), 
                                    max_value=float(str(lat.values.max())), 
                                    value=84.00,step=0.01, format='%2.2f',
                                    placeholder="Must be greater than Lat min.")
    lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=float(str(lon.values.min())), 
                                    max_value=float(str(lon.values.max())), 
                                    value=-174.00,step=0.01, format='%3.2f')
    lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=float(str(lon.values.min())), 
                                    max_value=float(str(lon.values.max())), 
                                    value=174.00,step=0.01, format='%3.2f',
                                    placeholder="Must be greater than Lon min.")
    return lat_min, lat_max, lon_min, lon_max
#%% [markdown]
# function to find lat, lon location for plotting
def user_input_loc(lat,lon):
    lat_loc = st.sidebar.number_input("Enter Latitude", min_value=float(str(lat.values.min())), 
                                          max_value=float(str(lat.values.max())), 
                                          value=16.5,step=0.01, format='%2.2f')
    lon_loc = st.sidebar.number_input("Enter Longitude", min_value=float(str(lon.values.min())), 
                                          max_value=float(str(lon.values.max())), 
                                          value=78.5,step=0.01, format='%3.2f')
    return lat_loc, lon_loc
#%%
def user_input_region(region):
    df = pd.DataFrame({
        'Region': ['India', 'USA', 'Europe', 'Russia', 'China', 'Japan', 'Australia', 'Africa', 'South America'],
        'Lat_min': [0.0, 22.0, 35.0, 38.0, 18.0, 20.0, -48.0, -40.0, -58.0, ],
        'Lat_max': [40.0, 55.0, 72.0, 82.0, 54.0, 48.0, -5.0, 40.0, 15.0 ],
        'Lon_min': [65.0, -130.0, -25.0, -5.0, 73.0, 122.0, 105.0, -25.0, -95.0],
        'Lon_max': [100.0, -65.0, 65.0, 169.0, 135.0, 153.0, 160.0, 60.0, -25.0  ]
    })
    if region in df['Region'].values:
    # Get the latitude and longitude values for the selected region
        lat_min = df.loc[df['Region'] == region, 'Lat_min'].values[0]
        lat_max = df.loc[df['Region'] == region, 'Lat_max'].values[0]
        lon_min = df.loc[df['Region'] == region, 'Lon_min'].values[0]
        lon_max = df.loc[df['Region'] == region, 'Lon_max'].values[0]
    else:
        lat_min = -84.00
        lat_max = 84.00
        lon_min = -174.00
        lon_max = 174.00
    return lat_min, lat_max, lon_min, lon_max
#%% [markdown] 
## Streamlit App
#col1, col2 = st.columns(2, gap='small', vertical_alignment='center')
#with col1:
#st.image("My_page_enhanced.png", width=150)
#with col2:
#    st.title("K Narender Reddy", anchor='False')
#    st.write("Early Career Scientist, research interests include land surface modeling, crop modeling, and associated surface fluxes.")
#st.button('Go to Met. Visualisation', )
st.markdown('''**Details of author**  
        Built by: K Narender Reddy   
        Email :email: : knreddyiitd@gmail.com  
        Web Page :globe_with_meridians: : https://knreddy.online  
        Version 1: November, 2024''')
st.write("---")
st.title("Met. Data Visualization")
st.write("_NOTE: All data shown here is the Climatology data (1991-2021)_")
st.logo('icon.png',size='large')
#st.logo('My_page.png', location='right')
ds_temp = convert_180_180(load_temp_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_u = convert_180_180(load_uwind_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_v = convert_180_180(load_vwind_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_pr = convert_180_180(load_pr_data()).sel(lat=slice(-85,85),lon=slice(-176,176))
ds_pr.reindex(lat=list(reversed(ds_pr.lat)))
ds_rh = convert_180_180(load_rh_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
lon = ds_temp['lon']
lat = ds_temp['lat']

regions = ['Global', 'India', 'USA', 'Europe', 'Russia', 'China', 'Japan', 'Australia', 'Africa', 'South America'] 
#%% [markdown]
# User Inputs  
var_type = option_menu("Choose the variable", ("Temp_2m", 
                                                "Wind", 
                                                "Precipitation",
                                                "Relative Humidity"),
                                                menu_icon="None",#default_index='None', 
                                                orientation="horizontal",
                        icons=['thermometer', 'tornado', "cloud", 'droplet'],) 
                        #menu_icon="cast", default_index=0, orientation="horizontal")
#var_type = st.sidebar.selectbox("Choose the variable", ("Temp_2m", "Wind", "Precipitation","Relative Humidity"))
st.write("---")
if var_type == 'Wind':
    st.markdown("*Viewing Wind data*")
    st.markdown('''**Wind data can be viewed as**:  
        (1) windrose at a location and pressure level,  
        (2) wind vectors for a selected region,  
        (3) monthly time series of speed and direction, and  
        (4) vertical profile at a location
             ''')
    st.write("select your choice of plot from the side bar:")
    plot_type = option_menu("", ("Wind Rose",
                                "Spatial Wind Vectors", 
                                "Time Series",
                                "Vertical Profile"), orientation="horizontal",
                                icons=['test','test','test','test'])
    #plot_type = st.sidebar.selectbox("Choose Plot Type", ("Wind Rose",
    #                                                      "Spatial Wind Vectors", 
    #                                                      "Time Series",
    #                                                      "Vertical Profile"))
    if plot_type == "Wind Rose":
        st.header("Wind Rose")
        st.write("Default location, and pressure level is shown here. Please select your region of interest using latitude and longitude and pressure level")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)

        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                          ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))

        plot_wind_rose(speed_loc.sel(level=level_sel).values, 
                       direction_loc.sel(level=level_sel).values,
                       speed_loc.lat.values,speed_loc.lon.values)

    elif plot_type == "Spatial Wind Vectors":
        st.header("Spatial Wind Vectors")    
        st.write("Default region is shown here. Please select your region of interest.")    
        region_sel = st.sidebar.selectbox('Select the region of interest', regions)
        [lat_min, lat_max, lon_min, lon_max] = user_input_region(region_sel)
        #[lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        
        ds_u_subset = ds_u['uwnd'].sel(lat=slice(lat_max, lat_min), 
                                       lon=slice(lon_min, lon_max),
                                       level=level_sel)
        ds_v_subset = ds_v['vwnd'].sel(lat=slice(lat_max, lat_min), 
                                       lon=slice(lon_min, lon_max),
                                       level=level_sel)
        plot_wind_vectors(ds_u_subset[time_sel-1,:,:], ds_v_subset[time_sel-1,:,:], 
                          lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)

    elif plot_type == "Time Series":
        st.header("Wind Speed Time Series")
        st.write("Default location and pressure level is shown here. Please select your location of interest using latitude and longitude, and the pressure level")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                  ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))
        plot_time_series(speed_loc.sel(level=level_sel),direction_loc.sel(level=level_sel))
    elif plot_type == 'Vertical Profile': 
        st.header("Wind Speed Vertical Profile")
        st.write("Default location and month is shown here. Please select your location of interest using latitude and longitude, and desired month")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        mon_sel = st.sidebar.selectbox("Select Month",np.arange(1,13))
        ds_u_loc = ds_u['uwnd'].sel(lat=lat_loc, lon=lon_loc, method='nearest')
        ds_v_loc = ds_v['vwnd'].sel(lat=lat_loc, lon=lon_loc, method='nearest')
        plot_vertical_wind(ds_u_loc,ds_v_loc,mon_sel-1)
########################### 2m Temperature
elif var_type == 'Temp_2m':
    st.markdown('''2m Temperature can be viewed as:  
                 (1) spatial plot for a selected region  
                 (2) monthly time series at a location
                 ''')
    st.write("select your choice of plot from the side bar:")
    #plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    plot_type = option_menu("", ("Spatial plot", 
                                "Time Series"), orientation="horizontal",
                                icons=['test','test'])
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")
        st.write("Default region is shown here. Please select your region of interest using latitude and longitude")       
        region_sel = st.sidebar.selectbox('Select the region of interest', regions)
        [lat_min, lat_max, lon_min, lon_max] = user_input_region(region_sel)
        #[lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        ds_temp_subset = ds_temp['air'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))

        
        anomaly_plt = st.sidebar.checkbox('Check Anomaly for the month selected \n in the year 2024 or 2025')
        if anomaly_plt:
            col1, col2 = st.columns(2)
            with col1:
                plot_spatial2(ds_temp_subset.sel(level=2), lat_min, lat_max, 
                      lon_min, lon_max,time_sel, region_sel)
            with col2:
                anomaly_plotting(ds_temp_subset.sel(level=2),time_sel, region_sel)
        else:
            plot_spatial2(ds_temp_subset.sel(level=2), lat_min, lat_max, 
                      lon_min, lon_max,time_sel, region_sel)
    else:
        st.header("Time series plot")
        st.write("Default location is shown here. Please select your location of interest using latitude and longitude") 
        lat_loc, lon_loc = user_input_loc(lat,lon)
        ds_air = ds_temp['air']
        lat_idx = (np.abs(lat-lat_loc)).argmin()
        lon_idx = (np.abs(lon-lon_loc)).argmin()
        temp_loc = ds_air.isel(level=0,lat=lat_idx,lon=lon_idx)
        plot_time_series2(temp_loc)
##################################### Precipitation        
elif var_type == 'Precipitation':
    st.markdown('''Monthly mean rainfall can be viewed as:  
                (1) spatial plot for a selected region  
                (2) a monthly time series at a loctaion
                ''')
    st.write("select your choice of plot from the side bar:")
    #plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    plot_type = option_menu("", ("Spatial plot", 
                                "Time Series"), orientation="horizontal",
                                icons=['test','test'])
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")
        st.write("Default region is shown here. Please select your region of interest.")
        region_sel = st.sidebar.selectbox('Select the region of interest', regions)
        [lat_min, lat_max, lon_min, lon_max] = user_input_region(region_sel)
        #[lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)        
        ds_pr_subset = ds_pr['precip'].sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        plot_spatial2(ds_pr_subset, lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
    else:
        st.header("Time series plot")
        st.write("Default location is shown here. Please select your location of interest using latitude and longitude") 
        lat_loc, lon_loc = user_input_loc(lat,lon)
        pr_loc = ds_pr['precip'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        plot_time_series2(pr_loc)
######################################## Relative Humidity
else: #Relative Humidity
    st.markdown('''Relative Humidity can be viewed as:  
                 (1) spatial plot for a selected region  
                 (2) monthly time series at a location  
                 (3) vertical profile at a location
                 ''')
    #plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", 
    #                                                      "Time Series",
    #                                                      "Vertical Profile"))
    plot_type = option_menu("", ("Spatial plot", 
                                "Time Series",
                                "Vertical Profile"), orientation="horizontal",
                                icons=['test','test','test'])
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")   
        st.write("Default region is shown here. Please select your region of interest.")     
        region_sel = st.sidebar.selectbox('Select the region of interest', regions)
        [lat_min, lat_max, lon_min, lon_max] = user_input_region(region_sel)
        #[lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        ds_rh_subset = ds_rh['rhum'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)

        plot_spatial2(ds_rh_subset.sel(level=level_sel), lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
    elif plot_type == 'Time Series':
        st.header("Time series plot")
        st.write("Default location and pressure level is shown here. Please select your location of interest using latitude and longitude, and the pressure level") 
        lat_loc, lon_loc = user_input_loc(lat,lon)
        rh_loc = ds_rh['rhum'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)
        plot_time_series2(rh_loc.sel(level=level_sel))
    elif plot_type == 'Vertical Profile':
        st.header("Relative Humidity Vertical Profile")
        st.write("Default location and month is shown here. Please select your location of interest using latitude and longitude, and desired month")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        mon_sel = st.sidebar.selectbox("Select Month",np.arange(1,13))
        rh_loc = ds_rh['rhum'].sel(lat=lat_loc, lon=lon_loc, method='nearest')
        plot_vertical_rh(rh_loc,mon_sel-1)  
# %%
# %% [markdown]
st.write("---")
st.markdown(f"Data Source: {ds_temp.attrs['source']}")
