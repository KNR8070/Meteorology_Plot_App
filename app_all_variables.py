#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 20:54:39 2024

@author: knreddy
"""
#%% [markdown]
## load modules

import glob
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
import plotly.graph_objects as go
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

def _find_era5_file(var_name):
    """Return the path to the most recent era5_*_<var_name>.nc file in data/."""
    files = sorted(glob.glob(f'data/era5_*_{var_name}.nc'))
    return files[-1] if files else f'data/era5_2024_2025_{var_name}.nc'


def _add_map_features(ax, region):
    if region in ('India', 'China'):
        India = load_shapefile()
        shape_feature = ShapelyFeature(India.geometry, ccrs.PlateCarree(),
                                       edgecolor='black', facecolor='none')
        ax.add_feature(shape_feature, linewidth=1.0, facecolor='none')
        if region == 'China':
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE)
    elif region == 'Global':
        ax.add_feature(cfeature.COASTLINE)
    else:
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE)


# %% Anomaly plotting
def anomaly_plotting(clim_var,mon, region):
    [lat_min, lat_max, lon_min, lon_max] = user_input_region(region)
    data_path = _find_era5_file(clim_var.name)

    try:
        var_data = xr.open_dataset(data_path)
    except FileNotFoundError:
        st.error(f"Current data file not found: {data_path}. Cannot compute anomaly.")
        return

    var_data = convert_180_180(var_data).sel(lat=slice(85,-85), lon=slice(-176,176))
    var_subset = var_data[clim_var.name].sel(lat=slice(lat_max, lat_min), 
                                             lon=slice(lon_min, lon_max))

    if 'level' in var_subset.dims:
        if 1000 in var_subset.level.values:
            var_subset = var_subset.sel(level=1000)
        else:
            var_subset = var_subset.isel(level=0)

    # Sort N→S so slice(lat_max, lat_min) works regardless of source lat order
    # (NCEP precip is S→N; temperature and wind are N→S)
    clim_var_ns = clim_var.sortby('lat', ascending=False)
    clim_var_subset = clim_var_ns.sel(lat=slice(lat_max, lat_min),
                                      lon=slice(lon_min, lon_max))

    var_plot_data, year = select_current_year_data(var_subset, mon)
    if var_plot_data is None:
        st.error(f"No data found for month {mon} in current dataset.")
        return

    # Interpolate climatology to the ERA5 grid (2.5° NCEP → 0.25° ERA5).
    # Use .values to pass bare numpy arrays — avoids a time-coordinate type conflict
    # (ERA5 uses datetime64 while NCEP climatology uses cftime.DatetimeGregorian).
    clim_interp = clim_var_subset.interp(lat=var_plot_data.lat.values, lon=var_plot_data.lon.values, method='linear')

    if clim_var.var_desc == 'Air temperature':
        clim_plot_data = np.squeeze(clim_interp.isel(time=mon-1))
        plot_data = var_plot_data - clim_plot_data  # both in K; difference equals °C anomaly
        cmap = cmc.vik
        vmin, vmax = -20, 20
        plot_label = '2m Temperature Anomaly (°C)'
    elif clim_var.var_desc == 'Precipitation':
        clim_plot_data = np.squeeze(clim_interp.isel(time=mon-1))
        plot_data = var_plot_data - clim_plot_data
        cmap = cmc.broc
        vmin, vmax = -15, 15
        plot_label = 'Precipitation Anomaly (mm/day)'
    elif clim_var.var_desc == 'Relative humidity':
        clim_plot_data = np.squeeze(clim_interp.isel(time=mon-1))
        plot_data = var_plot_data - clim_plot_data
        cmap = cmc.vik
        vmin, vmax = -30, 30
        plot_label = 'Relative Humidity Anomaly (%)'
    else:
        st.error('Anomaly plotting is available only for temperature, precipitation, and humidity.')
        return

###### Plotting
    x_size, y_size = calculate_x_y_size(lat_min, lat_max, lon_min, lon_max)
    fig, ax3 = plt.subplots(figsize=(x_size,y_size),
                            subplot_kw={'projection': ccrs.PlateCarree()})
    ax3.set_extent([var_subset.lon.values.min(),
                    var_subset.lon.values.max(),
                    var_subset.lat.values.min(),
                    var_subset.lat.values.max()],
                    crs=ccrs.PlateCarree())
    _add_map_features(ax3, region)

    # plot_data is on the ERA5 grid after interpolation — use its coords for the meshgrid
    lons, lats = np.meshgrid(var_plot_data.lon, var_plot_data.lat)
    s_plot = ax3.contourf(lons, lats, plot_data,
                          cmap=cmap,
                          vmin=vmin, vmax=vmax,
                          levels=np.linspace(vmin, vmax, 41),
                          extend='both')

    if x_size < y_size:
        cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.3)
    else:
        cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)
    cbar.set_label(plot_label, size='xx-small')
    ax3.set_title('Anomaly in '+calendar.month_name[mon]+' '+ str(year)+' compared to 1991-2021 clim.', size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    ax3.set_xlabel('Longitude',size='x-small')
    ax3.set_ylabel('Latitude',size='x-small')

    ax3.set_xticks(np.linspace(np.floor(var_plot_data.lon.values.min()),
                               np.floor(var_plot_data.lon.values.max()),
                               num=5,endpoint=True))
    ax3.set_yticks(np.linspace(np.floor(var_plot_data.lat.values.min()),
                               np.floor(var_plot_data.lat.values.max()),
                               num=5,endpoint=True))
    ax3.set_xticklabels(np.linspace(np.floor(var_plot_data.lon.values.min()),
                                    np.floor(var_plot_data.lon.values.max()),
                                    num=5,endpoint=True),
                                    size='xx-small')
    ax3.set_yticklabels(np.linspace(np.floor(var_plot_data.lat.values.min()),
                                    np.floor(var_plot_data.lat.values.max()),
                                    num=5,endpoint=True),
                                    size='xx-small')
    st.pyplot(fig)


def select_current_year_data(var_subset, mon):
    month_subset = var_subset.sel(time=var_subset.time.dt.month == mon)
    if len(month_subset.time) == 0:
        return None, None
    years = np.unique(month_subset.time.dt.year.values)
    year = int(years[-1])
    year_data = month_subset.sel(time=month_subset.time.dt.year == year)
    return year_data.isel(time=0), year


def plot_wind_speed_anomaly(ds_u_clim, ds_v_clim, lat_min, lat_max, lon_min, lon_max, level_sel, mon, region):
    data_path = _find_era5_file('wnd')
    try:
        ds_current = xr.open_dataset(data_path)
    except FileNotFoundError:
        st.error(f"Current wind data file not found: {data_path}. Cannot compute wind anomaly.")
        return

    ds_current = convert_180_180(ds_current).sel(lat=slice(85,-85), lon=slice(-176,176))
    u_current = ds_current['uwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max), level=level_sel)
    v_current = ds_current['vwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max), level=level_sel)

    u_current_sel, year = select_current_year_data(u_current, mon)
    v_current_sel, _ = select_current_year_data(v_current, mon)
    if u_current_sel is None or v_current_sel is None:
        st.error(f"No wind data found for month {mon} in current dataset.")
        return

    current_speed = np.sqrt(u_current_sel**2 + v_current_sel**2)
    clim_u = ds_u_clim['uwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max), level=level_sel).isel(time=mon-1)
    clim_v = ds_v_clim['vwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max), level=level_sel).isel(time=mon-1)
    clim_speed = np.sqrt(clim_u**2 + clim_v**2)
    clim_speed_interp = clim_speed.interp(lat=current_speed.lat.values, lon=current_speed.lon.values, method='linear')

    plot_data = current_speed - clim_speed_interp

    x_size, y_size = calculate_x_y_size(lat_min, lat_max, lon_min, lon_max)
    fig, ax3 = plt.subplots(figsize=(x_size,y_size), subplot_kw={'projection': ccrs.PlateCarree()})
    ax3.set_extent([plot_data.lon.values.min(), plot_data.lon.values.max(), plot_data.lat.values.min(), plot_data.lat.values.max()], crs=ccrs.PlateCarree())
    _add_map_features(ax3, region)

    lons, lats = np.meshgrid(plot_data.lon, plot_data.lat)
    s_plot = ax3.contourf(lons, lats, plot_data, cmap=cmc.vik, vmin=-10, vmax=10, levels=np.linspace(-10, 10, 41), extend='both')
    if x_size < y_size:
        cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.3)
    else:
        cbar = fig.colorbar(s_plot, ax=ax3, shrink=0.5)
    cbar.set_label('Wind Speed Anomaly (m/s)', size='xx-small')
    ax3.set_title('Wind speed anomaly in '+calendar.month_name[mon]+' '+ str(year)+' compared to 1991-2021 clim.', size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    ax3.set_xlabel('Longitude', size='x-small')
    ax3.set_ylabel('Latitude', size='x-small')
    ax3.set_xticks(np.linspace(np.floor(plot_data.lon.values.min()), np.floor(plot_data.lon.values.max()), num=5, endpoint=True))
    ax3.set_yticks(np.linspace(np.floor(plot_data.lat.values.min()), np.floor(plot_data.lat.values.max()), num=5, endpoint=True))
    ax3.set_xticklabels(np.linspace(np.floor(plot_data.lon.values.min()), np.floor(plot_data.lon.values.max()), num=5, endpoint=True), size='xx-small')
    ax3.set_yticklabels(np.linspace(np.floor(plot_data.lat.values.min()), np.floor(plot_data.lat.values.max()), num=5, endpoint=True), size='xx-small')
    st.pyplot(fig)

def _load_era5_point(var_name, lat_loc, lon_loc):
    """Load ERA5 file and return a 1-D time series at the nearest grid point.
    Uses isel+argmin to avoid InvalidIndexError from duplicate coords after lon wrapping."""
    data_path = _find_era5_file(var_name)
    try:
        ds = xr.open_dataset(data_path)
    except FileNotFoundError:
        return None, data_path
    ds = convert_180_180(ds).sel(lat=slice(85, -85), lon=slice(-176, 176))
    arr = ds[var_name]
    lat_idx = int(np.abs(arr.lat - lat_loc).argmin())
    lon_idx = int(np.abs(arr.lon - lon_loc).argmin())
    return arr.isel(lat=lat_idx, lon=lon_idx), None


def plot_time_series_with_anomaly(clim_var_loc, lat_loc, lon_loc):
    """Climatology time series (left axis) overlaid with ERA5 anomaly bars (right axis)."""
    var_loc, missing = _load_era5_point(clim_var_loc.name, lat_loc, lon_loc)
    if var_loc is None:
        st.warning(f"ERA5 file not found ({missing}). Showing climatology only.")
        plot_time_series2(clim_var_loc)
        return

    if 'level' in var_loc.dims:
        var_loc = var_loc.sel(level=1000) if 1000 in var_loc.level.values else var_loc.isel(level=0)

    months = np.arange(1, 13)
    anomalies, year_labels = [], []
    for mon in months:
        era5_val, year = select_current_year_data(var_loc, mon)
        if era5_val is None:
            anomalies.append(np.nan)
            year_labels.append(None)
        else:
            anomalies.append(float(era5_val.values) - float(clim_var_loc.isel(time=mon - 1).values))
            year_labels.append(year)

    year_str = '/'.join(str(y) for y in sorted(set(y for y in year_labels if y is not None)))

    if clim_var_loc.var_desc == 'Air temperature':
        clim_vals = clim_var_loc.values - 273.15
        clim_ylabel = '2m Temperature (°C)'
        anom_ylabel = 'Anomaly (°C)'
    elif clim_var_loc.var_desc == 'Precipitation':
        clim_vals = clim_var_loc.values
        clim_ylabel = 'Precipitation (mm/day)'
        anom_ylabel = 'Anomaly (mm/day)'
    else:
        clim_vals = clim_var_loc.values
        clim_ylabel = 'Relative Humidity (%)'
        anom_ylabel = 'Anomaly (%)'

    fig, ax1 = plt.subplots(figsize=(12, 5))
    clim_color = 'tab:blue'
    ax1.plot(months, clim_vals, color=clim_color, linewidth=2, marker='o', markersize=4,
             label='Climatology (1991-2021)')
    ax1.set_xlabel('Month')
    ax1.set_ylabel(clim_ylabel, color=clim_color, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=clim_color)
    ax1.set_xticks(months)
    ax1.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])

    ax2 = ax1.twinx()
    bar_colors = ['#d73027' if (not np.isnan(a) and a > 0) else '#4575b4' for a in anomalies]
    ax2.bar(months, np.nan_to_num(anomalies), color=bar_colors, alpha=0.55,
            edgecolor='black', linewidth=0.5, label=f'{year_str} Anomaly')
    ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax2.set_ylabel(anom_ylabel, fontsize=12)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
    ax1.set_title(f'Climatology & anomaly at ({lat_loc:.1f}°N, {lon_loc:.1f}°E) — {year_str} vs 1991-2021')
    ax1.text(0.7, -0.1, 'Data Source: ' + ds_temp.attrs['source'], fontsize=6, transform=ax1.transAxes)
    st.pyplot(fig)


def plot_wind_time_series_with_anomaly(speed_clim, direction_clim, lat_loc, lon_loc, level_sel):
    """Wind speed climatology (left axis) overlaid with ERA5 wind speed anomaly bars (right axis)."""
    data_path = _find_era5_file('wnd')
    try:
        ds_current = xr.open_dataset(data_path)
    except FileNotFoundError:
        st.warning(f"ERA5 wind file not found. Showing climatology only.")
        plot_time_series(speed_clim, direction_clim)
        return

    ds_current = convert_180_180(ds_current).sel(lat=slice(85, -85), lon=slice(-176, 176))
    lat_idx = int(np.abs(ds_current.lat - lat_loc).argmin())
    lon_idx = int(np.abs(ds_current.lon - lon_loc).argmin())
    u_loc = ds_current['uwnd'].isel(lat=lat_idx, lon=lon_idx).sel(level=level_sel)
    v_loc = ds_current['vwnd'].isel(lat=lat_idx, lon=lon_idx).sel(level=level_sel)

    months = np.arange(1, 13)
    anomalies, year_labels = [], []
    for mon in months:
        u_val, year = select_current_year_data(u_loc, mon)
        v_val, _ = select_current_year_data(v_loc, mon)
        if u_val is None or v_val is None:
            anomalies.append(np.nan)
            year_labels.append(None)
        else:
            era5_spd = float(np.sqrt(u_val**2 + v_val**2))
            clim_spd = float(speed_clim.isel(time=mon - 1).values)
            anomalies.append(era5_spd - clim_spd)
            year_labels.append(year)

    year_str = '/'.join(str(y) for y in sorted(set(y for y in year_labels if y is not None)))

    fig, ax1 = plt.subplots(figsize=(12, 5))
    clim_color = 'tab:blue'
    ax1.plot(months, speed_clim.values, color=clim_color, linewidth=2, marker='o', markersize=4,
             label='Climatology (1991-2021)')
    ax1.set_xlabel('Month')
    ax1.set_ylabel('Wind Speed (m/s)', color=clim_color, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=clim_color)
    ax1.set_xticks(months)
    ax1.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])

    ax2 = ax1.twinx()
    bar_colors = ['#d73027' if (not np.isnan(a) and a > 0) else '#4575b4' for a in anomalies]
    ax2.bar(months, np.nan_to_num(anomalies), color=bar_colors, alpha=0.55,
            edgecolor='black', linewidth=0.5, label=f'{year_str} Anomaly')
    ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax2.set_ylabel('Anomaly (m/s)', fontsize=12)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
    ax1.set_title(f'Wind speed climatology & anomaly at ({lat_loc:.1f}°N, {lon_loc:.1f}°E), {level_sel} hPa — {year_str} vs 1991-2021')
    st.pyplot(fig)


#%% [markdown]
# Calculate wind speed and direction from u10 and v10
def calculate_wind(uwnd,vwnd):
    speed = np.sqrt(uwnd**2 + vwnd**2)
    direction = (np.arctan2(vwnd, uwnd) * 180 / np.pi + 180) % 360
    return speed, direction
#%% [markdown]
# Function to plot windrose Plot
def plot_wind_rose(speed_pwr, direction_pwr,lat_l,lon_l):
    fig_ws = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax(fig=fig_ws)
    ax.bar(direction_pwr, speed_pwr, normed=True, opening=0.8, edgecolor='white')
    ax.set_title('Latitude = '+str(lat_l)+' and Longitude = '+str(lon_l))
    ax.set_legend(title="Wind Speed (m/s)",loc='best')
    st.pyplot(fig_ws)
#%% [markdown]
# calculating figsize for spatial plots
def calculate_x_y_size(lat_min, lat_max, lon_min, lon_max):
    ### Dynamic plot size
    lat_range = lat_max-lat_min
    lon_range = lon_max-lon_min
    
    x_ratio = lon_range/(lon_range+lat_range)
    y_ratio = lat_range/(lon_range+lat_range)
    
    if x_ratio > y_ratio:
        x_size = 7
        y_size = 5
    elif x_ratio == y_ratio:
        x_size = 6
        y_size = 6
    else:
        x_size = 5
        y_size = 7
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
    _add_map_features(ax, region)
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
    #ax.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax.transAxes)
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
    _add_map_features(ax3, region)
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
                              vmin=-40, vmax=40,
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
                              vmin=0, vmax=30,
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
    ax1.set_title('Latitude = '+str(speed.lat.values)+
                  ' and Longitude = '+str(speed.lon.values)+
                  '  Level:'+str(speed.level.values)+
                  ' '+speed.level.GRIB_name)
    #ax1.text(0.7,-0.1,'Data Source: '+ds_temp.attrs['source'],
    #         fontsize=6,transform=ax1.transAxes)
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
    # New CDS API (cdsapi >= 0.7) names the time dimension 'valid_time'; normalise to 'time'
    if 'valid_time' in ds_var.dims and 'time' not in ds_var.dims:
        ds_var = ds_var.rename({'valid_time': 'time'})
    if float(ds_var.lon.values.min()) >= 0:  # only convert 0-360 data; ERA5 is already -180/180
        ds_var.coords['lon'] = (ds_var.coords['lon'] + 180) % 360 - 180
    ds_var2 = ds_var.sortby(ds_var.lon)
    # Drop any duplicate lon values that can arise near the ±180 boundary
    _, unique_idx = np.unique(ds_var2.lon.values, return_index=True)
    if len(unique_idx) < len(ds_var2.lon):
        ds_var2 = ds_var2.isel(lon=sorted(unique_idx.tolist()))
    return ds_var2

# ── Globe plotting helpers ─────────────────────────────────────────────────────
_G_CMAPS      = {"air": "RdBu_r", "wnd": "Viridis",  "precip": "BrBG",   "rhum": "RdYlGn"}
_G_ANOM_CMAPS = {"air": "RdBu",   "wnd": "RdBu",     "precip": "BrBG",   "rhum": "RdYlGn"}
_G_CLIM_RANGE = {"air": (-30, 40), "wnd": (0, 20),   "precip": (0, 15),  "rhum": (0, 100)}
_G_ANOM_RANGE = {"air": (-10, 10), "wnd": (-5, 5),   "precip": (-8, 8),  "rhum": (-20, 20)}
_G_UNITS      = {"air": "°C",      "wnd": "m/s",     "precip": "mm/day", "rhum": "%"}


@st.cache_data
def _globe_coastlines():
    """Return NaN-separated x,y,z arrays for all Natural Earth coastlines on a unit sphere."""
    segments = []
    for geom in cfeature.COASTLINE.geometries():
        lines = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        for line in lines:
            coords = np.array(line.coords)
            if len(coords) < 2:
                continue
            lr    = np.radians(coords[:, 0])
            latr  = np.radians(coords[:, 1])
            theta = np.pi / 2 - latr
            x = 1.005 * np.sin(theta) * np.cos(lr)
            y = 1.005 * np.sin(theta) * np.sin(lr)
            z = 1.005 * np.cos(theta)
            segments.append(np.column_stack([x, y, z]))
    nan_row  = np.full((1, 3), np.nan)
    combined = np.vstack([np.vstack([s, nan_row]) for s in segments])
    return combined[:, 0], combined[:, 1], combined[:, 2]


def _da_to_sphere_mesh(da):
    """Convert a 2-D lat/lon DataArray to Cartesian sphere coordinates and surface values."""
    da = da.sortby("lat").sortby("lon")
    lon2d, lat2d = np.meshgrid(da.lon.values.astype(float), da.lat.values.astype(float))
    theta = np.radians(90.0 - lat2d)
    phi   = np.radians(lon2d)
    return (np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
            da.values.astype(float),
            lon2d, lat2d)


def _build_globe(da, title, cmap, vmin, vmax, units, height=450):
    """Return an interactive Plotly 3-D globe figure with auto-rotate animation."""
    x, y, z, surf, lon2d, lat2d = _da_to_sphere_mesh(da)
    cx, cy, cz = _globe_coastlines()
    fig = go.Figure([
        go.Surface(
            x=x, y=y, z=z,
            surfacecolor=surf,
            colorscale=cmap,
            cmin=vmin, cmax=vmax,
            colorbar=dict(title=dict(text=units, side="right"), thickness=12, len=0.5),
            lighting=dict(ambient=0.8, diffuse=0.4, specular=0.05),
            hovertemplate=(
                "lon: %{customdata[0]:.1f}°<br>"
                "lat: %{customdata[1]:.1f}°<br>"
                f"value: %{{customdata[2]:.2f}} {units}<extra></extra>"
            ),
            customdata=np.dstack([lon2d, lat2d, surf]),
        ),
        go.Scatter3d(
            x=cx, y=cy, z=cz,
            mode="lines",
            line=dict(color="white", width=0.8),
            showlegend=False, hoverinfo="none",
        ),
    ])

    # Animation frames: rotate camera eye 360° around the vertical axis
    n_frames = 72                       # 5° per step — smooth full rotation
    r_eye    = np.sqrt(1.4**2 + 0.5**2)  # keep same radial distance
    z_eye    = 0.5
    fig.frames = [
        go.Frame(
            layout=dict(scene_camera=dict(eye=dict(
                x=float(r_eye * np.cos(2 * np.pi * i / n_frames)),
                y=float(r_eye * np.sin(2 * np.pi * i / n_frames)),
                z=z_eye,
            ))),
            name=str(i),
        )
        for i in range(n_frames)
    ]

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=12, color="white")),
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            bgcolor="black", aspectmode="cube",
            camera=dict(eye=dict(
                x=float(r_eye * np.cos(0)),
                y=float(r_eye * np.sin(0)),
                z=z_eye,
            )),
        ),
        margin=dict(l=0, r=0, t=35, b=40),
        paper_bgcolor="black",
        height=height,
        updatemenus=[dict(
            type="buttons",
            showactive=False,
            direction="left",
            x=0.5, xanchor="center",
            y=0.0, yanchor="bottom",
            bgcolor="rgba(30,30,30,0.85)",
            bordercolor="rgba(120,120,120,0.5)",
            font=dict(color="white", size=12),
            buttons=[
                dict(
                    label="▶  Rotate",
                    method="animate",
                    args=[None, dict(
                        frame=dict(duration=50, redraw=True),
                        fromcurrent=True,
                        mode="immediate",
                        transition=dict(duration=0),
                        loop=True,
                    )],
                ),
                dict(
                    label="⏸  Pause",
                    method="animate",
                    args=[[None], dict(
                        frame=dict(duration=0, redraw=False),
                        mode="immediate",
                        transition=dict(duration=0),
                    )],
                ),
            ],
        )],
    )
    return fig


def plot_globe_comparison(var_key, month, level=None):
    """Render three interactive globe panels: climatology | ERA5 current | anomaly."""
    month_name = calendar.month_name[month]

    # climatology DataArray (uses module-level ds_* loaded at startup)
    if var_key == "air":
        da_clim = ds_temp["air"].isel(time=month - 1) - 273.15
        da_clim.name = "air"
    elif var_key == "wnd":
        u = ds_u["uwnd"].sel(level=level).isel(time=month - 1)
        v = ds_v["vwnd"].sel(level=level).isel(time=month - 1)
        da_clim = np.sqrt(u ** 2 + v ** 2)
        da_clim.name = "wnd"
    elif var_key == "precip":
        da_clim = ds_pr["precip"].sortby("lat").isel(time=month - 1)
    elif var_key == "rhum":
        da_clim = ds_rh["rhum"].sel(level=level).isel(time=month - 1)
    da_clim = da_clim.squeeze()

    # ERA5 DataArray
    path = _find_era5_file(var_key)
    try:
        ds_e5 = convert_180_180(xr.open_dataset(path)).sel(
            lat=slice(85, -85), lon=slice(-176, 176)
        )
    except (FileNotFoundError, OSError):
        st.error(f"No ERA5 file found for '{var_key}'. Run scripts/fetch_era5.py first.")
        return
    if var_key == "wnd":
        u5 = ds_e5["uwnd"].sel(level=level) if level is not None else ds_e5["uwnd"].isel(level=0)
        v5 = ds_e5["vwnd"].sel(level=level) if level is not None else ds_e5["vwnd"].isel(level=0)
        da_e5 = np.sqrt(u5 ** 2 + v5 ** 2)
        da_e5.name = "wnd"
    else:
        da_e5 = ds_e5[var_key]
        if "level" in da_e5.dims:
            da_e5 = da_e5.sel(level=level) if level is not None else da_e5.isel(level=0)
        if var_key == "air":
            da_e5 = da_e5 - 273.15
    da_e5_pt, year = select_current_year_data(da_e5, month)
    if da_e5_pt is None:
        st.error(f"No ERA5 data found for month {month}. Cannot show globe view.")
        return
    da_e5_pt = da_e5_pt.squeeze()

    # anomaly: interpolate climatology onto ERA5 grid
    da_clim_i = da_clim.sortby("lat").interp(
        lat=np.sort(da_e5_pt.lat.values),
        lon=np.sort(da_e5_pt.lon.values),
        method="linear",
    )
    anom = da_e5_pt.sortby(["lat", "lon"]) - da_clim_i

    units    = _G_UNITS[var_key]
    vmin_c, vmax_c = _G_CLIM_RANGE[var_key]
    vmin_a, vmax_a = _G_ANOM_RANGE[var_key]

    with st.spinner("Building globes …"):
        c1, c2, c3 = st.columns(3)
        with c1:
            st.plotly_chart(
                _build_globe(da_clim, f"Climatology – {month_name}",
                             _G_CMAPS[var_key], vmin_c, vmax_c, units),
                use_container_width=True)
            st.caption("NCEP 1991–2020 climatology")
        with c2:
            st.plotly_chart(
                _build_globe(da_e5_pt, f"ERA5 – {month_name} {year}",
                             _G_CMAPS[var_key], vmin_c, vmax_c, units),
                use_container_width=True)
            st.caption(f"ERA5 · {month_name} {year}")
        with c3:
            st.plotly_chart(
                _build_globe(anom, f"Anomaly – {month_name} {year}",
                             _G_ANOM_CMAPS[var_key], vmin_a, vmax_a, units),
                use_container_width=True)
            st.caption(f"ERA5 {year} − Climatology")

# ── end of globe helpers ───────────────────────────────────────────────────────
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
st.title("Met. Data Visualization")
st.write("_NOTE: All data shown here is the Climatology data (1991-2021)_")
st.logo('icon.png',size='large')
#st.logo('My_page.png', location='right')
ds_temp = convert_180_180(load_temp_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_u = convert_180_180(load_uwind_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_v = convert_180_180(load_vwind_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
ds_pr = convert_180_180(load_pr_data()).sel(lat=slice(-85,85),lon=slice(-176,176))
ds_rh = convert_180_180(load_rh_data()).sel(lat=slice(85,-85),lon=slice(-176,176))
lon = ds_temp['lon']
lat = ds_temp['lat']

anomaly_plt = False
globe_plt = False

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
        anomaly_plt = st.sidebar.checkbox('Show wind speed anomaly for selected region and month')
        globe_plt = st.sidebar.checkbox('View on Globe')
        if globe_plt:
            plot_globe_comparison('wnd', time_sel, int(level_sel))
        elif anomaly_plt:
            col1, col2 = st.columns(2)
            with col1:
                plot_wind_vectors(ds_u_subset[time_sel-1,:,:], ds_v_subset[time_sel-1,:,:],
                                  lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
            with col2:
                plot_wind_speed_anomaly(ds_u, ds_v, lat_min, lat_max, lon_min, lon_max, level_sel, time_sel, region_sel)
        else:
            plot_wind_vectors(ds_u_subset[time_sel-1,:,:], ds_v_subset[time_sel-1,:,:],
                              lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)

    elif plot_type == "Time Series":
        st.header("Wind Speed Time Series")
        st.write("Default location and pressure level is shown here. Please select your location of interest using latitude and longitude, and the pressure level")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                  ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))
        anomaly_plt = st.sidebar.checkbox('Show wind speed anomaly time series')
        if anomaly_plt:
            plot_wind_time_series_with_anomaly(speed_loc.sel(level=level_sel),
                                               direction_loc.sel(level=level_sel),
                                               lat_loc, lon_loc, level_sel)
        else:
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
        globe_plt = st.sidebar.checkbox('View on Globe')
        if globe_plt:
            plot_globe_comparison('air', time_sel)
        elif anomaly_plt:
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
        anomaly_plt = st.sidebar.checkbox('Show anomaly time series')
        if anomaly_plt:
            plot_time_series_with_anomaly(temp_loc, lat_loc, lon_loc)
        else:
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
        anomaly_plt = st.sidebar.checkbox('Check anomaly for the selected month')
        globe_plt = st.sidebar.checkbox('View on Globe')
        if globe_plt:
            plot_globe_comparison('precip', time_sel)
        elif anomaly_plt:
            col1, col2 = st.columns(2)
            with col1:
                plot_spatial2(ds_pr_subset, lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
            with col2:
                anomaly_plotting(ds_pr_subset, time_sel, region_sel)
        else:
            plot_spatial2(ds_pr_subset, lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
    else:
        st.header("Time series plot")
        st.write("Default location is shown here. Please select your location of interest using latitude and longitude")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        pr_loc = ds_pr['precip'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        anomaly_plt = st.sidebar.checkbox('Show anomaly time series')
        if anomaly_plt:
            plot_time_series_with_anomaly(pr_loc, lat_loc, lon_loc)
        else:
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
        anomaly_plt = st.sidebar.checkbox('Check anomaly for the selected month')
        globe_plt = st.sidebar.checkbox('View on Globe')
        if globe_plt:
            plot_globe_comparison('rhum', time_sel, int(level_sel))
        elif anomaly_plt:
            col1, col2 = st.columns(2)
            with col1:
                plot_spatial2(ds_rh_subset.sel(level=level_sel), lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
            with col2:
                anomaly_plotting(ds_rh_subset.sel(level=level_sel), time_sel, region_sel)
        else:
            plot_spatial2(ds_rh_subset.sel(level=level_sel), lat_min, lat_max, lon_min, lon_max,time_sel, region_sel)
    elif plot_type == 'Time Series':
        st.header("Time series plot")
        st.write("Default location and pressure level is shown here. Please select your location of interest using latitude and longitude, and the pressure level")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        rh_loc = ds_rh['rhum'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)
        anomaly_plt = st.sidebar.checkbox('Show anomaly time series')
        if anomaly_plt:
            plot_time_series_with_anomaly(rh_loc.sel(level=level_sel), lat_loc, lon_loc)
        else:
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
if anomaly_plt or globe_plt:
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"Climatology Data Source: {ds_temp.attrs['source']}")
    with col2:
        st.markdown(f"Current Data Source: Reanalysis-ERA5-Pressure-Levels-Monthly-Means")
else:
    st.markdown(f"Climatology Data Source: {ds_temp.attrs['source']}")

st.write("---")
st.markdown('''**Details of author**
Built by: K Narender Reddy
Email :email: : narender.kangari@ncas.ac.uk
Web Page :globe_with_meridians: : https://knreddy.online
Version 1.1: June, 2026''')