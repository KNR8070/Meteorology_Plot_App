#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 20:54:39 2024

@author: knreddy
"""
#%% [markdown]
### Meteorological Visualisation App.
#%% [markdown]
## load modules
import streamlit as st
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import calendar
#%% [markdown] 
## Load wind data from NetCDF
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
    ax.set_legend()
    ax.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax.transAxes)
    st.pyplot()  
#%% [markdown] 
## Spatial Wind Vector Plot
def plot_wind_vectors(ds_u,ds_v, lat_min, lat_max, lon_min, lon_max, time_s):
    # Select data within specified lat/lon box   
    speed_mean = np.sqrt(ds_u**2 + ds_v**2)
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
    fig, ax = plt.subplots(figsize=(x_size,y_size),
                           subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(ds_u.lon, ds_u.lat)
    speed_plot = ax.contourf(lons, lats, speed_mean, cmap='viridis', extend='both')
    if x_ratio>y_ratio:
        fig.colorbar(speed_plot, ax=ax, label="Wind Speed (m/s)",shrink=0.4)
    else:
        fig.colorbar(speed_plot, ax=ax, label="Wind Speed (m/s)",shrink=0.7)
    if (lat_max-lat_min)>60 and (lon_max-lon_min)>60:
        alt_num = 2
    else:
        alt_num = 1
    Q = ax.quiver(lons[::alt_num,::alt_num], 
              lats[::alt_num,::alt_num], 
              ds_u[::alt_num,::alt_num], 
              ds_v[::alt_num,::alt_num],pivot='mid', units='inches', alpha=0.6)
    if ds_u.level.values<800:
        qk = ax.quiverkey(Q, 0.8, 0.9, 15, r'$15 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    else:
        qk = ax.quiverkey(Q, 0.8, 0.9, 5, r'$5 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    ax.set_xticks(np.linspace(lon_min,lon_max,num=5,endpoint=True))
    ax.set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Month:'+calendar.month_name[time_s][:3]+'  Level:'+str(ds_u.level.values)+' '+ds_u.level.GRIB_name)
    ax.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax.transAxes)
    st.pyplot(fig)
#%% [markdown] 
# Function to plot wind vector plot
def plot_spatial(temp_subset, lat_min, lat_max, lon_min, lon_max):
    # Select data within specified lat/lon box
    fig, ax3 = plt.subplots(nrows=2,ncols=6,subplot_kw={'projection': ccrs.PlateCarree()},
                           sharex=True,sharey=True,figsize=(12,4))
    for i_row in np.arange(2):
        for i_col in np.arange(6):
            ax3[i_row,i_col].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax3[i_row,i_col].add_feature(cfeature.COASTLINE)
            # ax.add_feature(cfeature.BORDERS)
            
            # Plot contour fill for wind speed
            lons, lats = np.meshgrid(temp_subset.lon, temp_subset.lat)
            s_plot = ax3[i_row,i_col].contourf(lons, lats, 
                                              np.squeeze(temp_subset.isel(time=int(i_row*6+i_col)))-273.15, 
                                              cmap='viridis', extend='both')
            ax3[i_row,i_col].set_title(calendar.month_name[int(i_row*6+i_col)+1][:3])
            
            if i_row==1:
                ax3[i_row,i_col].set_xlabel('Longitude')
                ax3[i_row,i_col].set_xticks(np.linspace(lon_min,lon_max,num=2,endpoint=True))
                ax3[i_row,i_col].set_xticklabels([lon_min,lon_max],fontsize=6)
                if i_col==0:
                    ax3[i_row,i_col].set_ylabel('Latitude')

            else: #0th Row
                if i_col==0:
                    ax3[i_row,i_col].set_ylabel('Latitude')
                    ax3[i_row,i_col].set_yticks(np.linspace(lat_min,lat_max,num=2,endpoint=True))
                else:                    
                    ax3[i_row,i_col].set_yticks(np.linspace(lat_min,lat_max,num=2,endpoint=True))
    
    fig.colorbar(s_plot, ax=ax3, label="2m Temperature (degC)", shrink=0.75)
    ax3.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax3.transAxes)
    st.pyplot(fig)
#%% [markdown]
#Function to plot spatial variation in variables
def plot_spatial2(var_subset,lat_min, lat_max, lon_min, lon_max,time_s):
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig, ax3 = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})#,

    ax3.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.COASTLINE)
    
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(var_subset.lon, var_subset.lat)
    
    if var_subset.var_desc=='Air temperature':
        plot_data = np.squeeze(var_subset.isel(time=time_s))-273.15
        s_plot = ax3.contourf(lons,lats,plot_data,cmap='viridis', extend='both')
        fig.colorbar(s_plot, ax=ax3, label="2m Temperature (degC)", shrink=0.75)
        
    elif var_subset.var_desc=='Precipitation':
        plot_data = np.squeeze(var_subset.isel(time=time_s))
        s_plot = ax3.contourf(lons,lats,plot_data,cmap='viridis', extend='both')
        fig.colorbar(s_plot, ax=ax3, label="Mean Precipitation (mm/day)", shrink=0.75)
    else:
        plot_data = np.squeeze(var_subset.isel(time=time_s))
        s_plot = ax3.contourf(lons,lats,plot_data,cmap='viridis', extend='both')
        fig.colorbar(s_plot, ax=ax3, label="Relative humidity (%)", shrink=0.75)
        
    ax3.set_title(calendar.month_name[time_s][:3])
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax3.set_xticks(np.linspace(lon_min,lon_max,num=5,endpoint=True))
    ax3.set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    ax3.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax3.transAxes)
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
    ax1.text(0.7,-0.1,'Data Source: '+ds_temp.attrs['source'],fontsize=6,transform=ax1.transAxes)
    st.pyplot(fig)
#%% [markdown]
# Function to plot time series of variables
def plot_time_series2(var,lat_loc,lon_loc):
    # plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    fig, ax1 = plt.subplots(figsize=(12,6))
    #subfig, subax = fig.subfigures()
    s_lat_min = var.lat.values-15
    s_lat_max = var.lat.values+15
    s_lon_min = var.lon.values-15
    s_lon_max = var.lon.values+15
    subpos = [0.6,0.6,0.7,0.7]
    subax = add_subplot_axes(ax1,subpos)

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
        var_loc = var['rhum'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        ax1.plot(np.arange(1,13), var_loc.values,color=color)
        ax1.set_xlabel("Month")
        ax1.set_ylabel("Relative humidity (%)",color=color, fontsize=14)
        ax1.set_yticks(np.linspace(np.floor(min(var_loc.values)),
                                   np.floor(max(var_loc.values)+1),num=5))
        subax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        subax.add_feature(cfeature.COASTLINE)
        subax.scatter(lon_loc,lat_loc,marker='*')
    
    ax1.set_xlabel("Month")
    ax1.set_xticks(np.arange(1,13))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                        'Sep','Oct','Nov', 'Dec'])
    ax1.text(0.7,-0.2,'Data Source: '+ds_temp.attrs['source'],fontsize=4,transform=ax1.transAxes)
    st.pyplot(fig)

#%% [markdown] #URL: https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
# subax
def add_subplot_axes(ax,rect,facecolor='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    #subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
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
                                    value=-85.00,step=0.01, format='%2.2f')
    lat_max = st.sidebar.number_input("Enter Lat. max.", float(str(lat.values.min())), 
                                    max_value=float(str(lat.values.max())), 
                                    value=85.00,step=0.01, format='%2.2f',
                                    placeholder="Must be greater than Lat min.")
    lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=float(str(lon.values.min())), 
                                    max_value=float(str(lon.values.max())), 
                                    value=-178.00,step=0.01, format='%3.2f')
    lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=float(str(lon.values.min())), 
                                    max_value=float(str(lon.values.max())), 
                                    value=178.00,step=0.01, format='%3.2f',
                                    placeholder="Must be greater than Lon min.")
    return lat_min, lat_max, lon_min, lon_max
#%% [markdown]
# function to find lat, lon location for plotting
def user_input_loc(lat,lon):
    lat_loc = st.sidebar.number_input("Enter Latitude", min_value=float(str(lat.values.min())), 
                                          max_value=float(str(lat.values.max())), 
                                          value=0.00,step=0.01, format='%2.2f')
    lon_loc = st.sidebar.number_input("Enter Longitude", min_value=float(str(lon.values.min())), 
                                          max_value=float(str(lon.values.max())), 
                                          value=0.00,step=0.01, format='%3.2f')
    return lat_loc, lon_loc
#%% [markdown] 
## Streamlit App
st.title("Met. Data Visualization App")
ds_temp = convert_180_180(load_temp_data())
ds_u = convert_180_180(load_uwind_data())
ds_v = convert_180_180(load_vwind_data())
ds_pr = convert_180_180(load_pr_data())
ds_rh = convert_180_180(load_rh_data())
lon = ds_temp['lon']
lat = ds_temp['lat']
#%% [markdown]
# User Inputs  
var_type = st.sidebar.selectbox("Choose the variable", ("Temp_2m", "Wind", "Precipitation","Relative Humidity"))

if var_type == 'Wind':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Wind Rose", "Spatial Wind Vectors", "Time Series"))
    if plot_type == "Wind Rose":
        st.header("Wind Rose")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)

        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                          ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))

        plot_wind_rose(speed_loc.sel(level=level_sel).values, 
                       direction_loc.sel(level=level_sel).values,
                       speed_loc.lat.values,speed_loc.lon.values)

    elif plot_type == "Spatial Wind Vectors":
        st.header("Spatial Wind Vectors")        
        [lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        
        ds_u_subset = ds_u['uwnd'].sel(lat=slice(lat_max, lat_min), 
                                       lon=slice(lon_min, lon_max),
                                       level=level_sel)
        ds_v_subset = ds_v['vwnd'].sel(lat=slice(lat_max, lat_min), 
                                       lon=slice(lon_min, lon_max),
                                       level=level_sel)
        plot_wind_vectors(ds_u_subset[time_sel-1,:,:], ds_u_subset[time_sel-1,:,:], 
                          lat_min, lat_max, lon_min, lon_max,time_sel)

    elif plot_type == "Time Series":
        st.header("Wind Speed Time Series")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                  ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))
        plot_time_series(speed_loc.sel(level=level_sel),direction_loc.sel(level=level_sel))
########################### 2m Temperature
elif var_type == 'Temp_2m':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")        
        [lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        ds_temp_subset = ds_temp['air'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))

        plot_spatial2(ds_temp_subset.sel(level=2), lat_min, lat_max, lon_min, lon_max,time_sel)

    else:
        st.header("Time series plot")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        temp_loc = ds_temp['air'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        plot_time_series2(temp_loc.isel(level=0))
##################################### Precipitation        
elif var_type == 'Precipitation':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")
        [lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)        
        ds_pr_subset = ds_pr['precip'].sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        plot_spatial2(ds_pr_subset, lat_min, lat_max, lon_min, lon_max,time_sel)
    else:
        st.header("Time series plot")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        pr_loc = ds_pr['precip'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        plot_time_series2(pr_loc)
######################################## Relative Humidity
else:
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")        
        [lat_min, lat_max, lon_min, lon_max] = user_input_box(lat,lon)
        ds_rh_subset = ds_rh['rhum'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)

        plot_spatial2(ds_rh_subset.sel(level=level_sel), lat_min, lat_max, lon_min, lon_max,time_sel)
    else:
        st.header("Time series plot")
        lat_loc, lon_loc = user_input_loc(lat,lon)
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)
        plot_time_series2(ds_rh.sel(level=level_sel),lat_loc,lon_loc)  
# %%
