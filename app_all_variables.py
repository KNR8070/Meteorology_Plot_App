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
## select the spatial plot box
def select_box(lat,lon):
    # Latitude and Longitude Box Selection
    if ((lat_min > lat_max) or lon_min > lon_max):
        st.text("ERROR: Minimum value greater than Maximum")    
    return lat_min, lat_max, lon_min, lon_max
#%% [markdown] 
## Wind Rose Plot
def plot_wind_rose(speed_pwr, direction_pwr,lat_l,lon_l):
    # Ensure speed and direction are 1D arrays for the wind rose plot
    #speed_wr = np.array(speed_pwr)
    #direction_wr = np.array(direction_pwr)
    # if plt.fignum_exists(fig):
    #     fig=fig
    # else:
    fig_ws = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    ax.bar(direction_pwr, speed_pwr, normed=True, opening=0.8, edgecolor='white')
    ax.set_title('Latitude = '+str(lat_l)+' and Longitude = '+str(lon_l))
    ax.set_legend()
    st.pyplot()  

#%% [markdown] 
## Spatial Wind Vector Plot
def plot_wind_vectors(ds_u,ds_v, lat_min, lat_max, lon_min, lon_max, time_s):
    # Select data within specified lat/lon box   
    speed_mean = np.sqrt(ds_u**2 + ds_v**2)
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS)
    
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(ds_u.lon, ds_u.lat)
    speed_plot = ax.contourf(lons, lats, speed_mean, cmap='viridis', extend='both')
    fig.colorbar(speed_plot, ax=ax, label="Wind Speed (m/s)")

    ax.quiver(lons[::3], lats[::3], ds_u[::3], ds_v[::3])
    # if lon_max>180:
    #     ax.set_xticks(np.linspace(lon_min-180,lon_max-180,num=5,endpoint=True))
    # else:
    ax.set_xticks(np.linspace(lon_min,lon_max,num=5,endpoint=True))
    ax.set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # if climatology:
    #     ax.set_title('Plotted for time: Climatology (1991 to 2021)')
    # else:
    ax.set_title('Plotted for Month:'+calendar.month_name[time_s])
    st.pyplot(fig)
#%% [markdown] 
## Spatial Plot
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
          
    # ax[0,0].set_yticks(np.linspace(lat_min,lat_max,num=1,endpoint=True))
    # ax[1,0].set_yticks(np.linspace(lat_min,lat_max,num=1,endpoint=True))
    # ax[1,0].set_xticks([lon_min,lon_max])
    # ax[0,-1].set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    
    fig.colorbar(s_plot, ax=ax3, label="2m Temperature (degC)", shrink=0.75)
    st.pyplot(fig)
#%% [markdown]
##  Spatial plot test
def plot_spatial2(var_subset,lat_min, lat_max, lon_min, lon_max,time_s):
    # Select data within specified lat/lon box    
    # if lon_max>180:
    #     lon_max2 = lon_max
    #     lon_max2 = lon_max-180
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig, ax3 = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})#,
                           #sharex=True,sharey=True)#,figsize=(12,4))
    # if lon_max>180:
    #     ax3.set_extent([lon_min-(lon_max-180), 180, lat_min, lat_max], crs=ccrs.PlateCarree())
    # else:
    ax3.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS)
    
    # Plot contour fill for wind speed
    lons, lats = np.meshgrid(var_subset.lon, var_subset.lat)
    
    if var_subset.var_desc=='Air temperature':
        plot_data = np.squeeze(var_subset.isel(time=time_s))-273.15
        # if lon_max>180:
        #     # lons_new = np.hstack((temp_subset.lon[temp_subset.lon>180],lon
        #     plot_data2 = np.hstack([plot_data[:,plot_data.lon.values>180],plot_data[:,~(plot_data.lon.values>180)]])
        #     lons2 = np.hstack([lons[:,plot_data.lon.values>180],lons[:,~(plot_data.lon.values>180)]])
            
        #     s_plot = ax3.contourf(lons2, lats,plot_data2,cmap='viridis', extend='both')
        # else:
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
    # if lon_max>180:
    #     ax3.set_xticks(np.linspace(lon_min-(lon_max-180),180,num=5,endpoint=True))
    #     ax3.set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    # else:
    ax3.set_xticks(np.linspace(lon_min,lon_max,num=5,endpoint=True))
    ax3.set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
          
    # ax[0,0].set_yticks(np.linspace(lat_min,lat_max,num=1,endpoint=True))
    # ax[1,0].set_yticks(np.linspace(lat_min,lat_max,num=1,endpoint=True))
    # ax[1,0].set_xticks([lon_min,lon_max])
    # ax[0,-1].set_yticks(np.linspace(lat_min,lat_max,num=5,endpoint=True))
    
    
    st.pyplot(fig)
    
#%% [markdown] 
## Time Series Plot
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
    ax1.set_title('Latitude = '+str(speed_loc.lat.values)+' and Longitude = '+str(speed_loc.lon.values))
    # ax.legend('upper left')
    # ax2.legend()
    st.pyplot(fig)
#%% [markdown]
##  Time Series Plot
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
        ax1.plot(np.arange(1,13), var,color=color)
        ax1.set_xlabel("Month")
        ax1.set_ylabel("Relative humidity (%)",color=color, fontsize=14)
        ax1.set_yticks(np.linspace(np.floor(min(var.values)),np.floor(max(var.values)+1),num=5))
    
    ax1.set_xlabel("Month")
    ax1.set_xticks(np.arange(1,13))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                        'Sep','Oct','Nov', 'Dec'])
    # ax1.tick_params(axis='y', labelcolor=color)
    # ax.legend('upper left')
    # ax2.legend()
    st.pyplot(fig)
#%% [markdown] 
## Streamlit App
st.title("Met. Data Visualization App")
ds_temp = load_temp_data()
ds_u = load_uwind_data()
ds_v =load_vwind_data()
ds_pr =load_pr_data()
ds_rh = load_rh_data()
lon = ds_temp['lon']
lat = ds_temp['lat']

#%% [markdown]
## User Inputs  
var_type = st.sidebar.selectbox("Choose the variable", ("Temp_2m", "Wind", "Precipitation","Relative Humidity"))

if var_type == 'Wind':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Wind Rose", "Spatial Wind Vectors", "Time Series"))
    if plot_type == "Wind Rose":
        st.header("Wind Rose")
        # lat_loc = st.sidebar.selectbox("Select Latitude", lat.lat.values)
        # lon_loc = st.sidebar.selectbox("Select Longitude", lon.lon.values)
        lat_loc = st.sidebar.number_input("Enter Latitude", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lon_loc = st.sidebar.number_input("Enter Longitude", min_value=-180.00, max_value=180.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        # time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))

        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                          ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))

        # speed_at_loc = speed.sel(latitude=lat_loc, longitude=lon_loc,method='nearest')
        # direction_at_loc = direction.sel(latitude=lat_loc, longitude=lon_loc,method='nearest')

        plot_wind_rose(speed_loc.sel(level=level_sel), 
                       direction_loc.sel(level=level_sel),
                       speed_loc.lat.values,speed_loc.lon.values)

    elif plot_type == "Spatial Wind Vectors":
        st.header("Spatial Wind Vectors")        
        lat_min = st.sidebar.number_input("Enter Lat. min.", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lat_max = st.sidebar.number_input("Enter Lat. max.", min_value=-90.00, max_value=90.00, 
                                          value=40.00,step=0.01, format='%2.2f',placeholder="Must be greater than Lat min.")
        lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=0.00, max_value=360.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=0.00, max_value=360.00, 
                                       value=100.00,step=0.01, format='%3.2f',placeholder="Must be greater than Lon min.")
        
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        
        ds_u_subset = ds_u['uwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max),level=level_sel)
        ds_v_subset = ds_v['vwnd'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max),level=level_sel)
        plot_wind_vectors(ds_u_subset[time_sel-1,:,:], ds_u_subset[time_sel-1,:,:], lat_min, lat_max, lon_min, lon_max,time_sel)

    elif plot_type == "Time Series":
        st.header("Wind Speed Time Series")
        lat_loc = st.sidebar.number_input("Enter Latitude", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lon_loc = st.sidebar.number_input("Enter Longitude", min_value=0.00, max_value=360.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_u.level.values)
        speed_loc, direction_loc = calculate_wind(ds_u['uwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'),
                                  ds_v['vwnd'].sel(lat=lat_loc,lon=lon_loc,method='nearest'))
        plot_time_series(speed_loc.sel(level=level_sel),direction_loc.sel(level=level_sel))
########################### 2m Temperature
elif var_type == 'Temp_2m':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")        
        lat_min = st.sidebar.number_input("Enter Lat. min.", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lat_max = st.sidebar.number_input("Enter Lat. max.", min_value=-90.00, max_value=90.00, 
                                          value=40.00,step=0.01, format='%2.2f',placeholder="Must be greater than Lat min.")
        lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=-180.00, max_value=180.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=-180.00, max_value=180.00, 
                                       value=100.00,step=0.01, format='%3.2f',placeholder="Must be greater than Lon min.")
        ds_temp_subset = ds_temp['air'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        # plot_spatial(ds_temp_subset, lat_min, lat_max, lon_min, lon_max)
        plot_spatial2(ds_temp_subset.sel(level=2), lat_min, lat_max, lon_min, lon_max,time_sel)
        # plot_spatial2(ds_temp_subset)
        # ds_temp_subset.plot(col='time',col_wrap=6)
    else:
        st.header("Time series plot")
        lat_loc = st.sidebar.number_input("Enter Latitude", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lon_loc = st.sidebar.number_input("Enter Longitude", min_value=0.00, max_value=360.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        temp_loc = ds_temp['air'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        plot_time_series2(temp_loc.isel(level=0))
##################################### Precipitation        
elif var_type == 'Precipitation':
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")        
        lat_min = st.sidebar.number_input("Enter Lat. min.", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lat_max = st.sidebar.number_input("Enter Lat. max.", min_value=-90.00, max_value=90.00, 
                                          value=40.00,step=0.01, format='%2.2f',placeholder="Must be greater than Lat min.")
        lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=-180.00, max_value=180.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=-180.00, max_value=180.00, 
                                       value=100.00,step=0.01, format='%3.2f',placeholder="Must be greater than Lon min.")
        ds_pr_subset = ds_pr['precip'].sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        # plot_spatial(ds_temp_subset, lat_min, lat_max, lon_min, lon_max)
        plot_spatial2(ds_pr_subset, lat_min, lat_max, lon_min, lon_max,time_sel)
        # plot_spatial2(ds_temp_subset)
        # ds_temp_subset.plot(col='time',col_wrap=6)
    else:
        st.header("Time series plot")
        lat_loc = st.sidebar.number_input("Enter Latitude", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lon_loc = st.sidebar.number_input("Enter Longitude", min_value=0.00, max_value=360.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        pr_loc = ds_pr['precip'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        plot_time_series2(pr_loc)
######################################## Relative Humidity
else:
    plot_type = st.sidebar.selectbox("Choose Plot Type", ("Spatial plot", "Time Series"))
    if plot_type == 'Spatial plot':
        st.header("Spatial plot")        
        lat_min = st.sidebar.number_input("Enter Lat. min.", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lat_max = st.sidebar.number_input("Enter Lat. max.", min_value=-90.00, max_value=90.00, 
                                          value=40.00,step=0.01, format='%2.2f',placeholder="Must be greater than Lat min.")
        lon_min = st.sidebar.number_input("Enter Lon. min.", min_value=-180.00, max_value=180.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        lon_max = st.sidebar.number_input("Enter Lon. max.", min_value=-180.00, max_value=180.00, 
                                       value=100.00,step=0.01, format='%3.2f',placeholder="Must be greater than Lon min.")
        ds_rh_subset = ds_rh['rhum'].sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
        time_sel = st.sidebar.selectbox("Select Month", np.arange(1,13))
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)
        # plot_spatial(ds_temp_subset, lat_min, lat_max, lon_min, lon_max)
        plot_spatial2(ds_rh_subset.sel(level=level_sel), lat_min, lat_max, lon_min, lon_max,time_sel)
        # plot_spatial2(ds_temp_subset)
        # ds_temp_subset.plot(col='time',col_wrap=6)
    else:
        st.header("Time series plot")
        lat_loc = st.sidebar.number_input("Enter Latitude", min_value=-90.00, max_value=90.00, 
                                          value=0.00,step=0.01, format='%2.2f')
        lon_loc = st.sidebar.number_input("Enter Longitude", min_value=0.00, max_value=360.00, 
                                       value=0.00,step=0.01, format='%3.2f')
        rh_loc = ds_rh['rhum'].sel(lat=lat_loc,lon=lon_loc,method='nearest')
        level_sel = st.sidebar.selectbox("Select Level (hPa)", ds_rh.level.values)
        plot_time_series2(rh_loc.sel(level=level_sel))
    
    
    
    