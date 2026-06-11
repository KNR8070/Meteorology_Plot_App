#!/usr/bin/env python3
"""
3D Interactive Globe View for Met variables.
Uses Plotly go.Surface to render gridded data on a sphere — fully interactive
(rotate, zoom, pan) inside the browser via st.plotly_chart().

Run with:
    streamlit run app_globe.py
"""
import calendar

import cartopy.feature as cfeature
import numpy as np
import plotly.graph_objects as go
import streamlit as st
import xarray as xr

from app_all_variables import (
    _find_era5_file,
    convert_180_180,
    load_pr_data,
    load_rh_data,
    load_temp_data,
    load_uwind_data,
    load_vwind_data,
    select_current_year_data,
)

# ── page config ───────────────────────────────────────────────────────────────
st.set_page_config(layout="wide", page_title="Globe View")
st.title("Met. Data – 3D Globe View")

# ── constants ─────────────────────────────────────────────────────────────────
VAR_OPTS = {
    "Air Temperature (2m)": "air",
    "Wind Speed":           "wnd",
    "Precipitation":        "precip",
    "Relative Humidity":    "rhum",
}
CMAPS = {
    "air":    "RdBu_r",
    "wnd":    "Viridis",
    "precip": "BrBG",
    "rhum":   "RdYlGn",
}
ANOM_CMAPS = {
    "air":    "RdBu",
    "wnd":    "RdBu",
    "precip": "BrBG",
    "rhum":   "RdYlGn",
}
UNITS = {"air": "°C", "wnd": "m/s", "precip": "mm/day", "rhum": "%"}
CLIM_RANGE = {"air": (-30, 40), "wnd": (0, 20), "precip": (0, 15), "rhum": (0, 100)}
ANOM_RANGE = {"air": (-10, 10), "wnd": (-5, 5),  "precip": (-8, 8),  "rhum": (-20, 20)}


# ── sidebar ───────────────────────────────────────────────────────────────────
st.sidebar.title("Globe Settings")
var_label = st.sidebar.selectbox("Variable", list(VAR_OPTS))
var_key   = VAR_OPTS[var_label]
month     = st.sidebar.selectbox(
    "Month", list(range(1, 13)), format_func=lambda m: calendar.month_name[m]
)
view_mode = st.sidebar.radio(
    "View",
    ["Climatology", "ERA5 Current", "Anomaly (ERA5 − Climatology)",
     "Side-by-side comparison"],
)

level_sel = None
if var_key in ("wnd", "rhum"):
    src_ds    = load_uwind_data() if var_key == "wnd" else load_rh_data()
    level_sel = int(st.sidebar.selectbox("Pressure level (hPa)",
                                          src_ds["level"].values))

coastline_color = st.sidebar.color_picker("Coastline colour", "#ffffff")


# ── coastlines (cached) ───────────────────────────────────────────────────────
@st.cache_data
def _coastline_xyz(r: float = 1.005):
    """Return list of (x,y,z) arrays for coastline segments on a sphere of radius r."""
    segments = []
    for geom in cfeature.COASTLINE.geometries():
        lines = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        for line in lines:
            coords = np.array(line.coords)
            if len(coords) < 2:
                continue
            lons_c = np.radians(coords[:, 0])
            lats_c = np.radians(coords[:, 1])
            theta  = np.pi / 2 - lats_c
            x = r * np.sin(theta) * np.cos(lons_c)
            y = r * np.sin(theta) * np.sin(lons_c)
            z = r * np.cos(theta)
            # Insert NaN between segments so one Scatter3d trace covers all coasts
            segments.append(np.column_stack([x, y, z]))
    # Stack all segments separated by NaN rows
    nan_row   = np.full((1, 3), np.nan)
    combined  = np.vstack([np.vstack([seg, nan_row]) for seg in segments])
    return combined[:, 0], combined[:, 1], combined[:, 2]


# ── sphere mesh from DataArray ────────────────────────────────────────────────
def _da_to_sphere(da: xr.DataArray):
    """Map a 2-D lat/lon DataArray to sphere coordinates, wrapping lon to close the seam."""
    da        = da.sortby("lat").sortby("lon")
    lon_vals  = da.lon.values.astype(float)
    lat_vals  = da.lat.values.astype(float)
    data_vals = da.values.astype(float)

    # Wrap: append first lon column at end so the mesh closes at ±180°
    lon_step = float(lon_vals[1] - lon_vals[0]) if len(lon_vals) > 1 else 1.0
    lon_ext  = np.append(lon_vals, lon_vals[-1] + lon_step)
    surf_ext = np.hstack([data_vals, data_vals[:, :1]])

    lon2d, lat2d = np.meshgrid(lon_ext, lat_vals)
    theta = np.radians(90.0 - lat2d)
    phi   = np.radians(lon2d)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    lon2d_h, _ = np.meshgrid(np.append(lon_vals, lon_vals[0]), lat_vals)
    return x, y, z, surf_ext, lon2d_h, lat2d


# ── build a single Plotly figure ──────────────────────────────────────────────
def build_globe(da: xr.DataArray, title: str, cmap: str,
                vmin: float, vmax: float, height: int = 650) -> go.Figure:
    x, y, z, surf, lon2d_h, lat2d = _da_to_sphere(da)
    cx, cy, cz = _coastline_xyz()

    surface = go.Surface(
        x=x, y=y, z=z,
        surfacecolor=surf,
        colorscale=cmap,
        cmin=vmin, cmax=vmax,
        colorbar=dict(
            title=dict(text=f"{var_label}<br>({UNITS[var_key]})", side="right"),
            thickness=14, len=0.6,
        ),
        showscale=True,
        lighting=dict(ambient=0.8, diffuse=0.5, specular=0.1),
        hovertemplate="lon: %{customdata[0]:.1f}°<br>lat: %{customdata[1]:.1f}°<br>"
                      f"value: %{{customdata[2]:.2f}} {UNITS[var_key]}<extra></extra>",
        customdata=np.dstack([lon2d_h, lat2d, surf]),
    )
    coasts = go.Scatter3d(
        x=cx, y=cy, z=cz,
        mode="lines",
        line=dict(color=coastline_color, width=1),
        hoverinfo="none",
        showlegend=False,
    )

    fig = go.Figure(data=[surface, coasts])

    # Animation frames: rotate camera 360° around vertical axis
    n_frames = 72
    r_eye    = np.sqrt(1.4**2 + 0.5**2)
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
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=15)),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor="black",
            aspectmode="cube",
            camera=dict(eye=dict(
                x=float(r_eye * np.cos(0)),
                y=float(r_eye * np.sin(0)),
                z=z_eye,
            )),
        ),
        margin=dict(l=0, r=0, t=40, b=40),
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


# ── data helpers ──────────────────────────────────────────────────────────────
@st.cache_data
def get_clim(var_key: str, month_idx: int, level):
    if var_key == "air":
        ds = convert_180_180(load_temp_data()).sel(lat=slice(85, -85), lon=slice(-180, 180))
        da = ds["air"].isel(time=month_idx) - 273.15
        da.name = "air"
    elif var_key == "wnd":
        ds_u = convert_180_180(load_uwind_data()).sel(lat=slice(85, -85), lon=slice(-180, 180))
        ds_v = convert_180_180(load_vwind_data()).sel(lat=slice(85, -85), lon=slice(-180, 180))
        u  = ds_u["uwnd"].sel(level=level).isel(time=month_idx)
        v  = ds_v["vwnd"].sel(level=level).isel(time=month_idx)
        da = np.sqrt(u ** 2 + v ** 2)
        da.name = "wnd"
    elif var_key == "precip":
        ds = convert_180_180(load_pr_data()).sel(lat=slice(-85, 85), lon=slice(-180, 180))
        da = ds["precip"].sortby("lat").isel(time=month_idx)
    elif var_key == "rhum":
        ds = convert_180_180(load_rh_data()).sel(lat=slice(85, -85), lon=slice(-180, 180))
        da = ds["rhum"].sel(level=level).isel(time=month_idx)
    return da.squeeze()


def get_era5(var_key: str, month: int, level):
    path = _find_era5_file(var_key)
    try:
        ds = xr.open_dataset(path)
    except (FileNotFoundError, OSError):
        return None, None
    ds = convert_180_180(ds).sel(lat=slice(85, -85), lon=slice(-180, 180))
    if var_key == "wnd":
        u  = ds["uwnd"].sel(level=level) if level else ds["uwnd"].isel(level=0)
        v  = ds["vwnd"].sel(level=level) if level else ds["vwnd"].isel(level=0)
        da = np.sqrt(u ** 2 + v ** 2)
        da.name = "wnd"
    else:
        da = ds[var_key]
        if "level" in da.dims:
            da = da.sel(level=level) if level else da.isel(level=0)
        if var_key == "air":
            da = da - 273.15
    da_pt, year = select_current_year_data(da, month)
    if da_pt is None:
        return None, None
    return da_pt.squeeze(), year


def get_anomaly(var_key, month, level):
    da_era5, year = get_era5(var_key, month, level)
    if da_era5 is None:
        return None, None, None
    da_clim = get_clim(var_key, month - 1, level)
    da_clim_i = da_clim.sortby("lat").interp(
        lat=np.sort(da_era5.lat.values),
        lon=np.sort(da_era5.lon.values),
        method="linear",
    )
    anom = da_era5.sortby(["lat", "lon"]) - da_clim_i
    anom.name = var_key + "_anomaly"
    return anom, year, da_clim


# ── main ──────────────────────────────────────────────────────────────────────
month_name = calendar.month_name[month]

with st.spinner("Building globe ..."):

    if view_mode == "Climatology":
        da = get_clim(var_key, month - 1, level_sel)
        vmin, vmax = CLIM_RANGE[var_key]
        fig = build_globe(da, f"{var_label} Climatology — {month_name}",
                          CMAPS[var_key], vmin, vmax)
        st.plotly_chart(fig, use_container_width=True)
        st.caption("NCEP 1991–2020 climatology")

    elif view_mode == "ERA5 Current":
        da, year = get_era5(var_key, month, level_sel)
        if da is None:
            st.error(f"No ERA5 file found for '{var_key}'. Run scripts/fetch_era5.py first.")
        else:
            vmin, vmax = CLIM_RANGE[var_key]
            fig = build_globe(da, f"{var_label} ERA5 — {month_name} {year}",
                              CMAPS[var_key], vmin, vmax)
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"ERA5 monthly mean · {month_name} {year}")

    elif view_mode == "Anomaly (ERA5 − Climatology)":
        anom, year, _ = get_anomaly(var_key, month, level_sel)
        if anom is None:
            st.error(f"No ERA5 file found for '{var_key}'. Run scripts/fetch_era5.py first.")
        else:
            vmin, vmax = ANOM_RANGE[var_key]
            fig = build_globe(anom,
                              f"{var_label} Anomaly — {month_name} {year}",
                              ANOM_CMAPS[var_key], vmin, vmax)
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"ERA5 {year} − NCEP 1991–2020 climatology · {month_name}")

    else:  # Side-by-side comparison
        anom, year, da_clim = get_anomaly(var_key, month, level_sel)
        if anom is None:
            st.error(f"No ERA5 file found for '{var_key}'. Run scripts/fetch_era5.py first.")
        else:
            da_era5, _ = get_era5(var_key, month, level_sel)
            vmin_c, vmax_c = CLIM_RANGE[var_key]
            vmin_a, vmax_a = ANOM_RANGE[var_key]

            col1, col2, col3 = st.columns(3)
            with col1:
                fig1 = build_globe(da_clim,
                                   f"Climatology — {month_name}",
                                   CMAPS[var_key], vmin_c, vmax_c, height=450)
                st.plotly_chart(fig1, use_container_width=True)
                st.caption("NCEP 1991–2020 climatology")
            with col2:
                fig2 = build_globe(da_era5,
                                   f"ERA5 — {month_name} {year}",
                                   CMAPS[var_key], vmin_c, vmax_c, height=450)
                st.plotly_chart(fig2, use_container_width=True)
                st.caption(f"ERA5 · {month_name} {year}")
            with col3:
                fig3 = build_globe(anom,
                                   f"Anomaly — {month_name} {year}",
                                   ANOM_CMAPS[var_key], vmin_a, vmax_a, height=450)
                st.plotly_chart(fig3, use_container_width=True)
                st.caption(f"ERA5 {year} − Climatology")

st.write("---")
st.markdown(
    "**Built by:** K Narender Reddy &nbsp;|&nbsp;"
    "**Email:** knreddyiitd@gmail.com &nbsp;|&nbsp;"
    "**Web:** https://knreddy.online"
)
