#!/usr/bin/env python3
"""
Download ERA5 monthly mean data for the current and previous calendar year,
rename variables/coordinates to match NCEP naming used by the Streamlit app,
and write to data/era5_<start>_<end>_<var>.nc.

Variables downloaded:
  air  — 2 m temperature        (single-level)
  wnd  — u + v wind components  (pressure levels)
  rhum — relative humidity       (pressure levels)

CDS authentication: either set CDSAPI_URL / CDSAPI_KEY environment variables
(GitHub Actions workflow does this) or supply a ~/.cdsapirc file.
"""
import glob
import os
import sys
from datetime import date

import cdsapi
import numpy as np
import xarray as xr

# 1° grid keeps files well under GitHub's 100 MB limit while still being
# finer than the 2.5° NCEP climatology used for anomaly differencing.
GRID = '1.0/1.0'
PRESSURE_LEVELS = ['1000', '925', '850', '700', '500', '300', '200', '100']
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


def _year_range():
    today = date.today()
    end_year = today.year
    start_year = end_year - 1
    return start_year, end_year


def _all_months(start_year, end_year):
    years  = [str(y) for y in range(start_year, end_year + 1)]
    months = [f'{m:02d}' for m in range(1, 13)]
    return years, months


def _normalise_coords(ds):
    """Rename valid_time → time and latitude/longitude → lat/lon if present."""
    rename = {}
    if 'valid_time' in ds.dims and 'time' not in ds.dims:
        rename['valid_time'] = 'time'
    if 'latitude' in ds.dims:
        rename['latitude'] = 'lat'
    if 'longitude' in ds.dims:
        rename['longitude'] = 'lon'
    if 'pressure_level' in ds.dims:
        rename['pressure_level'] = 'level'
    return ds.rename(rename) if rename else ds


def _write(ds, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    ds.to_netcdf(path)
    print(f"Written: {path}")


def _remove_stale(var_name, keep_path):
    """Delete older era5_*_<var>.nc files that are no longer current."""
    for f in glob.glob(os.path.join(DATA_DIR, f'era5_*_{var_name}.nc')):
        if os.path.abspath(f) != os.path.abspath(keep_path):
            os.remove(f)
            print(f"Removed stale file: {f}")


def fetch_temperature(client, start_year, end_year):
    years, months = _all_months(start_year, end_year)
    tmp = '/tmp/era5_tmp_air.nc'
    client.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable':     '2m_temperature',
            'year':         years,
            'month':        months,
            'time':         '00:00',
            'grid':         GRID,
            'format':       'netcdf',
        },
        tmp,
    )

    ds = _normalise_coords(xr.open_dataset(tmp))
    if 't2m' in ds:
        ds = ds.rename({'t2m': 'air'})
    ds['air'].attrs['long_name'] = '2m Air Temperature'
    ds['air'].attrs['units']     = 'K'

    out = os.path.join(DATA_DIR, f'era5_{start_year}_{end_year}_air.nc')
    _write(ds, out)
    _remove_stale('air', out)
    return out


def fetch_wind(client, start_year, end_year):
    years, months = _all_months(start_year, end_year)
    tmp = '/tmp/era5_tmp_wnd.nc'
    client.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type':   'monthly_averaged_reanalysis',
            'variable':       ['u_component_of_wind', 'v_component_of_wind'],
            'pressure_level': PRESSURE_LEVELS,
            'year':           years,
            'month':          months,
            'time':           '00:00',
            'grid':           GRID,
            'format':         'netcdf',
        },
        tmp,
    )

    ds = _normalise_coords(xr.open_dataset(tmp))
    # ERA5 pressure-level wind variable names: 'u' and 'v'
    rename_map = {}
    if 'u' in ds:
        rename_map['u'] = 'uwnd'
    if 'v' in ds:
        rename_map['v'] = 'vwnd'
    if rename_map:
        ds = ds.rename(rename_map)
    ds['uwnd'].attrs['long_name'] = 'U-component of wind'
    ds['vwnd'].attrs['long_name'] = 'V-component of wind'

    out = os.path.join(DATA_DIR, f'era5_{start_year}_{end_year}_wnd.nc')
    _write(ds, out)
    _remove_stale('wnd', out)
    return out


def fetch_precipitation(client, start_year, end_year):
    years, months = _all_months(start_year, end_year)
    tmp = '/tmp/era5_tmp_precip.nc'
    client.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable':     'total_precipitation',
            'year':         years,
            'month':        months,
            'time':         '00:00',
            'grid':         GRID,
            'format':       'netcdf',
        },
        tmp,
    )

    ds = _normalise_coords(xr.open_dataset(tmp))
    # ERA5 monthly mean tp is in m/day; convert to mm/day to match NCEP
    if 'tp' in ds:
        ds = ds.rename({'tp': 'precip'})
    ds['precip'] = ds['precip'] * 1000.0
    ds['precip'].attrs['long_name'] = 'Total Precipitation'
    ds['precip'].attrs['units']     = 'mm/day'

    out = os.path.join(DATA_DIR, f'era5_{start_year}_{end_year}_precip.nc')
    _write(ds, out)
    _remove_stale('precip', out)
    return out


def fetch_rhum(client, start_year, end_year):
    years, months = _all_months(start_year, end_year)
    tmp = '/tmp/era5_tmp_rhum.nc'
    client.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type':   'monthly_averaged_reanalysis',
            'variable':       'relative_humidity',
            'pressure_level': PRESSURE_LEVELS,
            'year':           years,
            'month':          months,
            'time':           '00:00',
            'grid':           GRID,
            'format':         'netcdf',
        },
        tmp,
    )

    ds = _normalise_coords(xr.open_dataset(tmp))
    # ERA5 relative humidity variable name: 'r'
    if 'r' in ds:
        ds = ds.rename({'r': 'rhum'})
    ds['rhum'].attrs['long_name'] = 'Relative Humidity'
    ds['rhum'].attrs['units']     = '%'

    out = os.path.join(DATA_DIR, f'era5_{start_year}_{end_year}_rhum.nc')
    _write(ds, out)
    _remove_stale('rhum', out)
    return out


def main():
    start_year, end_year = _year_range()
    print(f"Fetching ERA5 data for {start_year}–{end_year} ...")

    # CDS client picks up credentials from env vars CDSAPI_URL / CDSAPI_KEY
    # or from ~/.cdsapirc
    client = cdsapi.Client()

    fetch_temperature(client, start_year, end_year)
    fetch_wind(client, start_year, end_year)
    fetch_rhum(client, start_year, end_year)
    fetch_precipitation(client, start_year, end_year)

    print("All ERA5 variables downloaded and saved.")


if __name__ == '__main__':
    main()
