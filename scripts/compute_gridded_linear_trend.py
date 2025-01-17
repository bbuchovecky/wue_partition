import time
import numpy as np
import xarray as xr
import xcdat as xc

def format_ds_coords(ds):
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

def calculate_annual_timeseries(da):
    """
    Calculates the annual timeseries, weighted by the number of days in each month
    """
    nyears = len(da.groupby('time.year'))
    month_length = da.time.dt.days_in_month

    weights = month_length.groupby('time.year') / month_length.astype(float).groupby('time.year').sum()        
    np.testing.assert_allclose(weights.groupby('time.year').sum().values, np.ones(nyears)) 

    return (da * weights).groupby('time.year').sum(dim='time')

global_metadata = {
    # 'description': 'Gridded to (time,lat,lon,pft) from (time,pft). All global attributes below are copied from the original history file.',
    'title': 'CLM History file information',
    'comment': 'NOTE: None of the variables are weighted by land fraction!',
    'Conventions': 'CF-1.0',
    'history': 'created on 04/27/21 13:27:32',
    'source': 'Community Terrestrial Systems Model',
    'hostname': 'cheyenne',
    'username': 'oleson',
    'version': 'cesm2_3_alpha02c',
    'revision_id': '$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $',
    'case_title': 'UNSET',
    'case_id': 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist',
    'Surface_dataset': 'surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr1850_c190214.nc',
    'Initial_conditions_dataset': 'finidat_interp_dest.nc',
    'PFT_physiological_constants_dataset': 'clm50_params.c210217.nc',
    'time_period_freq': 'month_1',
    'Time_constant_3Dvars_filename': './clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist.clm2.h3.1850-01-01-00000.nc',
    'Time_constant_3Dvars': 'ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY',
}

variable_metadata = {
    'GSSUNLN': {
        'long_name': 'sunlit leaf stomatal conductance at local noon',
        'units': 'umol H20/m2/s',
        'cell_methods': 'time: mean',
    },
    'GPP': {
        'long_name': 'gross primary production',
        'units': 'gC/m^2/s',
        'cell_methods': 'time: mean',
    },
    'FCTR': {
        'long_name': 'canopy transpiration',
        'units': 'W/m^2',
        'cell_methods': 'time: mean',
    },
    'TLAI': {
        'long_name': 'total projected leaf area index',
        'units': 'm^2/m^2',
        'cell_methods': 'time: mean',
    },
    'WUE': {
        'long_name': 'water use efficiency',
        'units': 'molC/molH2O',
        'description': 'GPP/GSSUNLN',
        'cell_methods': 'time: mean',
    },
    'BTRANMN': {
        'long_name': 'daily minimum of transpiration beta factor',
        'units': 'unitless',
        'cell_methods': 'time: mean',
    },
    'TOTVEGC': {
        'long_name': 'total vegetation carbon, excluding cpool',
        'units': 'gC/m^2',
        'cell_methods': 'time: mean',
    },
    'TSA': {
        'long_name': '2m air temperature',
        'units': 'K',
        'cell_methods': 'time: mean',
    },
    'VPD_CAN': {
        'long_name': 'canopy vapor pressure deficit',
        'units': 'kPa',
        'cell_methods': 'time: mean',
    },
    'RAIN_FROM_ATM': {
        'long_name': 'atmospheric rain received from atmosphere (pre-repartitioning)',
        'units': 'mm/s',
        'cell_methods': 'time: mean',
    },
    'SOILLIQ': {
        'long_name': 'soil liquid water (natural vegetated and crop landunits only)',
        'units': 'kg/m2',
        'cell_methods': 'time: mean',
    },
    'RH2M': {
        'long_name': '2m relative humidity',
        'units': '%',
        'cell_methods': 'time: mean',
    }
}

hist_type = 'h0'
case_name = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'

if hist_type == 'h1':
    main_directory_in = '/glade/work/bbuchovecky/WUE_analysis'
if hist_type == 'h0':
    main_directory_in = '/glade/campaign/cgd/tss/common/Land_Only_Simulations/CTSM51_DEV/CLM50_CTSM51_LAND_ONLY_RELEASE'

main_directory_out = '/glade/work/bbuchovecky/WUE_analysis'

variable_list = [
    # 'GPP',
    # 'GSSUNLN',
    # 'FCTR',
    'TLAI',
    # 'TSA',
    # 'WUE',
    # 'BTRANMN',
    # 'VPD_CAN',
    # 'RAIN_FROM_ATM',
    # 'SOILLIQ',
    # 'RH2M',
    # 'TSA',
]

time_period_in_years = [
    [1901, 2014],
    [1901, 1964],
    [1965, 2014],
]

for var in variable_list:
    print(var)

    # Load variable dataset
    if hist_type == 'h1':
        ds = xc.open_dataset(f'{main_directory_in}/timeseries_from_gridded/{case_name}.clm2.{hist_type}.{var}.185001-201412_gridded_annts.nc')
        ds = ds[var]
    if hist_type == 'h0':
        ds = xc.open_dataset(f'{main_directory_in}/{case_name}/lnd/proc/tseries/month_1/{case_name}.clm2.{hist_type}.{var}.185001-201412.nc')
        ds = format_ds_coords(ds)
        ds = ds.sel(time=slice('1901-01', '2014-12'))[var]
        ds = calculate_annual_timeseries(ds)

    for start_year, end_year in time_period_in_years:
        # Select variable and time period
        da = ds.sel(year=slice(start_year, end_year))

        # Compute linear trend
        print('-> computing linear trend', end=' ')
        start = time.time()
        coeff = da.polyfit(dim='year', deg=1, full=True, skipna=False)
        end = time.time()
        print(f'{(end-start): 0.4f}s')

        # Format linear trend
        coeff['polyfit_coefficients'].attrs = variable_metadata[var]
        coeff['polyfit_coefficients'].attrs['long_name'] = 'polyfit coefficients for '+coeff['polyfit_coefficients'].attrs['long_name']
        coeff['polyfit_coefficients'].attrs['description'] = f'polyfit of degree 1 over the year dimension {start_year}-{end_year} on the annual mean timeseries'
        coeff.attrs = global_metadata

        # Save linear trend to NetCDF
        coeff.to_netcdf(f'{main_directory_out}/trends/gridded/{var}_coeff.annts.{start_year}-{end_year}.nc')

        # Compute predicted values
        print('-> computing predicted values', end=' ')
        start = time.time()
        predicted = xr.polyval(da.year, coeff.polyfit_coefficients.chunk({'lon': 100, 'lat': 100}))
        end = time.time()
        print(f'{(end-start): 0.4f}s')

        # Format predicted values
        predicted = predicted.to_dataset(name=f'{var}_pred')
        predicted[f'{var}_pred'].attrs = variable_metadata[var]
        predicted[f'{var}_pred'].attrs['long_name'] = 'predicted '+predicted[f'{var}_pred'].attrs['long_name']
        predicted[f'{var}_pred'].attrs['description'] = f'polyfit of degree 1 over the year dimension {start_year}-{end_year} on the annual mean timeseries'
        predicted.attrs = global_metadata

        # Save predicted values to NetCDF
        predicted.to_netcdf(f'{main_directory_out}/trends/gridded/{var}_pred.annts.{start_year}-{end_year}.nc')