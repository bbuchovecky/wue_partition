import sys
import time
import numpy as np
import xarray as xr

def calculate_annual_timeseries(da):
    """
    Calculates the annual timeseries, weighted by the number of days in each month
    """
    nyears = len(da.groupby('time.year'))
    month_length = da.time.dt.days_in_month

    weights = month_length.groupby('time.year') / month_length.astype(float).groupby('time.year').sum()        
    np.testing.assert_allclose(weights.groupby('time.year').sum().values, np.ones(nyears)) 

    return (da * weights).groupby('time.year').sum(dim='time')

main_directory = '/glade/work/bbuchovecky/WUE_analysis'
case_name = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'

variable_list = [
    # 'GPP',
    # 'GSSUNLN',
    # 'FCTR',
    # 'TLAI',
    # 'TSA',
    # 'WUE',
    'BTRANMN',
]

time_slice = slice('1851-01', '2014-12')

# Create slice objects to select the corresponding ang/gym tree PFTs
angtreepft_index_slice = slice(4,9)
gymtreepft_index_slice = slice(1,4)
natpft_index_slice = slice(0,14)

for variable in variable_list:
    print(variable)

    # Load gridded timeseries
    print('\tstart loading gridded monthly timeseries')
    start = time.time()
    file_path = f'{main_directory}/{case_name}.clm2.h1.{variable}.185001-201412_gridded.nc'
    monthly_ts = xr.open_dataset(file_path)
    monthly_ts = monthly_ts[variable]
    end = time.time()
    print(f'\tdone loading gridded monthly timeseries ({(end-start): 0.4f}s)')

    # Select PFTs corresponding to angiosperm and gymnosperm classifications
    # print('\tstart computing PFT mean')
    # start = time.time()
    # monthly_ts_ang = monthly_ts.sel(vegtype=angtreepft_index_slice).mean(dim='vegtype')
    # monthly_ts_gym = monthly_ts.sel(vegtype=gymtreepft_index_slice).mean(dim='vegtype')
    monthly_ts = monthly_ts.sel(vegtype=natpft_index_slice).rename({'vegtype': 'pft', 'vegtype_name': 'pft_name'})
    # end = time.time()
    # print(f'\tdone computing PFT mean ({(end-start): 0.4f}s)')

    # Compute annual mean
    print('\tstart computing annual mean')
    start = time.time()
    # annual_ts_ang = calculate_annual_timeseries(monthly_ts_ang.sel(time=time_slice))
    # annual_ts_gym = calculate_annual_timeseries(monthly_ts_gym.sel(time=time_slice))
    annual_ts = calculate_annual_timeseries(monthly_ts.sel(time=time_slice))
    end = time.time()
    print(f'\tdone computing annual mean ({(end-start): 0.4f}s)')

    # Add attributes
    # annual_ts_ang.attrs['description'] = 'average over natural broadleaf PFTs'
    # annual_ts_gym.attrs['description'] = 'average over natural needleleaf PFTs'

    # Add name
    # annual_ts_ang = annual_ts_ang.rename(variable)
    # annual_ts_gym = annual_ts_gym.rename(variable)
    annual_ts = annual_ts.rename(variable)

    # Save to NetCDF file
    print('\tstart saving to NetCDF')
    start = time.time()
    # output_file_path_ang = f'{main_directory}/timeseries_from_gridded/{case_name}.clm2.h1.{variable}.185001-201412_gridded_ang_annts.nc'
    # annual_ts_ang.to_netcdf(output_file_path_ang)

    # output_file_path_gym = f'{main_directory}/timeseries_from_gridded/{case_name}.clm2.h1.{variable}.185001-201412_gridded_gym_annts.nc'
    # annual_ts_gym.to_netcdf(output_file_path_gym)

    output_file_path = f'{main_directory}/timeseries_from_gridded/{case_name}.clm2.h1.{variable}.185001-201412_gridded_annts.nc'
    annual_ts.to_netcdf(output_file_path)

    end = time.time()
    print(f'\tdone saving to NetCDF ({(end-start): 0.4f}s)')