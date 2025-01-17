import time
import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc

def open_files_with_overlap(glob_path, select_files=0):
    """
    Opens multiple files into a single dataset correcting
    for overlapping time coordinate dimensions
    """
    import glob

    # Get the files that match the glob path - sorting is crucial!
    matching_files = np.array(sorted(glob.glob(glob_path)))[select_files:]

    # Create empty array to hold the start/end year for each file
    file_times = np.empty((len(matching_files), 2), dtype=int)

    # Iterate through each file and save the start/end year
    for i, file in enumerate(matching_files):
        file_times[i,:] = [ int(x[:4]) for x in file.split('.')[-2].split('-') ]

    # Flatten and shift the array of years to check for overlap
    file_times_flat = file_times.flatten()
    overlap_times = file_times_flat[2::2] - file_times_flat[1:-1:2]
    overlap_files_index = np.where(overlap_times < 0)[0]

    # Return the time slice to eliminate the overlap
    if overlap_files_index.size > 0:
        print(f'{overlap_files_index.size} overlapping file(s)')
        correct_end_years = file_times[overlap_files_index][:,1] + overlap_times[overlap_files_index]# - 1  # Note: these runs start on Feb!!
        correct_end_times = [ f'{x:05.0f}-01'[1:] for x in correct_end_years ]  # Note: these runs start on Feb!!

        # Open each file and select the correct time periods to eliminate overlap
        these_ds = []
        overlap_index_counter = 0
        for i, file in enumerate(matching_files):
            ds = xr.open_dataset(file)
            if i in overlap_files_index:
                ds = ds.sel(time=slice(None, correct_end_times[overlap_index_counter]))
                overlap_index_counter += 1
            these_ds.append(ds)
        # Concatenate the individual datasets into one
        return xr.concat(these_ds, dim='time')

    # Open normally!
    return xc.open_mfdataset(glob_path)  

def check_time_concat(da):
    """
    Checks that the number of months in the DataArray is equal
    to the time dimension. Only works for monthly timeseries.
    """
    start_date = [ int(x) for x in da.time[ 0].item().strftime().split(' ')[0].split('-') ]
    end_date   = [ int(x) for x in da.time[-1].item().strftime().split(' ')[0].split('-') ]
    length_in_months = ( (end_date[0] - start_date[0]) * 12 ) + (end_date[1] - start_date[1]) + 1
    return da.time.size == length_in_months

def format_ds_coords(ds):
    """
    Formats a Dataset by centering the time coordinate to midpoints,
    adding bounds for the x-axis (lon) and y-axis (lat), and shifting
    the longitude coordinates to [-180,180].
    """
    ds = xc.center_times(ds)
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

def get_medslope_runs_var(this_var, case, dom='lnd', freq='month_1', htype='h0', is_coupled=True, remove_first_years=0, select_files=0):
    """
    Add docstring
    """
    if is_coupled:
        maindir = '/glade/campaign/cgd/tss/people/czarakas/StomatalSlope2021/coupled_experiments'
    if not is_coupled:
        maindir = '/glade/campaign/cgd/tss/people/czarakas/StomatalSlope2021/offline_experiments'    

    # Open file(s) with variable of interest

    ds = open_files_with_overlap(f'{maindir}/{case}_0*/{dom}/proc/tseries/{freq}/{case}_*.{htype}.{this_var}*.nc', select_files=select_files)
    ds = format_ds_coords(ds)

    # Select the variable of interest
    da = ds[this_var]
    assert check_time_concat(da)

    # Select the equilibrium period
    da = da.isel(time=slice(remove_first_years*12, None))

    return da

def get_medslope_runs_vardict(this_var):
    """
    Add docstring
    """
    this_medslope = { 'DEF':{}, 'LOW':{}, 'HIGH':{} }
    for this_slope, this_slope_dict in this_medslope.items():
        print(this_slope)

        print(' -> 1xCO2')
        select_files = 0
        if this_slope == 'HIGH':
            select_files = 1

        this_slope_dict['1xCO2'] = get_medslope_runs_var(
            this_var, f'coupled_{this_slope}medslope_1xCO2',
            select_files=select_files,
            remove_first_years=40)


        print(' -> 2xCO2')
        this_slope_dict['2xCO2'] = get_medslope_runs_var(
            this_var, f'coupled_{this_slope}medslope_2xCO2',
            remove_first_years=40)

    return this_medslope

def select_sites_for_medslope_runs(xr_grid_data, df_site_data):
    """
    Add docstring
    """
    nsite = df_site_data.iloc[:,0].size
    np_site_data = np.empty((nsite), dtype=xr.DataArray)

    # Select the grid box nearest to the coordinates for each site
    for i, row in df_site_data.iterrows():
        gridbox_timeseries = xr_grid_data.sel(lat=row['lat'], lon=row['lon'], method='nearest')
        np_site_data[i] = gridbox_timeseries

    xr_site_data = xr.concat(np_site_data, dim='site')
    xr_site_data = xr_site_data.assign_coords({'site': np.arange(nsite)})

    return xr_site_data

# Load tree ring coordinates
tree_ring_coords = pd.read_csv(
    '/glade/u/home/bbuchovecky/projects/wue_partition/latlon_gym_ang.csv',
    usecols=[0,1])

variables = ['GSSUNLN', 'GPP', 'FCTR']

for var in variables:
    print(f'loading {var}')
    medslope = get_medslope_runs_vardict(var)

    medslope_site = {}
    for slope, slope_dict in medslope.items():
        print(f'selecting 1xCO2 and 2xCO2 for {slope}')
        medslope_site_1xco2 = select_sites_for_medslope_runs(medslope[slope]['1xCO2'], tree_ring_coords)
        medslope_site_2xco2 = select_sites_for_medslope_runs(medslope[slope]['2xCO2'], tree_ring_coords)

        medslope_site_1xco2 = medslope_site_1xco2.mean(dim='time')
        medslope_site_2xco2 = medslope_site_2xco2.mean(dim='time')

        print(f'computing delta for {slope}')
        medslope_site_delta = medslope_site_2xco2 - medslope_site_1xco2
        
        start = time.time()
        medslope_site_delta.to_netcdf(f'/glade/work/bbuchovecky/WUE_analysis/medslope/{var}.coupled_{slope}medslope_2xCO2-1xCO2.nc')
        end = time.time()

        print(f'{end-start: 0.5} s')
        print(f'/glade/work/bbuchovecky/WUE_analysis/medslope/{var}.coupled_{slope}medslope_2xCO2-1xCO2.nc')
