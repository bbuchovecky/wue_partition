import time
import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc
import scipy.stats as stats

## Functions
def format_ds_coords(ds):
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

def timeseries_gridded_pft_linregress(gridded_data, time_period, pft_subset, pft_subset_is_ang, get_timing=False):
    """
    site_data:   DataArrays for the given variable
    time_period: a slice object for the time period to perform the linear regression over
    pft_subset:  the DataArray of PFTs, a subset of the full set of PFTs (only trees in this case)
    """

    if get_timing:
        start = time.time()

    # Select the data over the specified time period
    time_subset_data = this_site.sel(time=time_period)

    # Create a time array for the linear regression
    month = np.arange(0, len(time_subset_data['time']), 1)

    # Iterate through the specified PFTs
    for j, this_pft in enumerate(pft_subset):
        # Check if the PFT exists for the given site
            if ~np.isnan(time_subset_data.sel(vegtype=j).isel(time=0)):
                # Update the site-level PFT mask to indicate that the PFT exists
                site_pft_mask[i,j] = 1

                # Perform the linear regression
                lr = stats.linregress(
                    month,
                    time_subset_data.sel(vegtype=j)
                )

                # Save the linear regression statistics
                site_linregress[i,j,0] = lr.slope * 12 * 10  # Convert slope from X/1month to X/10yr
                site_linregress[i,j,1] = lr.rvalue
                site_linregress[i,j,2] = lr.pvalue

    # (lat, lon, vegtype, stats)

    # create DataArray enumerating months from start date in shape of gridded_data
    # perform xs.linslope(a=gridded_month, b=gridded_data)
    # perform xs.pearson_r(a=gridded_month, b=gridded_data)
    # perform xs.pearson_r_p_value(a=gridded_month, b=gridded_data)

    # concatenate stats into a single DataArray, along new stats dimension


## Variables to compute trends over
variables = [
    'GSSUNLN',
    'GPP',
    'FCTR',
    'TLAI',
    'TSA',
    'WUE',
]

## File management
directory = '/glade/work/bbuchovecky/WUE_analysis/'
casename = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'

## Time periods of compute trends over
time_period = [
    ['1901-01', '2014-12'],
    ['1901-01', '1964-12'],
    ['1965-01', '2014-12'],
]

for var in variables:
    print(var)

    ## Load variable
    this_data = xc.open_dataset(directory+casename+'.clm2.h1.'+var+'.185001-201412_gridded.nc')
    this_data = format_ds_coords(this_data)
    this_data = this_data[var]

    ## Subselect sites
    site_data = select_site_from_gridded_data(this_data, tree_ring_coords, do_time=False)

    ## Select PFTs to use
    # needleleaf~gymnosperm and broadleaf~angiosperm

    # All non-crop PFTs
    noncrop_pft = this_data['vegtype'][1:15]
    noncrop_pft_names = this_data['vegtype_name'][1:15]
    noncrop_pft_is_ang = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    print('noncrop_pft_names:')
    print(noncrop_pft_names.values)

    ## Iterate through each time period
    for this_time_period in time_period:
        ## Create time slice
        this_start_time = this_time_period[0]
        this_end_time = this_time_period[1]
        this_time_slice = slice(this_time_period[0], this_time_period[1])
        print(f'{this_time_period[0]} {this_time_period[1]}')

        ## Perform linear regression
        print('/glade/work/bbuchovecky/WUE_analysis/'+var+'_linregress.'+this_start_time.replace('-', '')+'-'+this_end_time.replace('-', '')+'.nc')
        site_linregress, site_pft_mask = timeseries_site_pft_linregress(
            site_data=site_data,
            time_period=this_time_slice,
            pft_subset=noncrop_pft,
            pft_subset_is_ang=noncrop_pft_is_ang,
            get_timing=True)
        site_linregress.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/'+var+'_linregress.'+this_start_time.replace('-', '')+'-'+this_end_time.replace('-', '')+'.nc')
        site_pft_mask.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/'+var+'_pftmask.'+this_start_time.replace('-', '')+'-'+this_end_time.replace('-', '')+'.nc')
