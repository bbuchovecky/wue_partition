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

def select_site_from_gridded_data(xr_data, df_site_data, do_time=True, do_coords=True):
    assert do_time or do_coords
    nsite = df_site_data.iloc[:,0].size
    np_site_data = np.empty((nsite), dtype=xr.DataArray)

    for i, row in df_site_data.iterrows():
        # Select the grid box nearest to the site coordinates
        if do_coords:
            gridbox_timeseries = xr_data.sel(lat=row['lat'], lon=row['lon'], method='nearest')

        # Select the time period of the site (if observational data)
        if do_time:
            gridbox_timeseries = gridbox_timeseries.sel(time=slice(row['start_date'], row['end_date']))

        np_site_data[i] = gridbox_timeseries
    
    return np_site_data

def timeseries_site_pft_linregress(site_data, time_period, pft_subset, pft_subset_is_ang, get_timing=False):
    """
    site_data:   the array of DataArrays for each site
    time_period: a slice object for the time period to perform the linear regression over
    pft_subset:  the DataArray of PFTs, a subset of the full set of PFTs (only trees in this case)
    """
    nsite = site_data.size
    npft = pft_subset.size
    nstats = 3

    site_linregress = np.zeros((nsite, npft, nstats))
    site_pft_mask = np.zeros((nsite, npft))

    for i, this_site in enumerate(site_data):
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

            # Fill with NaN if the PFT does not exist for the given site
            else:
                site_linregress[i,j,:] = [np.nan, np.nan, np.nan]
            
        if get_timing:
            end = time.time()
            print(f'time for site {i} with {site_pft_mask[i].sum()} active pfts:{end-start: 0.5} s')

    site_linregress = xr.DataArray(
        data=site_linregress,
        dims=['site', 'vegtype', 'stats'],
        coords=dict(
            site=np.arange(nsite),
            vegtype=pft_subset,
            vegtype_name=(('vegtype'), pft_subset['vegtype_name'].values),
            is_ang=(('vegtype'), pft_subset_is_ang),
            stats=np.arange(3),
            stats_descr=(('stats'), np.array(['slope [X/10yr]', 'rvalue', 'pvalue'])),
        ),
        name=site_data[0].name,
    )

    site_pft_mask = xr.DataArray(
        data=site_pft_mask,
        dims=['site', 'vegtype'],
        coords=dict(
            site=np.arange(nsite),
            vegtype=pft_subset,
            vegtype_name=(('vegtype'), pft_subset['vegtype_name'].values),
            is_ang=(('vegtype'), pft_subset_is_ang),
        ),
        name=str(site_data[0].name)+'_pftmask'
    )
    
    return site_linregress, site_pft_mask

## Load tree ring sites
tree_ring_coords = pd.read_csv('/glade/u/home/bbuchovecky/projects/wue_partition/latlon_gym_ang.csv', usecols=[0,1,2])

# Tree type classification: 0=gym, 1=ang, 2=none
tree_ring_coords = tree_ring_coords.rename(columns={'gym=0/ang=1': 'isAng'})
tree_ring_coords = tree_ring_coords.replace({'isAng': np.nan}, 2)

ang_coords = tree_ring_coords.loc[tree_ring_coords['isAng']==1].drop(columns='isAng')
gym_coords = tree_ring_coords.loc[tree_ring_coords['isAng']==0].drop(columns='isAng')

## Variables to compute trends over
variables = [
    # 'GSSUNLN',
    # 'GPP',
    # 'WUE',
    # 'FCTR',
    'TSA',
    'TLAI',
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
        print(this_time_period)

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
