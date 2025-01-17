import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc


def format_ds_coords(ds):
    """
    Formats the coordinates of an xarray Dataset.
    
    This function performs the following operations on the input Dataset:
        1. Adds bounds to the 'X' coordinate.
        2. Adds bounds to the 'Y' coordinate.
        3. Swaps the longitude axis to the range (-180, 180).
    
    Parameters:
    ds (xarray.Dataset):
        The input Dataset to be formatted.
    
    Returns:
    ds_formatted (xarray.Dataset):
        The formatted Dataset with updated coordinates.
    """
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds_formatted = xc.swap_lon_axis(ds, (-180, 180))
    return ds_formatted


def select_sites_from_gridded_data(xr_grid_data, df_site_data):
    """
    Selects grid data for specific sites from a gridded dataset.

    Parameters:
    xr_grid_data (xarray.Dataset or xarray.DataArray):
        The gridded dataset from which to select data.
    df_site_data (pandas.DataFrame):
        A DataFrame containing site information with columns
        'lat' and 'lon' for latitude and longitude.

    Returns:
    xr_site_data (xarray.DataArray):
        An xarray DataArray containing the selected site data,
        with a new dimension 'site' corresponding to each site in df_site_data.
    """
    
    # Get the number of sites and create an empty array to store the selected data
    nsite = df_site_data.iloc[:,0].size
    np_site_data = np.empty((nsite), dtype=xr.DataArray)

    # Select the grid box nearest to the coordinates for each site
    for i, row in df_site_data.iterrows():
        np_site_data[i] = xr_grid_data.sel(lat=row['lat'], lon=row['lon'], method='nearest')

    # Concatenate the selected data along a new 'site' dimension
    xr_site_data = xr.concat(np_site_data, dim='site')
    xr_site_data = xr_site_data.assign_coords({'site': np.arange(nsite)})
    xr_site_data['site'].attrs = {'long_name': 'tree ring site, arbitrary numbering from csv file'}

    return xr_site_data


global_metadata = {
    'latlon_site_file': '/glade/u/home/bbuchovecky/projects/wue_trend/lat_lon_pft_all.csv',
    'description': 'Gridded to (time,lat,lon,pft) then only grid boxes corresponding to tree ring sites were selected. All global attributes below are copied from the original history file.',
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
    }
}

# Load tree ring site data
site_data = pd.read_csv('./lat_lon_pft_all.csv', usecols=['lat', 'lon', 'PFT'])

# Load CLM5 gridded data
directory = '/glade/work/bbuchovecky/WUE_analysis'
casename = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'
variables = ['GSSUNLN', 'GPP', 'FCTR', 'TLAI']
tperiod = slice('1901-01', '2014-12')
sitepft = sorted(site_data['PFT'].unique())

for var in variables:
    # Load gridded variable and format coordinates
    gridded = xc.open_dataset(f'{directory}/gridded/{casename}.clm2.h1.{var}.185001-201412_gridded.nc')
    gridded = format_ds_coords(gridded)

    # Format PFT coordinates
    gridded = gridded.rename({'vegtype': 'pft', 'vegtype_name': 'pft_name'})
    gridded['pft'].attrs = {'long_name': 'plant functional type'}
    gridded['pft_name'].attrs = {'long_name': 'plant functional type name'}

    # Select PFTs and time period
    gridded = gridded[var].sel(pft=sitepft, time=tperiod)

    # Select sites from gridded variable
    indvsites = select_sites_from_gridded_data(gridded, site_data)

    # Iterate through PFTs
    for pft in sitepft:
        # Select site indices with PFT
        pft_site = site_data[site_data['PFT'] == pft]

        # Select corresponding sites
        pft_data = indvsites.sel(site=pft_site.index, pft=pft)
        pft_data = pft_data.drop('pft')

        # Convert to dataset
        pft_data = pft_data.to_dataset(name=var)

        # Add metadata
        pft_data.attrs = global_metadata
        pft_data[var].attrs = variable_metadata[var]

        # Save to NetCDF
        pft_data.to_netcdf(f'{directory}/for/marja/clm50_1deg_GSWP3V1.{var}.190101-201412.sites.pft{pft}.nc')

    # Convert to dataset
    indvsites = indvsites.to_dataset(name=var)

    # Add metadata
    indvsites.attrs = global_metadata
    indvsites[var].attrs = variable_metadata[var]

    # Save to NetCDF
    indvsites.to_netcdf(f'{directory}/for/marja/clm50_1deg_GSWP3V1.{var}.190101-201412.sites.pftall.nc')
