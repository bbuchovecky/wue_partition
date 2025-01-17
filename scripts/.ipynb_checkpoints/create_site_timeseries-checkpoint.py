import numpy as np
import xarray as xr
import pandas as pd
import xcdat as xc

def format_ds_coords(ds):
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

def select_sites_from_gridded_data(xr_grid_data, df_site_data):
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

global_metadata = {
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
    # 'GSSUNLN': {
    #     'long_name': 'sunlit leaf stomatal conductance at local noon',
    #     'units': 'umol H20/m2/s',
    #     'cell_methods': 'time: mean',
    # },
    # 'GPP': {
    #     'long_name': 'gross primary production',
    #     'units': 'gC/m^2/s',
    #     'cell_methods': 'time: mean',
    # },
    # 'FCTR': {
    #     'long_name': 'canopy transpiration',
    #     'units': 'W/m^2',
    #     'cell_methods': 'time: mean',
    # },
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
    }
}

directory = '/glade/work/bbuchovecky/WUE_analysis'
casename = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'
variables = [
    # 'GSSUNLN',
    # 'GPP',
    # 'FCTR',
    'TLAI',
    'WUE',
]

# Load tree ring coordinates
tree_ring_coords = pd.read_csv('/glade/u/home/bbuchovecky/projects/wue_partition/latlon_gym_ang.csv', usecols=[0,1])

for this_var in variables:
    print(this_var)

    # Load gridded variable and format coordinates
    print('loading gridded variable')
    this_data = xc.open_dataset(f'{directory}/{casename}.clm2.h1.{this_var}.185001-201412_gridded.nc')
    this_data = format_ds_coords(this_data)
    this_data = this_data[this_var]

    # Select sites from gridded variable
    print('selecting sites')
    this_site_data = select_sites_from_gridded_data(this_data, tree_ring_coords)

    # Convert to dataset
    this_site_data = this_site_data.to_dataset(name=this_var)

    # Format coordinates
    this_site_data = this_site_data.rename({'vegtype': 'pft', 'vegtype_name': 'pft_name'})
    this_site_data['pft'].attrs = {'long_name': 'plant functional type'}
    this_site_data['pft_name'].attrs = {'long_name': 'plant functional type name'}
    this_site_data['site'].attrs = {'long_name': 'tree ring site'}

    # Add metadata
    this_site_data.attrs = global_metadata
    this_site_data[this_var].attrs = variable_metadata[this_var]

    # Save to NetCDF
    print('saving')
    this_site_data.to_netcdf(f'{directory}/sites_from_gridded/{casename}.clm2.h1.{this_var}.185001-201412_sites.nc')