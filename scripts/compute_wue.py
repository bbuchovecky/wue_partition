import numpy as np
import xarray as xr
import xcdat as xc

## Load GSSUNLN and GPP
print('loading variables')
def format_ds_coords(ds):
    """
    Add docstring
    """
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

gssunln = xc.open_dataset('/glade/work/bbuchovecky/WUE_analysis/clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist.clm2.h1.GSSUNLN.185001-201412_gridded.nc')
gssunln = format_ds_coords(gssunln)
gssunln = gssunln['GSSUNLN']

gpp = xc.open_dataset('/glade/work/bbuchovecky/WUE_analysis/clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist.clm2.h1.GPP.185001-201412_gridded.nc')
gpp = format_ds_coords(gpp)
gpp = gpp['GPP']

## Convert GPP from gC to molC
g_to_mol_C = 12.011  # g / mol
gpp = gpp / g_to_mol_C

## Convert GSSUNLN from umol H2O to mol H2O
umol_to_mol = 1e-6
gssunln = gssunln * umol_to_mol

## Compute WUE
print('computing WUE')
wue = gpp / gssunln
wue.name = 'WUE'
wue.attrs = {
    'long_name': 'water use efficiency',
    'units': 'molC/molH2O',
    'description': 'GPP/GSSUNLN'
}

## Save to NetCDF
print('saving WUE')
wue.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist.clm2.h1.WUE.185001-201412_gridded.nc')