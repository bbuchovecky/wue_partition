import xarray as xr
import xcdat as xc

def format_ds_coords(ds):
    ds = ds.bounds.add_bounds('X')
    ds = ds.bounds.add_bounds('Y')
    ds = xc.swap_lon_axis(ds, (-180, 180))
    return ds

time_slice = slice('1901-01', '2014-12')

min_lai_threshold = [
    0.5,
    0.25,
    0.1,
    0.01
]

pft_constants = xr.open_dataset('/glade/p/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm50_params.c210217.nc')
pftnames = pft_constants.pftname
pftnames = [ str(name.strip())[2:-1] for name in pftnames[0:15].values ]

# gridbox-level
main_directory = '/glade/campaign/cgd/tss/common/Land_Only_Simulations/CTSM51_DEV/CLM50_CTSM51_LAND_ONLY_RELEASE'
case_name = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'
tlai_h0 = xc.open_dataset(f'{main_directory}/{case_name}/lnd/proc/tseries/month_1/{case_name}.clm2.h0.TLAI.185001-201412.nc')
tlai_h0 = format_ds_coords(tlai_h0)
tlai_h0 = tlai_h0['TLAI'].sel(time=time_slice)

# PFT-level
main_directory = '/glade/work/bbuchovecky/WUE_analysis'
case_name = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'
tlai_h1 = xc.open_dataset(f'{main_directory}/{case_name}.clm2.h1.TLAI.185001-201412_gridded.nc')
tlai_h1 = tlai_h1.rename({'vegtype': 'pft', 'vegtype_name': 'pft_name'}).sel(pft=slice(0,14))
tlai_h1['pft_name'] = pftnames
tlai_h1 = format_ds_coords(tlai_h1)
tlai_h1 = tlai_h1['TLAI'].sel(time=time_slice)


# create mask 
# -> compare the maximum monthly LAI for each year to the threshold
# -> mask grid boxes where the LAI drops below the threshold at any year in the timeseries
# -> save to netcdf since the file is the computation takes a while
lai_masks = []
pftlai_masks = []
for min_lai in min_lai_threshold:
    print(min_lai)

    this_lai_mask = (tlai_h0.resample(time='Y').max() > min_lai).min(dim='time')
    lai_masks.append(this_lai_mask.rename(f'min_lai_{str(min_lai)}'))

    this_pftlai_mask = (tlai_h1.resample(time='Y').max() > min_lai).min(dim='time')
    pftlai_masks.append(this_pftlai_mask.rename(f'min_lai_{str(min_lai)}'))

lai_masks = xr.merge(lai_masks)
pftlai_masks = xr.merge(pftlai_masks)

lai_masks.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/lai_masks.nc')
pftlai_masks.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/pftlai_masks.nc')