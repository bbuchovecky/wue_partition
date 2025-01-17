import numpy as np
import sparse
import xarray as xr
import time

## Functions copied from https://ncar.github.io/esds/posts/2022/sparse-PFT-gridding/#community-land-model-clm-output
def to_sparse(data, vegtype, jxy, ixy, shape):
    """
    Takes an input numpy array and converts it to a sparse array.

    Parameters
    ----------
    data: numpy.ndarray
        1D or 2D Data stored in compressed form.
    vegtype: numpy.ndarray

    jxy: numpy.ndarray
        Latitude index
    ixy: numpy.ndarray
        Longitude index
    shape: tuple
        Shape provided as sizes of (vegtype, jxy, ixy) in uncompressed
        form.

    Returns
    -------
    sparse.COO
        Sparse nD array
    """
    import sparse

    # This constructs a list of coordinate locations at which data exists
    # it works for arbitrary number of dimensions but assumes that the last dimension
    # is the "stacked" dimension i.e. "pft"
    if data.ndim == 1:
        coords = np.stack([vegtype, jxy - 1, ixy - 1], axis=0)
    elif data.ndim == 2:
        # generate some repeated time indexes
        # [0 0 0 ... 1 1 1... ]
        itime = np.repeat(np.arange(data.shape[0]), data.shape[1])
        # expand vegtype and friends for all time instants
        # by sequentially concatenating each array for each time instants
        tostack = [np.concatenate([array] * data.shape[0]) for array in [vegtype, jxy - 1, ixy - 1]]
        coords = np.stack([itime] + tostack, axis=0)
    else:
        raise NotImplementedError

    return sparse.COO(
        coords=coords,
        data=data.ravel(),
        shape=data.shape[:-1] + shape,
        fill_value=np.nan,
    )


def convert_pft_variables_to_sparse(dataset, pftnames):
    """
    Convert 2D PFT variables in dataset to 4D sparse arrays.

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset with DataArrays that have a `pft` dimension.

    Returns
    -------
    xarray.Dataset
        Dataset whose "PFT" arrays are now sparse arrays
        with `pft` dimension expanded out to (type, lat, lon)
    """

    import sparse
    import xarray as xr

    # extract PFT variables
    pfts = xr.Dataset({k: v for k, v in dataset.items() if "pft" in v.dims})

    # extract coordinate index locations
    ixy = dataset.pfts1d_ixy.astype(int)
    jxy = dataset.pfts1d_jxy.astype(int)
    vegtype = dataset.pfts1d_itype_veg.astype(int)
    npft = len(pftnames.data)

    # expected shape of sparse arrays to pass to `to_sparse` (excludes time)
    output_sizes = {
        "vegtype": npft,
        "lat": dataset.sizes["lat"],
        "lon": dataset.sizes["lon"],
    }

    result = xr.Dataset()
    # we loop over variables so we can specify the appropriate dtype
    for var in pfts:
        result[var] = xr.apply_ufunc(
            to_sparse,
            pfts[var],
            vegtype,
            jxy,
            ixy,
            kwargs=dict(shape=tuple(output_sizes.values())),
            input_core_dims=[["pft"]] * 4,
            output_core_dims=[["vegtype", "lat", "lon"]],
            dask="parallelized",
            dask_gufunc_kwargs=dict(
                meta=sparse.COO(np.array([], dtype=pfts[var].dtype)),
                output_sizes=output_sizes,
            ),
            keep_attrs=True,
        )

    # copy over coordinate variables lat, lon
    result = result.update(dataset[["lat", "lon"]])
    result["vegtype"] = pftnames.data
    # save the dataset attributes
    result.attrs = dataset.attrs
    return result

##

# File management
data_dir = '/glade/campaign/cgd/tss/common/Land_Only_Simulations/CTSM51_DEV/CLM50_CTSM51_LAND_ONLY_RELEASE/'
case_name = 'clm50_cesm23a02cPPEn08ctsm51d030_1deg_GSWP3V1_hist'
variables = [
    # 'GSSUNLN',
    # 'GPP',
    # 'TLAI',
    # 'FCTR',
    # 'TSA',
]

do_area = True
if do_area:
    variables = ['TSA', 'TSA']
    area_variables = ['area', 'landfrac']

# Load PFT names
pft_constants = xr.open_dataset('/glade/p/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm50_params.c210217.nc')
pftnames = pft_constants.pftname

for i, var in enumerate(variables):
    if do_area:
        var = area_variables[i]
    print(var)

    # Open dataset with PFT-level data
    ds = xr.open_dataset(
        data_dir+case_name+'/lnd/proc/tseries/month_1/'+case_name+'.clm2.h1.'+var+'.185001-201412.nc',
        decode_times=True,
    )
    ds_var_attrs = ds[var].attrs

    # Create PFT dimension array
    pftdim = ds.pft
    npftdim = pftdim.size

    # Create PFT type array
    npft = ds.pfts1d_itype_veg.max().astype(int)
    pftlist = np.arange(0, npft+1)

    # Convert 1d PFT-level array to gridded PFT-level array
    start = time.time()
    sparse_gridded_ds = convert_pft_variables_to_sparse(ds, pftlist)
    end = time.time()
    print('time to generate sparse array:', (end-start), 's')

    # Assign PFT name coordinate
    sparse_gridded_ds = sparse_gridded_ds.assign_coords({'vegtype_name': pftnames.data})

    # Create dense DataArray
    start = time.time()
    dense_data = sparse_gridded_ds[var].data.todense()
    dense_da = xr.DataArray(
        data=dense_data,
        dims=['time', 'vegtype', 'lat', 'lon'],
        coords=dict(
            time=sparse_gridded_ds.time,
            vegtype=sparse_gridded_ds.vegtype,
            vegtype_name=(('vegtype'), pftnames.values),
            lat=sparse_gridded_ds.lat,
            lon=sparse_gridded_ds.lon,
        )
    )
    end = time.time()
    print('time to densify sparse array:', (end-start), 's')

    # Format dense DataArray
    dense_da = dense_da.rename(var)
    dense_da.attrs = ds_var_attrs

    # Save as NetCDF file
    print('/glade/work/bbuchovecky/WUE_analysis/'+case_name+'.clm2.h1.'+var+'.185001-201412_gridded.nc')
    dense_da.to_netcdf('/glade/work/bbuchovecky/WUE_analysis/'+case_name+'.clm2.h1.'+var+'.185001-201412_gridded.nc')
