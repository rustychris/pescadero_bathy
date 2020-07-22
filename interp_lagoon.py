from stompy.grid import unstructured_grid
import stompy.grid.quad_laplacian as quads
from stompy.spatial import wkb2shp
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap
from stompy.spatial import field
from stompy import utils

import xarray as xr
import six

turbo=scmap.load_gradient('turbo.cpt')
## 
six.moves.reload_module(unstructured_grid)
six.moves.reload_module(quads)

gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid7.pkl')
qgs=quads.QuadsGen(gen,cells=[1])

qgs.plot_result()

##

ds=xr.open_dataset('lagoon-quads.nc')
ds.close()

g=unstructured_grid.UnstructuredGrid.read_ugrid('lagoon-quads.nc')

##

all_points=wkb2shp.shp2geom('../data/cbec/cbec-survey/All_Points.shp',
                            target_srs='EPSG:26910')


xyz=np.array( [ np.array(p) for p in all_points['geom']] )
xyz[:,2]=all_points['Z_m']

##

zoom=(552131, 552509., 4124173., 4124631.)

as_built_dem=field.GdalGrid('../data/cbec/to_RCD/asbuilt/merged.tif',
                            geo_bounds=zoom)

## 

# HERE
#  Boundary nodes of g get data from the raster.
#  The rest of AllPoints data takes an IJ from the grid.

# Do some version of anisotropic inverse distance weighting in IJ
# space. ==> refer to HOR ADCP interpolation code.

# Pick up samples along the boundary from the DEM

boundaries=g.boundary_linestrings()

xy_boundary=boundaries[0]

z_boundary=as_built_dem(xy_boundary)
xyz_boundary=np.c_[ xy_boundary,z_boundary]

all_xyz=np.concatenate( [ xyz_boundary,
                          xyz],axis=0 )


##

plt.figure(1).clf()
g.plot_edges(color='k',lw=0.5)

scat=plt.scatter(all_xyz[:,0],all_xyz[:,1],20,all_xyz[:,2],cmap=turbo)
img=as_built_dem.plot(cmap=turbo,alpha=0.3)

img.set_clim(scat.get_clim())

## 
# Aniso interpolation:

# Inputs:
#   grid
#   phi field defined on grid nodes.
#   psi field defined on grid nodes.

# samples in x,y,z
samp_xy=all_xyz[:,:2]
samp_z=all_xyz[:,2]

phi=g.nodes['ij'][:,0]
psi=g.nodes['ij'][:,1]
tri=g.mpl_triangulation()

from matplotlib.tri import LinearTriInterpolator

#from scipy.interpolate import LinearNDInterpolator

phi_interp=LinearTriInterpolator(tri,phi)
psi_interp=LinearTriInterpolator(tri,psi)


##
samp_phi= phi_interp(samp_xy[:,0],samp_xy[:,1])
samp_psi= psi_interp(samp_xy[:,0],samp_xy[:,1])

##

tgt_x=np.arange(zoom[0],zoom[1],1.0)
tgt_y=np.arange(zoom[2],zoom[3],1.0)
tgt_X,tgt_Y=np.meshgrid(tgt_x,tgt_y)

tgt_phi=phi_interp(tgt_X,tgt_Y)
tgt_psi=psi_interp(tgt_X,tgt_Y)

##

plt.figure(2).clf()
fig,axs=plt.subplots(1,2,num=2)
for ax,scal in zip(axs,[samp_phi,samp_psi]):
    g.plot_edges(color='k',lw=0.5,ax=ax)
    scat=ax.scatter(samp_xy[:,0],samp_xy[:,1],20,scal,cmap=turbo)

##

plt.figure(3).clf()
fig,ax=plt.subplots(num=3)

ax.scatter(samp_phi,samp_psi,30,samp_z,cmap=turbo)

##

# 0.2 seemed streaky in the longitudinal direction.
# 0.4 allows some weird fills on the south.
from scipy.interpolate import griddata
psi_scale=0.25
valid=np.isfinite(samp_phi*samp_psi*samp_z)

gridded=griddata( np.c_[ samp_phi, psi_scale*samp_psi][valid],
                  samp_z[valid],
                  np.c_[tgt_phi.ravel(),psi_scale*tgt_psi.ravel()] )
gridded=gridded.reshape( tgt_X.shape )


z_fld=field.SimpleGrid( F=gridded,
                        extents=zoom )

z_fld.smooth_by_convolution(kernel_size=5,iterations=3)
## 
plt.figure(4).clf()

img=z_fld.plot(cmap=turbo)
scat=plt.scatter(all_xyz[:,0],all_xyz[:,1],20,all_xyz[:,2],cmap=turbo)

## 
scat.set_clim(img.get_clim())
plt.axis(zoom)


## 

z_fld.write_gdal('lagoon-interp-1m.tif')

