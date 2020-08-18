from stompy.grid import unstructured_grid
import stompy.grid.quad_laplacian as quads
from stompy.spatial import wkb2shp
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap
from stompy.spatial import field
from stompy import utils
from matplotlib.tri import LinearTriInterpolator
from scipy.interpolate import griddata
import numpy as np
import xarray as xr
import six
import os

turbo=scmap.load_gradient('turbo.cpt')
##

fig_dir="figs-20200805"
os.path.exists(fig_dir) or os.makedirs(fig_dir)

## 
all_points=wkb2shp.shp2geom('../data/cbec/cbec-survey/All_Points.shp',
                            target_srs='EPSG:26910')


xyz=np.array( [ np.array(p) for p in all_points['geom']] )
xyz[:,2]=all_points['Z_m']

##

def ortho_interp(g,coords='ij',aniso=1.0,dx=1.0,dy=1.0,
                 plot_geo=False,plot_mapped=False):
    """
    Wrap up the orthogonal interpolation code for Pescadero
    g: a grid with g.nodes[coords] defined, which sets the
    orthogonal coordinate system.
    
    aniso: a global anisotropy factor, multiplicative with whatever
     anisotropy is baked into g.nodes[coords]

    dx, dy: resolution of the output raster

    plot_geo: Show the samples in geographic coordinates
    plot_mapped: Show the samples in mapped coordinates
    """
    zoom=g.bounds()

    as_built_dem=field.GdalGrid('../data/cbec/to_RCD/asbuilt/merged.tif',
                                geo_bounds=zoom)

    #  Boundary nodes of g get data from the raster.
    #  The rest of AllPoints data takes an IJ from the grid.

    # Pick up samples along the boundary from the DEM
    boundaries=g.boundary_linestrings()

    xy_boundary=boundaries[0]

    z_boundary=as_built_dem(xy_boundary)
    xyz_boundary=np.c_[ xy_boundary,z_boundary]

    all_xyz=np.concatenate( [ xyz_boundary,
                              xyz],axis=0 )

    # samples in x,y,z
    samp_xy=all_xyz[:,:2]
    samp_z=all_xyz[:,2]

    # Get the samples into psi/phi space
    phi=g.nodes[coords][:,0]
    psi=g.nodes[coords][:,1]
    # likey thereare better methods that utilize the properties of g
    # (quadrilateral, approx orthogonal)
    tri=g.mpl_triangulation()

    phi_interp=LinearTriInterpolator(tri,phi)
    psi_interp=LinearTriInterpolator(tri,psi)

    samp_phi=phi_interp(samp_xy[:,0],samp_xy[:,1])
    samp_psi=psi_interp(samp_xy[:,0],samp_xy[:,1])

    tgt_x=np.arange(zoom[0],zoom[1],dx)
    tgt_y=np.arange(zoom[2],zoom[3],dy)
    tgt_X,tgt_Y=np.meshgrid(tgt_x,tgt_y)

    tgt_phi=phi_interp(tgt_X,tgt_Y)
    tgt_psi=psi_interp(tgt_X,tgt_Y)

    if plot_geo:
        plt.figure(2).clf()
        fig,axs=plt.subplots(1,2,num=2)
        for ax,scal in zip(axs,[samp_phi,samp_psi]):
            g.plot_edges(color='k',lw=0.5,ax=ax)
            scat=ax.scatter(samp_xy[:,0],samp_xy[:,1],20,scal,cmap=turbo)

    if plot_mapped:
        plt.figure(3).clf()
        fig,ax=plt.subplots(num=3)
        ax.scatter(samp_phi,samp_psi,30,samp_z,cmap=turbo)

    valid=np.isfinite(samp_phi*samp_psi*samp_z)

    print("Running griddata")
    gridded=griddata( np.c_[ samp_phi, aniso*samp_psi][valid],
                      samp_z[valid],
                      np.c_[tgt_phi.ravel(),aniso*tgt_psi.ravel()] )
    print("Done")

    gridded=gridded.reshape( tgt_X.shape )

    z_fld=field.SimpleGrid( F=gridded,
                            extents=zoom )

    # ehh -- probably ought to split this out
    return z_fld

##
from stompy.grid import front
six.moves.reload_module(unstructured_grid)
six.moves.reload_module(front)
six.moves.reload_module(quads)

gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid9.pkl')

##

# HERE  - write out summary plots for these.

qg=quads.QuadGen(gen,cell=0,execute=True,anisotropic=False,nom_res=5)
# aniso: smaller is more longitudinal
z_fld=ortho_interp(qg.g_final,aniso=0.3)
z_fld.smooth_by_convolution(kernel_size=5,iterations=5)

plt.figure(5).clf()
z_fld.plot(cmap=turbo)

# z_fld.write_gdal('lagoon-interp-1m.tif')

##---

# HERE -This is now updated, re-rendered.  Try building the field again in compare-lidar.
gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid9.pkl')
qg=quads.QuadGen(gen,cell=1,final='triangle',nom_res=5,execute=True)

z_fld=ortho_interp(qg.g_final,aniso=2.0)
z_fld.smooth_by_convolution(kernel_size=3,iterations=3)

plt.figure(4).clf()
z_fld.plot(cmap=turbo)

z_fld.write_gdal('north_channel-interp-1m.tif')

##

qg=quads.QuadGen(gen,cell=2,anisotropic=False,nom_res=15)


##
z_fld=ortho_interp(qg.g_final,aniso=2.0)
z_fld.smooth_by_convolution(kernel_size=3,iterations=3)

plt.figure(6).clf()
z_fld.plot(cmap=turbo)

##

# Southeast marsh channel
gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid9.pkl')

qg=quads.QuadGen(gen,cell=3,anisotropic=False,nom_res=3)


##
z_fld=ortho_interp(qg.g_final,aniso=0.2)
z_fld.smooth_by_convolution(kernel_size=3,iterations=3)

plt.figure(6).clf()
z_fld.plot(cmap=turbo)

z_fld.write_gdal('marsh_south_channel-interp-1m.tif')
