"""
Interpolate the cbec AllPoints data onto smooth surfaces
generates:
 - north_channel_and_ditch-1m.tif
 - lagoon-1m.tif

"""
from stompy.grid import unstructured_grid
import stompy.grid.quad_laplacian as quads
from stompy.spatial import wkb2shp, interp_orthogonal
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap
from stompy.spatial import field
from stompy import utils
import numpy as np
import six
import os
import pandas as pd

turbo=scmap.load_gradient('turbo.cpt')
##

fig_dir="figs-20210820"
os.path.exists(fig_dir) or os.makedirs(fig_dir)

## 
all_points=wkb2shp.shp2geom('../data/cbec/cbec-survey/All_Points.shp',
                            target_srs='EPSG:26910')

xyz=np.array( [ np.array(p) for p in all_points['geom']] )
xyz[:,2]=all_points['Z_m']

##

gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid15.pkl')


##
samples=pd.DataFrame(dict(x=xyz[:,0],
                          y=xyz[:,1],
                          value=xyz[:,2]))

sqg=quads.SimpleQuadGen(gen,cells=[1,2,3,4,5,6],nom_res=2.0)

##


from shapely import geometry
# There are some levee samples that don't get into the channel
# but are a real challenge for the anistropy.
exclude1=geometry.Polygon([[ 552331., 4124673.],
                           [ 552362., 4124678.],
                           [ 552489., 4124514.],
                           [ 552455., 4124483.],
                           [ 552359., 4124587.]])
exclude2=geometry.Polygon([[ 552575., 4124461.],
                           [ 552570., 4124434.],
                           [ 552490., 4124455.],
                           [ 552501., 4124483.]])

# But I also want to drop any lateral diffusivity in this reach
oink=interp_orthogonal.OrthoInterpolator(sqg,samples,overrides=[ (exclude1,[1,0]), 
                                                                  (exclude2,[1,0])])

oink.plot_result()
oink.grid.plot_edges(color='k',lw=0.4,alpha=0.3)

fld=oink.rasterize(dx=1.0,dy=1.0)

fld.write_gdal('north_channel_and_ditch-1m.tif')

##

sqg=quads.SimpleQuadGen(gen,cells=[9,11,12],nom_res=5)
oink=interp_orthogonal.OrthoInterpolator(sqg,samples)
oink.plot_result()
oink.grid.plot_edges(color='k',lw=0.4,alpha=0.3)

fld=oink.rasterize(dx=1.0,dy=1.0)

# Before this was used as is:
#fld.write_gdal('lagoon-1m.tif')

##

# That approach for the lagoon yields an overly narrow,
# shallow lagoon.
# Bring parts of the UCB model
epmp=unstructured_grid.RgfGrid(grd_fn="../data/ucb-model/MouthVariation/EpMp/PDO_EpMp.grd")
# Rasterize to match the lagoon fld from above:
epmp_fld=epmp.node_values_to_field('depth_node')
epmp_dem=epmp_fld.extract_tile(match=fld)
epmp_dem.F *=-1 # depth => elevation

##
if 1:
    plt.figure(3).clf()
    fig,axs=plt.subplots(1,3,num=3)
    img0=fld.plot(ax=axs[0],cmap=turbo,vmin=-1,vmax=2)
    img1=epmp_dem.plot(ax=axs[1],cmap=turbo,vmin=-1,vmax=2)

    delta=fld.copy()
    delta.F-=epmp_dem.F

    img2=delta.plot(ax=axs[2],cmap='coolwarm',vmin=-1,vmax=1)

    for ax in axs:
        plt.colorbar(ax.images[0],ax=ax,orientation='horizontal')
        ax.axis('off')
    axs[0].set_title("Interpolated cbec survey")
    axs[1].set_title("UCB Delft3D grid")
    axs[2].set_title("cbec - UCB")
    #fig.savefig("cbec_vs_UCB.png")

## 
fld.F=np.where( epmp_dem.F<fld.F, epmp_dem.F, fld.F)

##

fld.write_gdal('lagoon-1m.tif')

##---

# Southeast marsh channel:
# Used to be separate, as marsh_south_channel-interp-1m.tif, but now
# rolled into north_channel_and_ditch-1m.tif

