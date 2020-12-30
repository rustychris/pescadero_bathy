"""
Interpolate the cbec AllPoints data onto smooth surfaces
generates:
 - north_channel_and_ditch-1m.tif
 - lagoon-1m.tif

"""
from stompy.grid import unstructured_grid
import stompy.grid.quad_laplacian as quads
from stompy.spatial import wkb2shp
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap
from stompy.spatial import field
from stompy import utils
import numpy as np
import six
import os

turbo=scmap.load_gradient('turbo.cpt')
##

fig_dir="figs-20201230"
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

oink=interp_orthogonal.OrthoInterpolator(sqg,samples)

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
fld.write_gdal('lagoon-1m.tif')

##---

# Southeast marsh channel:
# Used to be separate, as marsh_south_channel-interp-1m.tif, but now
# rolled into north_channel_and_ditch-1m.tif

