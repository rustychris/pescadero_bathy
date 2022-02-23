# Load in Lidar and pre-compiled DEMs, particularly in N Marsh
#
import os,shutil
import pandas as pd
from stompy.spatial import field, wkb2shp
import six
import matplotlib.pyplot as plt
import numpy as np
import stompy.plot.cmap as scmap
from stompy.plot import plot_wkb
from stompy import utils

from stompy.spatial import interp_orthogonal
six.moves.reload_module(interp_orthogonal)

turbo=scmap.load_gradient('turbo.cpt')

turbo_gray=scmap.load_gradient('turbo.cpt')
turbo_gray.set_over( "0.8")

sst=scmap.load_gradient('oc-sst.cpt')
##

fig_dir="figs-20210913"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

##

version="existing" # "existing" or "asbuilt"
dredge=False
date_str="20210928"
render_05m=False # 

import sys

if __name__=='__main__' and not sys.argv[0].endswith('python'):
    import argparse
    parser = argparse.ArgumentParser(description='Compile DEM for Pescadero/Butano.')
    parser.add_argument("--version", help="Select asbuilt or existing", 
                        default=version)
    parser.add_argument("--dredge", help="Enable dredging the mouth",action='store_true')
    parser.add_argument("--date", help="Label for output, typically YYYYMMDD",
                        default=date_str)

    args=parser.parse_args()
    version=args.version
    date_str=args.date
    # print(args)

clip=zoom=(551800, 553290, 4124100, 4125400)

# cbec surfaces:
# For starters load just the N marsh and pond portion.
# For compilation at the end load the whole thing

options=""
if dredge:
    options+="_dredge"

if version=='asbuilt':
    cbec_dem_fn='../data/cbec/to_RCD/asbuilt/merged.tif'
    cbec_dem_name="cbec As Built"
elif version=='existing':
    cbec_dem_fn='../data/cbec/to_RCD/existing_grade/merged.tif'
    cbec_dem_name="cbec Existing Grade"

cbec_dem = field.GdalGrid(cbec_dem_fn,geo_bounds=clip)
cbec_full_dem=field.GdalGrid(cbec_dem_fn)
cbec_full_dem.name=cbec_dem.name=cbec_dem_name

# 2017 CoNED LIDAR
lidar2017=field.GdalGrid('../data/noaa/lidar/2017CoNED-north_marsh-lidar/Job558295_cent_ca_coned.tif',geo_bounds=clip)
lidar2017.name="CoNED 2017 LiDaR"
lidar2011=field.GdalGrid('../data/noaa/lidar/2011CoastCons-north_marsh-lidar/Job558340_CA2009_coastal_DEM.tif',geo_bounds=clip)
lidar2011.name="Coastal Cons. 2009-2011 LiDaR"

##
if 0:
    plt.figure(1).clf()
    fig,axs=plt.subplots(2,2,num=1)
    fig.set_size_inches([7.5,8.],forward=True)
    raise Exception("Naming for cbec layers factored out. This needs updating to render figure")

    for ax,dem in zip(axs.ravel(),[as_built_dem,existing_dem,lidar2017,lidar2011]):
        img=dem.plot(ax=ax,cmap=sst,vmin=0.5,vmax=3.5)
        ax.axis('tight')
        ax.axis('equal')
        ax.axis(zoom)
        ax.axis('off')
        ax.text(0.95,0.9,dem.name,ha='right',transform=ax.transAxes,
                fontsize=10)

    cax=fig.add_axes([0.02,0.60,0.03,0.25])
    plt.colorbar(img,cax=cax, label="Elevation (m NAVD88)")

    fig.subplots_adjust(left=0.12,right=0.99,wspace=0.01,bottom=0.02,hspace=0.01,top=0.99)

    fig.savefig(os.path.join(fig_dir,"compare-4-lidar-north_marsh_pond.png"),dpi=200)

##

# Show difference between cbec as-built and the 2011 lidar.
cbec_on_2011 = cbec_dem.extract_tile(match=lidar2011)
cbec_on_2011.name=cbec_dem.name

delta=field.SimpleGrid(F=cbec_on_2011.F - lidar2011.F,
                       extents=cbec_on_2011.extents)
delta.name=f"[cbec {version}] - [2011 Lidar]"
delta.smooth_by_convolution(iterations=4)

if 0:
    plt.figure(2).clf()
    fig,axs=plt.subplots(1,3,num=2)
    fig.set_size_inches([8.5,4.5],forward=True)

    for ax,dem in zip(axs.ravel(),[cbec_on_2011,delta,lidar2011]):
        img=dem.plot(ax=ax,cmap=sst,vmin=0.5,vmax=3.5)
        ax.axis('tight')
        ax.axis('equal')
        ax.axis(zoom)
        ax.axis('off')
        ax.set_title(dem.name)
        if dem==delta:
            img.set_clim([-0.5,0.5])
            img.set_cmap('seismic')
            delta_img=img

    cax=fig.add_axes([0.3,0.12,0.4,0.03])
    plt.colorbar(delta_img,cax=cax, label="Difference (m)",orientation='horizontal')
    fig.subplots_adjust(left=0.02,right=0.99,wspace=0.01,bottom=0.10,hspace=0.01,top=0.90)

    fig.savefig(os.path.join(fig_dir,"cbec_vs_lidar2011-north_marsh_pond.png"),dpi=200)


## And identify peak of histogram:
if 0:
    plt.figure(3).clf()
    valid=np.isfinite(delta.F)
    plt.hist( delta.F[valid], bins=np.linspace(-0.2,0.1,500))
    plt.grid()
    plt.savefig(os.path.join(fig_dir,"cbec_vs_lidar2011-histogram.png"),dpi=200)

##

# Bring in the point data from the later sample:
all_points=wkb2shp.shp2geom('../data/cbec/cbec-survey/All_Points.shp',
                            target_srs='EPSG:26910')

xyz=np.array( [ np.array(p) for p in all_points['geom']] )
xyz[:,2]=all_points['Z_m']

##
if 0:
    plt.figure(4).clf()
    fig,axs=plt.subplots(1,3,num=4)
    fig.set_size_inches([8.5,5.],forward=True)

    for ax,dem in zip(axs.ravel(),[cbec_dem,lidar2017,lidar2011]):
        img=dem.plot(ax=ax,cmap='gray',vmin=-1,vmax=5.)
        dem_at_xyz=dem( xyz[:,:2] )
        delta=xyz[:,2] - dem_at_xyz

        scat=ax.scatter(xyz[:,0],xyz[:,1],20,delta,cmap='seismic')
        scat.set_clim([-1,1])

        ax.axis('tight')
        ax.axis('equal')
        ax.axis(zoom)
        ax.axis('off')
        ax.set_title(dem.name)

    cax=fig.add_axes([0.3,0.12,0.4,0.03])
    plt.colorbar(scat,cax=cax,label="[cbec AllPoints] - [DEM]",orientation='horizontal')

    fig.subplots_adjust(left=0.02,right=0.99,wspace=0.01,bottom=0.13,hspace=0.01,top=0.90)

    fig.savefig(os.path.join(fig_dir,"NorthMarsh-3dems-error.png"))

##

# Break up the comparison to AllPoints by polygon:
adj_regions=wkb2shp.shp2geom("dem-polygons.shp")

## 

# east_marsh:
# no inundated areas to worry about.
# 2011 and 2017 are the same, and cbec is offset down 0.073m
# Use the 2017 data: most recent, and no round-trip through TIN
from shapely import geometry
region_name='east_marsh'
region=adj_regions['geom'][ adj_regions['name']==region_name ][0]

survey_in_region=np.array( [ region.contains(geometry.Point(p[:2])) for p in xyz] )

region_xyz=xyz[ survey_in_region ]

region_2017 = lidar2017(region_xyz[:,:2])

error=region_2017 - region_xyz[:,2]

print("Region: %s"%region_name)

print("  mean error:   %.3f m"%np.mean(error)) # 0.197 m
print("  median error: %.3f m"%np.median(error)) # 0.197 m

##

# Center Marsh, dry areas.
# 2017 seems to put this area lower than 2011 or cbec,
# though 2011 and cbec include some additional small
# pockets of inundation.

for region_name in [
        'east_marsh','center_marsh',
        'center_marsh_inundated_east',
        'center_marsh_inundated_west',
        'west_marsh_dry'
]:
    region=adj_regions['geom'][ adj_regions['name']==region_name ][0]

    survey_in_region=np.array( [ region.contains(geometry.Point(p[:2])) for p in xyz] )
    region_xyz=xyz[ survey_in_region ]

    print("---- Region: %s ----"%region_name)
    for dem in [lidar2017,lidar2011]:
        from_dem=dem(region_xyz[:,:2])
        error=from_dem - region_xyz[:,2]
        corr=np.corrcoef(from_dem,region_xyz[:,2])[0,1]

        print(dem.name)
        print("  mean error:   %.3f m"%np.mean(error)) # 0.197 m
        print("  median error: %.3f m"%np.median(error)) # 0.197 m
        print("  R^2:          %.3f"%corr**2) 
        print()

##

# At this point we render an intermediate DEM, ignoring the ponded areas.
from stompy.spatial import field

# Compile corrections so far, starting with the same shapefile
# Reload so we can iterate with QGIS
shp_data=wkb2shp.shp2geom("dem-polygons.shp")

new_fields=[ ('src_name',shp_data['name']),
             ('src',np.zeros(len(shp_data),'O')),
             ('priority',np.zeros(len(shp_data))),
             ('data_mode',np.zeros(len(shp_data),'S200')),
             ('alpha_mode',np.zeros(len(shp_data),'S200')) ]
shp_data=utils.recarray_add_fields(shp_data, new_fields)

shp_data['data_mode']=''
shp_data['alpha_mode']=''
shp_data['priority']=-1 # default is disabled


def params(name,delete=False,**kw):
    global shp_data
    idx=np.nonzero( shp_data['src_name']==name )[0]

    if delete:
        assert len(idx)==1
        sel=np.arange(len(shp_data))!=idx[0]
        shp_data=shp_data[sel]
        return shp_data
    if len(idx)==1:
        idx=idx[0]
    elif len(idx)==0:
        shp_data=utils.array_append(shp_data)
        idx=-1
        shp_data['name'][idx]=shp_data['src_name'][idx]=name
    else:
        raise Exception('Name matches multiple?')
    for k in kw:
        shp_data[k][idx]=kw[k]
    return shp_data[idx]

def factory(feat):
    return params(feat['name'])['src']


def field_bounds(f):
    # This had been hardwired to cbec_dem
    return geometry.box( *[f.extents[i] for i in [0,2,1,3] ] )
    
params('cbec_dem',
       geom=field_bounds(cbec_dem),
       data_mode='fill(3.0)', # existing grade has a few missing pixels.
       priority=20,
       src=cbec_dem)

params('center_marsh',src=lidar2011 - 0.20, priority=40,alpha_mode='feather_out(5.0)')

west_inundated_src=lidar2011 - 0.32 # temporary, while we estimate it below
north_pond_src=lidar2017 # will be updated

params('east_marsh',src=lidar2011 - 0.20,priority=50,alpha_mode='feather_out(5.0)')
params('center_marsh_inundated_east',src=lidar2011 - 0.11,priority=50,
       alpha_mode='feather_out(5.0)')
# Maybe this is the one causing the seam in the west inundation
# was -0.06. But with 0.03, it causes a little ring at the
# boundary.
# feather in to reduce ring
params('center_marsh_inundated_west',src=lidar2011 - 0.06,priority=50,alpha_mode='feather_in(5.0)')

# west_marsh_dryand west_marsh_inundated are combined below
params('west_marsh_dry',priority=-50)
# src=lidar2011 - 0.32,alpha_mode='buffer(-5.0),feather_in(5.0)'

# These will both get updated below 
params('west_marsh_inundated',src=lidar2011-0.32,priority=-50,
       alpha_mode='feather_out(5.0)')
params('north_pond',src=lidar2017,priority=50,
       alpha_mode='feather_out(5.0)')


# Merge west_marsh_dry and west_marsh_inundated to avoid creating
# extra seams.  Better to just blend layers
west_marsh=params('west_marsh_dry')['geom'].union(params('west_marsh_inundated')['geom'])
params('west_marsh',geom=west_marsh,
       # alpha_mode='feather_in(4.0)',
       alpha_mode='blur_alpha(4.0)',
       priority=60,src=lidar2011 - 0.32)

# Marsh--pond channel
# (552227.2488704546, 552470.7314864796, 4124655.0209947866, 4124934.6201988547)
# Seems that west_marsh pushes that low constant area out too far?  then west_marsh_inundated
# is probably pulling too much from it.
# N channel interp was a bit wide => Update footprint in cbec-survey-interp-grid9.pkl
# west_marsh_inundated is buffering and feathering out a lot => Adjust that boundary
#   in north_marsh_pond_adjustment_polygons.shp, and slight adjustments to parameters.

north_channel_dem=field.GdalGrid("north_channel_and_ditch-1m.tif")
params('north_channel',
       priority=75,src=north_channel_dem,
       alpha_mode='valid(), feather_in(2.0)',
       data_mode='min()',
       geom=field_bounds(north_channel_dem))

# Lagoon mouth (551996.7690783864, 552359.0075459109, 4124370.5881716525, 4124786.5586785264)
#   As-built seems to TIN-interpolate over the channel here.
#   Are any of the other DEMs good?  Best is maybe the lidar2017 data.
#   This area will have to get carved out regardless, but that's a later step.
# This is a big improvement, but this is going to evolve as Dane's QCM is incorporated,
# and according to any additional survey data.  As it stands, probably reasonable with
# respect to the cbec survey points and the berm survey points.
params('lagoon-lidar2017',alpha_mode='blur_alpha(10.0)',
       data_mode='overlay()',priority=55,
       src=lidar2017)

# Lagoon  (552129.3226207336, 552521.6926489467, 4124153.228958399, 4124603.800540797)
# The main issue was actually lagoon-lidar2017 above, which had too sharp of a transition
# Also decreased the feather on lagoon, and use min()
lagoon_dem=field.GdalGrid("lagoon-1m.tif")

params('lagoon',
       priority=85,src=lagoon_dem,
       data_mode='min()',
       alpha_mode='valid(),feather_in(10.0)',
       geom=field_bounds(lagoon_dem))

params('feeder_culvert',
       priority=95,src=field.ConstantField(0.6),
       alpha_mode='feather_in(3.0),feather_out(4.0)',
       data_mode='min()')

params('south_channel_culvert',
       priority=95,src=field.ConstantField(1.0),
       alpha_mode='feather_in(3.0)',
       data_mode='min()')

# a few survey points are 1.0.  With the blur_alpha(), as the polygon
# gets narrower the min depth is going to decrease, too.  The value
# of the ConstantField doesn't necessarily dictate the min depth.
params('north_ditch',
       priority=95,src=field.ConstantField(0.7),
       alpha_mode='blur_alpha(4.0)',
       data_mode='min()')


south_channel_dem=field.GdalGrid("north_channel_and_ditch-1m.tif")
# And now include that channel

params('marsh_south_channel',
       priority=90,src=south_channel_dem,
       data_mode='min()',
       alpha_mode='valid(),feather_in(2.0)',
       geom=field_bounds(lagoon_dem))

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0

# This can then be used below with the ContourInterpolator
dem_dry=comp_field.to_grid(dx=res,dy=res,bounds=clip)
# 
# ##
# # Tuning north_ditch
# 
# zoom=(552314., 552519., 4124655., 4124902.) 
# 
# dem_local,stack=comp_field.to_grid(dx=1,dy=1,bounds=zoom,
#                                    stackup='return')
# fig=comp_field.plot_stackup(dem_local, stack,cmap=turbo,num=1,z_factor=1.5)
# fig.tight_layout()
# 
# dem_test=comp_field.to_grid(dx=res,dy=res,bounds=zoom)
# plt.figure(1).clf()
# dem_test.plot(cmap=turbo,clim=[0.7,1.5])

##

# Assume the western pan is sort of radially similar.
from stompy.grid import unstructured_grid, exact_delaunay
from stompy.spatial import interp_concentric
        
west=adj_regions['geom'][ adj_regions['name']=='west_marsh_inundated' ][0]
center=adj_regions['geom'][ adj_regions['name']=='center_marsh_inundated_west' ][0]

# Go out about 10.0m to get some dry land to pull from the DEM 
region = west.union(center).buffer(10.0)
survey_in_region=np.array( [ region.contains(geometry.Point(p[:2])) for p in xyz] )
region_xyz=xyz[ survey_in_region ]

pnt_samples=pd.DataFrame()
pnt_samples['x']=region_xyz[:,0]
pnt_samples['y']=region_xyz[:,1]
pnt_samples['value']=region_xyz[:,2]
pnt_samples['weight']=1.0

# The data above is actually pretty good in this subset of the
# inundated pan (from the 2011 lidar)
good_dem=dem_dry
good_poly=region.intersection( adj_regions['geom'][ adj_regions['name']=='center_marsh_inundated_west'][0] )
good_mask=good_dem.polygon_mask(good_poly)

X,Y = good_dem.XY()

dem_samples=pd.DataFrame()
dem_samples['x']=X[good_mask]
dem_samples['y']=Y[good_mask]
dem_samples['value']=good_dem.F[good_mask]
dem_samples['weight']=0.1

samples=pd.concat([pnt_samples,dem_samples])

# For boundary cells additionally pull from DEM

ci=interp_concentric.ConcentricInterpolator(region,samples,dx=3.,
                                            anisotropy=0.02,
                                            background_weight=0.01,
                                            alpha=0.05,
                                            background_field=dem_dry)

# And convert to raster:
dem_west_inundated=field.rasterize_grid_cells(ci.grid,values=ci.result,dx=2,dy=2,
                                              stretch=True)
dem_west_inundated.smooth_by_convolution(iterations=1)
if 0:
    plt.figure(1).clf()
    dem_west_inundated.plot(cmap=turbo,vmin=0.5,vmax=2.50)
    dem_west_inundated.plot_hillshade()

# the -0.05 here gets a nice transition with the existing
# good pan data.  Might be more 
params('west_marsh_inundated',src=dem_west_inundated-0.05,
       priority=70,alpha_mode='valid(),buffer(5.0),feather_out(10.0)')


# Testing out blend with inundation
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0

# Testing out blend with inundation
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)


# North Pond
# Go out about 10.0m to get some dry land to pull from the DEM 
region = adj_regions['geom'][ adj_regions['name']=='north_pond' ][0]
survey_in_region=np.array( [ region.contains(geometry.Point(p[:2])) for p in xyz] )
region_xyz=xyz[ survey_in_region ]

pnt_samples=pd.DataFrame()
pnt_samples['x']=region_xyz[:,0]
pnt_samples['y']=region_xyz[:,1]
pnt_samples['value']=region_xyz[:,2]
pnt_samples['weight']=1.0

dem_so_far=comp_field.to_grid(dx=2,dy=2,
                              bounds=[region.bounds[i] for i in [0,2,1,3]])

ci=interp_concentric.ConcentricInterpolator(region,pnt_samples,dx=3.,
                                            anisotropy=0.05,
                                            alpha=0.5,
                                            background_field=dem_so_far,
                                            background_weight=0.001)

fig,ax,items=ci.plot_result(clim=[0.5,2.5])
items[1].set_lw(0.5) ; items[1].set_edgecolor('k')

# And convert to raster:
dem_north_pond=ci.rasterize(dx=2)
dem_north_pond.smooth_by_convolution()

params('north_pond',src=dem_north_pond,
       priority=80,data_mode='min()',
       alpha_mode='feather_in(5.0)')

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0
# This can then be used below with the ContourInterpolator
dem_wetdry=comp_field.to_grid(dx=res,dy=res,bounds=clip)

# --
plt.figure(1).clf()
fig,ax=plt.subplots(1,1,num=1)
fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
ax.axis('off')
ax.axis('tight')
ax.axis('equal')
img=dem_wetdry.plot(cmap=turbo,vmin=0,vmax=3.5)

dem_wetdry.plot_hillshade(ax=ax)
plt.colorbar(img)

# And now include that channel
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)

# This can then be used below with the ContourInterpolator
dem_wetdry=comp_field.to_grid(dx=1,dy=1,bounds=clip)

## 
# Debug and adjust small areas here.

# N. Marsh internal levee:
# (552383.5041504747, 552529.5834280611, 4124606.5168135753, 4124743.4661363126)
# west_marsh looks decent here, but then west_marsh_inundated is blurred out
# and obscures probably too much of west_marsh
#   it is set to valid(),buffer(5.0),feather_out(10.0)

params('west_marsh',
       # alpha_mode='blur_alpha(4.0)'
       data_mode='src_data.F[src_data.F<1.60]=np.nan ; fill(10.0)',
       alpha_mode='valid(),blur_alpha(4.0)'
)

params('west_marsh_inundated',
       #alpha_mode='valid(),buffer(5.0),feather_out(10.0)'
       alpha_mode='valid(),buffer(2.0),feather_in(2.0)'
       )

params('n_marsh_pan_connector',
       # 2021-09-13: Make this 1ft shallower than it had been
       # 2021-09-20: And even shallower.
       src=field.ConstantField(1.9),
       priority=99,
       data_mode='min()',
       # And a bit smoother on the burn in 
       alpha_mode='buffer(1.0), feather_out(2.0)')

# N Marsh fill
six.moves.reload_module(field)
                        
params('west_marsh_fill',
       src=field.ConstantField(10),
       priority=105,
       data_mode='diffuser(),min()',
       alpha_mode='feather_in(2)')

##

# And do something a bit better for the mouth -- pull from
# the UCB grid.
pdo_grid=unstructured_grid.RgfGrid("../data/ucb-model/MouthVariation/EpMp/PDO_EpMp.grd")
pdo_grid.add_node_field('z_bed_node',-pdo_grid.nodes['depth_node'])

# Convert that to a Field
pdo_field=field.XYZField(X=pdo_grid.nodes['x'],F=pdo_grid.nodes['z_bed_node'])
pdo_field._tri=pdo_grid.mpl_triangulation()

##
params('mouth',
       src=pdo_field,
       priority=100,
       data_mode='blur(5)',
       alpha_mode='feather_in(10.0),blur_alpha(4)')

if dredge:
    params('lag_thalweg',
           src=field.ConstantField(0.0),
           priority=110,
           data_mode='min()',
           alpha_mode='feather_in(5.0),blur_alpha(5.0)')

# The lidar gets the sediment plugs reasonably well, so nothing
# to do but keep the north marsh channel from dredging too much
# here.

##

# See if USGS topo is better than PDO

usgs_dem=field.GdalGrid("../data/usgs/CentCA_Topobathy_DEM_1m_auLleA1aqXot7vR4sWse.tiff")

# Seems to be largely similar to lidar2017. How about just overlay in the same areas?

params('usgs_dem',alpha_mode='blur_alpha(10.0)',
       data_mode='overlay()',priority=57, # just above lagoon lidar
       src=usgs_dem)

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)

##

# Add in the cut off Pescadero onto Delta Marsh.
# The layout here is digitized roughly from waterways
# surface, but the match from that surface to cbec
# is not very good, so I'm not pulling the surface
# directly, and rather just prescribing some values.

params('pescadero_cut',
       src=field.ConstantField(2.75),
       priority=100,
       data_mode='min()',
       alpha_mode='feather(4.0)')

if 0:
    params('cbec_dem',
           geom=geometry.box( *[cbec_full_dem.extents[i] for i in [0,2,1,3] ] ),
           src=cbec_full_dem)

    comp_field=field.CompositeField(shp_data=shp_data,
                                    factory=factory)

    dem_local,stack=comp_field.to_grid(dx=1,dy=1,bounds=(553220., 553407.,
                                                         4123793., 4123912.),
                                       stackup='return')
    fig=comp_field.plot_stackup(dem_local, stack,cmap=turbo,num=3,z_factor=1.5)
    fig.tight_layout()

##

ww_pesc_dem=field.GdalGrid("../data/waterways/Pescadero-Model-Terrain/"
                           "Existing_Conditions_HECRAS_Surface_2020-10-27/"
                           "18-057-EG-V5-20201027_EPSG2610_m.tif")
# The generous feather is to overwrite the polygonal channel from the cbec surface
params('waterways_pesc_channel',
       priority=100,
       src=ww_pesc_dem,
       data_mode='overlay()',
       alpha_mode='feather_out(10.0)')

## 

params('cbec_dem',
       geom=geometry.box( *[cbec_full_dem.extents[i] for i in [0,2,1,3] ] ),
       src=cbec_full_dem)

blend_poly=params('pesc_blend_fiction')['geom']
params('pesc_blend_fiction',priority=-1)
xyxy=blend_poly.bounds
comp_field=field.CompositeField(shp_data=shp_data,factory=factory)
blend_src=comp_field.to_grid(dx=1,dy=1,bounds=[xyxy[0],xyxy[2],xyxy[1],xyxy[3]])

##
blend_int=interp_orthogonal.OrthoInterpolator(blend_poly,nom_res=2.5,background_field=blend_src,
                                              anisotropy=1e-5)
blend_fld=blend_int.field()

params('pesc_blend_fiction',priority=105,src=blend_fld,
       data_mode='min()',alpha_mode='valid(),feather(1.0)')

##


if 1:
    # For the spot fixes below, be sure I'm using the cropped as_built to keep
    # things speedy.
    params('cbec_dem',
           geom=geometry.box( *[cbec_full_dem.extents[i] for i in [0,2,1,3] ] ),
           src=cbec_full_dem)

    comp_field=field.CompositeField(shp_data=shp_data,
                                    factory=factory)

    dem_local,stack=comp_field.to_grid(dx=1,dy=1,bounds=( 551900., 552250.,
                                                         4124450, 4124850 ),
                                       stackup='return')
    fig=comp_field.plot_stackup(dem_local, stack,cmap=turbo,num=3,z_factor=1.5)
    fig.tight_layout()
    plt.setp( fig.axes, adjustable='datalim')
    
    print()
    print(f"{'pri':4} {'src':30}  {'data_mode':12} {'alpha_mode':12}")
    print("---------------------------------------------------------------------")

    for idx in np.argsort(shp_data['priority']):
        row=shp_data[idx]
        print(f"{row['priority']:4.0f} {row['src_name']:30}  {row['data_mode'].decode():12} {row['alpha_mode'].decode():12}")

    # for ax in plt.gcf().axes:
    #     plot_wkb.plot_wkb(blend_poly,ax=ax,fc='none')

##

import stompy.plot.cmap as scmap
sst=scmap.load_gradient('oc-sst.cpt')

if 1: # Final render
    # bring in the full as_built DEM
    params('cbec_dem',
           geom=geometry.box( *[cbec_full_dem.extents[i] for i in [0,2,1,3] ] ),
           src=cbec_full_dem)

    comp_field=field.CompositeField(shp_data=shp_data,
                                    factory=factory)

    # Render a tile to match as built
    dem_final=comp_field.to_grid(dx=1,dy=1,bounds=cbec_full_dem.extents)

    plt.figure(1).clf()
    fig,ax=plt.subplots(1,1,num=1)
    fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
    ax.axis('off')
    ax.axis('tight')
    ax.axis('equal')

    # trim out high land that's not useful for the model
    dem_final.F[ dem_final.F> 25.0]=np.nan
    
    img=dem_final.plot(cmap=sst,vmin=0,vmax=3.5)

    dem_final.plot_hillshade(ax=ax,z_factor=3)
    plt.colorbar(img)

    dem_final.write_gdal(f'compiled-dem-{version}-{date_str}-1m.tif',overwrite=True)

##

# Write a shapefile of the composite regions
wkb2shp.wkb2shp(f"composite-input-{version}-{date_str}.shp",shp_data['geom'],
                fields=shp_data,overwrite=True)

##

if render_05m:
    # Final render at 0.5m
    comp_field=field.CompositeField(shp_data=shp_data,
                                    factory=factory)

    # use a post processing step to separately save rgb versions of the tiles
    rgb_tiles=[]
    def post(tile,**kw):
        tile_rgb=tile.to_rgba(cmap=sst,vmin=0,vmax=5.0)
        shade_rgba=tile.hillshade_shader(z_factor=3)
        # Manually composite:
        tile_rgb.overlay_rgba(shade_rgba)

        dem_fn=kw['output_fn']
        rgb_fn=dem_fn.replace('.tif','-rgb.tif')
        tile_rgb=tile_rgb.crop(kw['bounds'])
        tile_rgb.assign_projection('EPSG:26910')

        os.path.exists(rgb_fn) and os.unlink(rgb_fn)
        tile_rgb.write_gdal_rgb(rgb_fn)
        rgb_tiles.append(rgb_fn)
        return tile

    tm=field.TileMaker(comp_field,dx=0.5,dy=0.5,tx=500,ty=500,
                       output_dir='rendered-0.5m')
    tm.post_render=post
    tm.force=True
    extents=cbec_full_dem.extents
    tm.tile(extents[0],extents[2],extents[1],extents[3])
    tm.merge()
    ##

    extents = cbec_full_dem.extents

    # For more general use, 1m is fine.
    tm=field.TileMaker(comp_field,dx=1.0,dy=1.0,tx=1000,ty=1000,
                       output_dir='rendered-1.0m')
    tm.force=True
    tm.tile(extents[0],extents[2],extents[1],extents[3])
    tm.merge()

    ##

    # And merge the rgbs:
    merged_rgb_fn='rendered-0.5m/merged-rgb.tif'
    cmd=f"gdal_merge.py -o {merged_rgb_fn} " + " ".join(rgb_tiles)
    import subprocess
    subprocess.call(cmd,shell=True)

    ##
    #  Render with oc-sst and hillshade to an RGB geotiff.
    #  This is super slow, and not scalable.  Would be better to
    #  (a) fix TileMaker to include pads, and
    #  (b) provide a per-tile post-processing step in TileMaker.

    if 1:
        demc=field.GdalGrid('rendered-0.5m/merged.tif',geo_bounds=clip)

        # Generate a colorbar
        fig=plt.figure(1)
        fig.clf()
        fig.set_size_inches([6.4, 4.8],forward=True)
        img=demc.plot(cmap=sst,vmin=0,vmax=5.0)
        plt.colorbar(img,label='m NAVD88')
        fig.axes[0].set_visible(0)
        fig.savefig('topo-colorbar.png',bbox_inches='tight')
    ## 
    #  Warp to 4326 and split to KML tiles
    merc_rgb_fn=merged_rgb_fn.replace('.tif','-EPSG4326.tif')
    assert merc_rgb_fn!=merged_rgb_fn

    subprocess.call(f'gdalwarp -s_srs EPSG:26910 -t_srs EPSG:4326 {merged_rgb_fn} {merc_rgb_fn}',
                    shell=True)

    tile_dir=f'pescadero-{version}-{date_str}'
    subprocess.call(f'gdal2tiles.py -k -z 12-18 {merc_rgb_fn} {tile_dir}',
                    shell=True)

    ##

    # Add the colorbar into the kml
    colorbar_kml="""
     <ScreenOverlay>
         <name>
             Legend: Topography
         </name>
         <Icon>
           <href>topo-colorbar.png</href>
         </Icon>
         <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
         <screenXY x="25" y="95" xunits="pixels" yunits="pixels"/>
         <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
         <size x="0" y="0" xunits="pixels" yunits="pixels"/>
     </ScreenOverlay>
    """

    shutil.copyfile('topo-colorbar.png',
                    os.path.join(tile_dir,'topo-colorbar.png'))

    with open(os.path.join(tile_dir,'doc.kml')) as fp:
        orig_kml=fp.read()

    with open(os.path.join(tile_dir,f'{version}_updated.kml'),'wt') as fp:
        new_kml=orig_kml.replace('</Document>',
                                 colorbar_kml+"\n</Document>")
        fp.write(new_kml)

    ##

    dem=field.GdalGrid('rendered-1.0m/merged.tif')

    fig=plt.figure(1)
    fig.clf()
    fig.set_size_inches((7,7),forward=True)

    ax=fig.add_axes([0,0,1,1])

    img=dem.plot(cmap=sst,vmin=0,vmax=5.0,ax=ax)
    ax.axis('tight')
    ax.axis('equal')
    ax.axis( (551954, 553326., 4124051., 4125441))

    cax=fig.add_axes([0.7,0.65,0.03,0.25])

    plt.colorbar(img,cax=cax,label='m NAVD88')
    fig.savefig('overview-lagoon-marsh.png')
