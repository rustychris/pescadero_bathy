# Load in Lidar and pre-compiled DEMs, particularly in N Marsh
#
import os
import pandas as pd
from stompy.spatial import field, wkb2shp
import six
import matplotlib.pyplot as plt
import numpy as np
import stompy.plot.cmap as scmap
from stompy.plot import plot_wkb
from stompy import utils

turbo=scmap.load_gradient('turbo.cpt')

turbo_gray=scmap.load_gradient('turbo.cpt')
turbo_gray.set_over( "0.8")

sst=scmap.load_gradient('oc-sst.cpt')
##

fig_dir="figs-20200806"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

##

clip=zoom=(551800, 553290, 4124100, 4125400)

# cbec surfaces:
# For starters load just the N marsh and pond portion.
# For compilation at the end load the whole thing
as_built_dem = field.GdalGrid('../data/cbec/to_RCD/asbuilt/merged.tif',geo_bounds=clip)
as_built_dem.name="cbec As Built"

as_built_extents,res = field.GdalGrid.metadata('../data/cbec/to_RCD/asbuilt/merged.tif')

existing_dem =  field.GdalGrid('../data/cbec/to_RCD/existing_grade/merged.tif',geo_bounds=clip)
existing_dem.name="cbec Existing Grade"

# 2017 CoNED LIDAR
lidar2017=field.GdalGrid('../data/noaa/lidar/2017CoNED-north_marsh-lidar/Job558295_cent_ca_coned.tif',geo_bounds=clip)
lidar2017.name="CoNED 2017 LiDaR"
lidar2011=field.GdalGrid('../data/noaa/lidar/2011CoastCons-north_marsh-lidar/Job558340_CA2009_coastal_DEM.tif',geo_bounds=clip)
lidar2011.name="Coastal Cons. 2009-2011 LiDaR"

##

plt.figure(1).clf()
fig,axs=plt.subplots(2,2,num=1)
fig.set_size_inches([7.5,8.],forward=True)

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
cbec_on_2011 = as_built_dem.extract_tile(match=lidar2011)
cbec_on_2011.name=as_built_dem.name

delta=field.SimpleGrid(F=cbec_on_2011.F - lidar2011.F,
                       extents=cbec_on_2011.extents)
delta.name="[As Built] - [2011 Lidar]"
delta.smooth_by_convolution(iterations=4)

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
plt.figure(4).clf()
fig,axs=plt.subplots(1,3,num=4)
fig.set_size_inches([8.5,5.],forward=True)

for ax,dem in zip(axs.ravel(),[as_built_dem,lidar2017,lidar2011]):
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
adj_regions=wkb2shp.shp2geom("north_marsh_pond_adjustment_polygons.shp")


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
        #'east_marsh','center_marsh',
        #'center_marsh_inundated_east',
        #'center_marsh_inundated_west',
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
shp_data=wkb2shp.shp2geom("north_marsh_pond_adjustment_polygons.shp")

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
    return geometry.box( *[as_built_dem.extents[i] for i in [0,2,1,3] ] )
    
params('as_built',
       geom=geometry.box( *[as_built_extents[i] for i in [0,2,1,3] ] ),
       priority=20,
       src=as_built_dem)

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

# there's a big step, and it's better to just avoid the step.
params('west_marsh_dry',src=lidar2011 - 0.32,priority=-50,
       alpha_mode='buffer(-5.0),feather_in(5.0)')

# will get updated below
params('west_marsh_inundated',src=lidar2011-0.32,priority=-50,
       alpha_mode='feather_out(5.0)') HERE
params('west_marsh_inundated',alpha_mode='valid(),buffer(5.0),feather_out(10.0)') # ORIG


# will get updated below
params('north_pond',src=lidar2017,priority=50,
       alpha_mode='feather_out(5.0)')

## 
# Merge west_marsh_dry and west_marsh_inundated to avoid creating
# extra seams.  Better to just blend layers

west_marsh=params('west_marsh_dry')['geom'].union(params('west_marsh_inundated')['geom'])
params('west_marsh',geom=west_marsh,
       # alpha_mode='feather_in(4.0)',
       alpha_mode='blur_alpha(4.0)'),
       priority=60,src=lidar2011 - 0.32)


# Marsh--pond channel
# (552227.2488704546, 552470.7314864796, 4124655.0209947866, 4124934.6201988547)
# Seems that west_marsh pushes that low constant area out too far?  then west_marsh_inundated
# is probably pulling too much from it.
# N channel interp was a bit wide => Update footprint in cbec-survey-interp-grid9.pkl
# west_marsh_inundated is buffering and feathering out a lot => Adjust that boundary
#   in north_marsh_pond_adjustment_polygons.shp, and slight adjustments to parameters.
params('north_channel',data_mode='min()',alpha_mode='valid(),feather_in(2.0)')

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





## 



# -----------------------
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0

    
##     
# This can then be used below with the ContourInterpolator
dem_dry=comp_field.to_grid(dx=res,dy=res,bounds=clip)

## 
#----
# Assume the western pan is sort of radially similar.
from stompy.grid import unstructured_grid, exact_delaunay
from stompy.spatial import interp_concentric
six.moves.reload_module(unstructured_grid)
six.moves.reload_module(exact_delaunay)
six.moves.reload_module(interp_concentric)
        
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

# fig,ax,items=ci.plot_result(clim=[0.5,2.5])

# And convert to raster:
dem_west_inundated=field.rasterize_grid_cells(ci.grid,values=ci.result,dx=2,dy=2,
                                              stretch=True)
dem_west_inundated.smooth_by_convolution(iterations=1)

plt.figure(1).clf()
dem_west_inundated.plot(cmap=turbo,vmin=0.5,vmax=2.50)
dem_west_inundated.plot_hillshade()

##
# the -0.05 here gets a nice transition with the existing
# good pan data.  Might be more 
params('west_marsh_inundated',src=dem_west_inundated-0.05,
       priority=70,alpha_mode='valid(),buffer(5.0),feather_out(10.0)')

# Testing out blend with inundation
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0

# --
if 1: # For debugging the inundated area, before the concentric interp
    dem_local,stack=comp_field.to_grid(dx=res,dy=res,
                                       bounds=(552300, 552766, 4124530, 4124853),
                                       stackup='return')

    plt.figure(3).clf()
    fig,axs=plt.subplots(3,2,num=3)

    for ax,(name,data,alpha) in zip( axs.ravel(), stack ):
        data.plot(ax=ax,vmin=0,vmax=3.5,cmap=turbo)
        data.plot_hillshade(ax=ax)
        ax.axis('off')
        ax.set_title(name)
    fig.subplots_adjust(left=0,right=1,top=0.95,bottom=0,hspace=0.08)


##

# And now include that channel
north_channel_dem=field.GdalGrid("north_channel-interp-1m.tif")

params('north_channel',
       priority=75,src=north_channel_dem,
       alpha_mode='valid(), feather_in(2.0)',
       data_mode='min()',
       geom=field_bounds(north_channel_dem))

# Testing out blend with inundation
comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)


if 0: 
    dem_wetdry=comp_field.to_grid(dx=1,dy=1,bounds=clip)
    plt.figure(1).clf()
    fig,ax=plt.subplots(1,1,num=1)
    fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
    ax.axis('off')
    ax.axis('tight')
    ax.axis('equal')
    img=dem_wetdry.plot(cmap=turbo,vmin=0,vmax=3.5)

    dem_wetdry.plot_hillshade(ax=ax,z_factor=3)
    plt.colorbar(img)
    ax.axis( (552486.9713107223, 552835.1304622293, 4124335.4499081606, 4124673.0755821546) )

## 
if 0: # For debugging the inundated area, before the concentric interp
    # xxyy=(552302.3362373341, 552394.1411568073, 4124665.499731275, 4124729.251936685)
    dem_local,stack=comp_field.to_grid(dx=2.,dy=2,
                                       bounds=clip,
                                       stackup='return')
##
if 0:
    plt.figure(3).clf()
    fig,axs=plt.subplots(2,2,num=3,sharex=True,sharey=True)

    for ax,(name,data,alpha) in zip( axs.ravel(), stack[::-1] ):
        data.plot(ax=ax,vmin=0,vmax=3.5,cmap=turbo)
        # data.plot_hillshade(ax=ax,z_factor=2)
        ax.axis('off')
        ax.set_title(name)
    fig.subplots_adjust(left=0,right=1,top=0.95,bottom=0,hspace=0.08)
##

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
## 
params('north_pond',src=dem_north_pond,
       priority=80,data_mode='min()',
       alpha_mode='feather_in(5.0)')


##
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
# ax.axis( (552510, 552766, 4124530, 4124853))

##
# p0=np.array([552625., 4124706.])
# p1=np.array([552631., 4124756.])
# 
# transect=p0 + np.linspace(0,1,50)[:,None]*(p1-p0)
# 
# d=utils.dist_along(transect)
# z=dem_wetdry(transect)
# 
# # And the sources:
# pad=10.0
# tran_bounds=( transect[:,0].min()-pad, transect[:,0].max()+pad,
#               transect[:,1].min()-pad, transect[:,1].max()+pad )
# dem_tran,stack=comp_field.to_grid(dx=res,dy=res,bounds=tran_bounds,return_stack=True)


# getting close, but the transition is still a little bump
# Might have hit the limit on just buffering -- need to
# pass a bigger region to InterpConcentric.

# plt.figure(2).clf()
# fig,axs=plt.subplots(3,1,num=2)
# axs[0].plot(d,dem_wetdry(transect),label='result')
# axs[0].plot(d,dem_tran(transect),label='transect')
# 
# for name,data,alpha in stack:
#     axs[1].plot(d, data(transect), label=name)
#     axs[2].plot(d, alpha(transect), lw=3.5,label=name)
# axs[1].legend(loc=[1.05,0])
# axs[2].legend(loc=[1.05,0])
# 
# axs[1].plot( d, srcs['west_marsh_dry'](transect), label='west marsh dry source',
#              lw=3,alpha=0.3)
# 
# fig.subplots_adjust(right=0.7)


##

lagoon_dem=field.GdalGrid("lagoon-interp-1m.tif")
# And now include that channel

params('lagoon',
       priority=85,src=lagoon_dem,
       alpha_mode='valid(),feather_in(20.0)',
       geom=field_bounds(lagoon_dem))

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)
res=2.0
# This can then be used below with the ContourInterpolator

dem_wetdry=comp_field.to_grid(dx=1,dy=1,bounds=clip)


##

south_channel_dem=field.GdalGrid("marsh_south_channel-interp-1m.tif")
# And now include that channel

params('marsh_south_channel',
       priority=90,src=south_channel_dem,
       data_mode='min()',
       alpha_mode='valid(),feather_in(2.0)',
       geom=field_bounds(lagoon_dem))

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)

dem_wetdry=comp_field.to_grid(dx=1,dy=1,bounds=clip)

# Decent cut for the N marsh and pond.
# Missing two channels: one that goes east after the culverts
# along the S side of the marsh, and one that goes along the
# north side of the marsh.

# Paste in the lagoon DEM, generate a full 1m DEM, and
# see where we stand.  That might be the thing to pass on
# for review at this point, and come back to the small
# ditches later.

# --
plt.figure(1).clf()
fig,ax=plt.subplots(1,1,num=1)
fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
ax.axis('off')
ax.axis('tight')
ax.axis('equal')
img=dem_wetdry.plot(cmap=turbo,vmin=0,vmax=3.5)

dem_wetdry.plot_hillshade(ax=ax,z_factor=3)
plt.colorbar(img)
ax.axis( (552486.9713107223, 552835.1304622293, 4124335.4499081606, 4124673.0755821546) )

##

# Debug and adjust small areas here.

if 0:
    # For the spot fixes below, be sure I'm using the cropped as_built to keep
    # things speedy.
    params('as_built',
           geom=geometry.box( *[as_built_dem.extents[i] for i in [0,2,1,3] ] ),
           src=as_built_dem)
## 

comp_field=field.CompositeField(shp_data=shp_data,
                                factory=factory)

dem_local,stack=comp_field.to_grid(dx=1,dy=1,bounds=(552415.9434701396, 552736.8874662898, 4124555.6072339383, 4124924.1579228495),
                                   stackup='return')
fig=comp_field.plot_stackup(dem_local, stack,cmap=turbo,num=3,z_factor=1.5)
fig.tight_layout()

print()
print(f"{'pri':4} {'src':30}  {'data_mode':12} {'alpha_mode':12}")
print("---------------------------------------------------------------------")

for idx in np.argsort(shp_data['priority']):
    row=shp_data[idx]
    print(f"{row['priority']:4.0f} {row['src_name']:30}  {row['data_mode'].decode():12} {row['alpha_mode'].decode():12}")

## 

if 1:
    # For the final render, bring in the full as_built DEM
    full_as_built_dem=field.GdalGrid('../data/cbec/to_RCD/asbuilt/merged.tif')

    params('as_built',
           geom=geometry.box( *[full_as_built_dem.extents[i] for i in [0,2,1,3] ] ),
           src=full_as_built_dem)

    comp_field=field.CompositeField(shp_data=shp_data,
                                    factory=factory)

    extents,res = field.GdalGrid.metadata('../data/cbec/to_RCD/asbuilt/merged.tif')

    # Render a tile to match as built
    dem_final=comp_field.to_grid(dx=1,dy=1,bounds=extents)

    # 
    plt.figure(1).clf()
    fig,ax=plt.subplots(1,1,num=1)
    fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
    ax.axis('off')
    ax.axis('tight')
    ax.axis('equal')
    img=dem_final.plot(cmap=turbo,vmin=0,vmax=3.5)

    dem_final.plot_hillshade(ax=ax,z_factor=3)
    plt.colorbar(img)

    dem_final.write_gdal('compiled-dem-20200813-1m.tif',overwrite=True)
