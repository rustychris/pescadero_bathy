from stompy.spatial import field
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap

##


## 

dem=field.GdalGrid('compiled-dem-existing-20210608-1m.tif')

##
#cmap=scmap.load_gradient('hot_desaturated.cpt')
cmap=scmap.load_gradient('BlueYellowRed.cpt')

fig=plt.figure(num=1)
fig.set_size_inches((6,6.5),forward=True)
fig.clf()

ax=fig.add_axes([0,0,1,1])

ax.axis('off')
zoom=(552000, 554343., 4.1227e6, 4.12545e6)
dc=dem.crop(zoom)
img=dc.plot(ax=ax,cmap=cmap,vmin=-0.5,vmax=4.5)
dc.plot_hillshade(z_factor=4.)
cax=fig.add_axes([0.7,0.68,0.03,0.3])
plt.colorbar(img,cax=cax,label="Elevation (m)")
ax.axis('tight')
ax.axis('equal')

fig.savefig('bathy-overview.png',dpi=120)
