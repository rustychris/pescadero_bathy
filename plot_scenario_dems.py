import six
import numpy as np
import matplotlib.pyplot as plt
from stompy.spatial import generate_dem, field
from stompy import utils

# Plots of the difference
base=field.GdalGrid("compiled-dem-asbuilt-20210920-1m.tif")
scen1=field.GdalGrid("scen1/output_res1.tif")
scen2=field.GdalGrid("scen2/output_res1.tif")
scen3=field.GdalGrid("scen3/output_res1.tif")

##

for s_num,scen in enumerate([scen1, scen2, scen3]):
    if s_num<2: continue # DBG
    delta=base.extract_tile(match=scen)
    delta.F=scen.F-delta.F

    changed=np.isfinite(delta.F) & (delta.F!=0.0)
    if not np.any(changed):
        print(f"No changes for scenario {1+s_num}")
        continue

    fignum=1+s_num
    plt.figure(fignum).clf()
    fig,axs=plt.subplots(1,2,num=fignum)
    if scen==scen2:
        fig.set_size_inches([10.5,5],forward=True)
    else: # scen==scen3
        fig.set_size_inches([10.5,3.75],forward=True)
    axs[0].set_adjustable('datalim')
    axs[1].set_adjustable('datalim')

    img1=scen.plot(ax=axs[0],cmap='turbo')
    img1hs=scen.plot_hillshade(ax=axs[0],z_factor=2)
    img1.set_clim([-0.5,4.5])
    
    img2=delta.plot(ax=axs[1],cmap='coolwarm')
    img2.set_clim([-1,1])
    img2hs=scen.plot_hillshade(ax=axs[1],z_factor=2)

    plt.colorbar(img2,#cax=fig.add_axes([0.9,0.2,0.015,0.35]),
                 ax=axs[1],
                 label="Change in elevation (m)")
    plt.colorbar(img1,
                 ax=axs[0],
                 label="Elevation (m)")

    axs[0].xaxis.set_visible(0)
    axs[0].yaxis.set_visible(0)
    axs[1].xaxis.set_visible(0)
    axs[1].yaxis.set_visible(0)


    if 0:
        zoom=(553016.3736891494, 553276.0100527861, 4124131.5179354837, 4124358.523870968)
    else:
        # Choose zoom based on where delta is nonzero
        x,y=delta.xy()
        if not np.any(changed):
            print(f"No changes for scenario {1+s_num}")
            continue
        x_changed=np.nonzero( np.any(changed,axis=0) )[0]
        y_changed=np.nonzero( np.any(changed,axis=1) )[0]
        xmin=x[x_changed[0]]
        xmax=x[x_changed[-1]]
        ymin=y[y_changed[0]]
        ymax=y[y_changed[-1]]
        zoom=utils.expand_xxyy([xmin,xmax,ymin,ymax],0.05)
        if zoom[0]<delta.extents[0]:
            zoom[0]=delta.extents[0]
        if zoom[1]>delta.extents[1]:
            zoom[1]=delta.extents[1]
        if zoom[2]<delta.extents[2]:
            zoom[2]=delta.extents[2]
        if zoom[3]>delta.extents[3]:
            zoom[3]=delta.extents[3]

    fig.tight_layout()
    axs[0].axis(zoom)
    axs[1].axis(zoom)
    fig.savefig(f'scen{s_num+1}-delta-zoomed.png',dpi=200)
##


# These have been inserted into the initial report document.
