import six
import numpy as np
import matplotlib.pyplot as plt
from stompy.spatial import generate_dem, field
from stompy import utils
six.moves.reload_module(field)
six.moves.reload_module(generate_dem)

##

generate_dem.main(["-s","dem-scenarios.shp",
                   "-o","scen1",
                   "-r","1",
                   "-f",
                   "-p",".",
                   "-b","552250","553320","4124040","4124900",
                   "-q","tags like '%nm1%'"])
## 
generate_dem.main(["-s","dem-scenarios.shp",
                   "-o","scen2",
                   "-r","1",
                   "-f",
                   "-p",".",
                   "-b","552250","553320","4124040","4124900",
                   "-q","tags like '%nm2%'"])
## 
generate_dem.main(["-s","dem-scenarios.shp",
                   "-o","scen3",
                   "-r","1",
                   "-f",
                   "-p",".",
                   "-b","552250","553320","4124040","4124900",
                   "-q","tags like '%nm3%'"])
## 

# Plots of the difference
base=field.GdalGrid("compiled-dem-asbuilt-20210920-1m.tif")
scen1=field.GdalGrid("scen1/output_res1.tif")
scen2=field.GdalGrid("scen2/output_res1.tif")
scen3=field.GdalGrid("scen3/output_res1.tif")

##

for s_num,scen in enumerate([scen1, scen2, scen3]):
    delta=base.extract_tile(match=scen)
    delta.F=scen.F-delta.F

    changed=np.isfinite(delta.F) & (delta.F!=0.0)
    if not np.any(changed):
        print(f"No changes for scenario {1+s_num}")
        continue

    fignum=1+s_num
    plt.figure(fignum).clf()
    fig,axs=plt.subplots(1,2,num=fignum)
    fig.set_size_inches([10.5,5],forward=True)
    axs[0].set_adjustable('datalim')
    axs[1].set_adjustable('datalim')

    img1=scen.plot(ax=axs[0],cmap='turbo')
    img1hs=scen.plot_hillshade(ax=axs[0],z_factor=2)
    img1.set_clim([-0.5,4.5])
    
    img2=delta.plot(ax=axs[1],cmap='coolwarm')
    img2.set_clim([-1,1])

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
        zoom=utils.expand_xxyy([xmin,xmax,ymin,ymax],0.15)

    axs[0].axis(zoom)
    axs[1].axis(zoom)
    fig.tight_layout()
    fig.savefig(f'scen{s_num+1}-delta-zoomed.png',dpi=200)
##


# These have been inserted into the initial report document.
