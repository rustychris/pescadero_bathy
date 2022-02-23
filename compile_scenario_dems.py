import six
from stompy.spatial import generate_dem, field
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

# Plots of the difference
base=field.GdalGrid("compiled-dem-asbuilt-20210920-1m.tif")
scen1=field.GdalGrid("scen1/output_res1.tif")
scen2=field.GdalGrid("scen2/output_res1.tif")

