import os
import shapely, fiona
from shapely.geometry import mapping, shape


# function return geometries as geoJSON
def getGeometries_planet(shpfile):
    
    with fiona.open(shpfile, "r") as shapefile:
        geoms = [feature["geometry"] for feature in shapefile]
        
    return geoms
    
    
shp = r"C:\Projects\rd\planet\proj_cc_sites.shp"
polys = getGeometries_planet(shp)

# convert to centroids
centroids = [shape(p).centroid for p in polys]

# buffer out 100 meters. square capping.
buffer_distance = 100.
buffer_p = [c.buffer(buffer_distance, cap_style=3) for c in centroids]

