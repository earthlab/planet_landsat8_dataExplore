import planet
import json
import os
import requests
from requests.auth import HTTPBasicAuth
import pprint

# create a dictionary for planet labs item types
planet_itemtype_dict = {}
planet_itemtype_dict['ps3'] = 'PSScene3Band'	#PlanetScope Scenes
planet_itemtype_dict['ps4'] = 'PSScene4Band'	# PlanetScope Scenes
planet_itemtype_dict['psortho'] = 'PSOrthoTile'	# PlanetScope OrthoTiles
planet_itemtype_dict['reortho'] = 'REOrthoTile'	# RapidEye OrthoTiles
planet_itemtype_dict['reraw'] = 'REScene'	# RapidEye Scenes (unorthorectified strips)
planet_itemtype_dict['skysat'] = 'SkySatScene'	# SkySat Scenes
planet_itemtype_dict['L8_MS'] = 'Landsat8L1G'	# Landsat8 Scenes
planet_itemtype_dict['sentinel2_MS'] = 'Sentinel2L1C'	# Copernicus Sentinel-2 Scenes

boulder_geometry = {"type": "Polygon",
        "coordinates": [
          [
            [
              -105.4848861694336,
              39.98645910523746
            ],
            [
              -105.42755126953125,
              39.98645910523746
            ],
            [
              -105.42755126953125,
              39.99947896972128
            ],
            [
              -105.4848861694336,
              39.99947896972128
            ],
            [
              -105.4848861694336,
              39.98645910523746
            ]
          ]
        ]
      }
      

# filter for items the overlap with our chosen geometry
geometry_filter = {
  "type": "GeometryFilter",
  "field_name": "geometry",
  "config": boulder_geometry
}

# filter images acquired in a certain date range
start_date = "2017-08-01T00:00:00.000Z"
end_date = "2017-10-31T00:00:00.000Z"
date_range_filter = {
  "type": "DateRangeFilter",
  "field_name": "acquired",
  "config": {
    "gte": start_date,
    "lte": end_date
  }
}

# filter any images which are more than 50% clouds
cloud_c = 0.05
cloud_cover_filter = {
  "type": "RangeFilter",
  "field_name": "cloud_cover",
  "config": {
    "lte": cloud_c
  }
}

# create a filter that combines our geo and date filters
# could also use an "OrFilter"
boulder_ridge_road = {
  "type": "AndFilter",
  "config": [geometry_filter, date_range_filter, cloud_cover_filter]
}




# Stats API request object
stats_endpoint_request = {
  "interval": "day",
  "item_types": [planet_itemtype_dict['ps4'], planet_itemtype_dict['sentinel2_MS']],
  "filter": boulder_ridge_road
}

# fire off the POST request
result = \
  requests.post(
    'https://api.planet.com/data/v1/stats',
    auth=HTTPBasicAuth(os.environ['PL_API_KEY'], ''),
    json=stats_endpoint_request)

#print result.text


# Search API request object
search_endpoint_request = {
  #"item_types": [planet_itemtype_dict['ps4'], planet_itemtype_dict['sentinel2_MS']],
  "item_types": [planet_itemtype_dict['ps4']],
  "filter": boulder_ridge_road
}

result = \
  requests.post(
    'https://api.planet.com/data/v1/quick-search',
    auth=HTTPBasicAuth(os.environ['PL_API_KEY'], ''),
    json=search_endpoint_request)
    
print(result.text)
res_j = json.loads(result.text)




# ## optionally activate and download
# import os
# import requests

# #item_id = "20160707_195147_1057916_RapidEye-1"
# item_id = str(res_j['features'][0]['id'])
# item_type = planet_itemtype_dict['sentinel2_MS']
# asset_type = "analytic"

# # setup auth
# session = requests.Session()
# session.auth = (os.environ['PL_API_KEY'], '')

# # request an item
# item = \
  # session.get(
    # ("https://api.planet.com/data/v1/item-types/" +
    # "{}/items/{}/assets/").format(item_type, item_id))

# # extract the activation url from the item for the desired asset
# # if getting Sentinel 2 data, need asset_type + 'b<n>', n 1-10
# item_activation_url = item.json()[asset_type]["_links"]["activate"]

# # request activation
# response = session.post(item_activation_url)

# print response.status_code