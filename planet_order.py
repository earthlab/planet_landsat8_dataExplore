import planet
import os, sys
import requests

ids = [ 'PSScene3Band:20170813_170429_100e', 
        'PSScene4Band:20170810_170719_101a', 
        'PSScene4Band:20170810_170718_101a', 
        'PSScene3Band:20170810_170719_101a', 
        'PSScene3Band:20170810_170718_101a', 
        'PSScene4Band:20170812_170427_1032', 
        'PSScene4Band:20170812_170426_1032', 
        'PSScene3Band:20170812_170427_1032', 
        'PSScene3Band:20170812_170426_1032', 
        'PSScene4Band:20170812_170517_0f3f', 
        'PSScene3Band:20170812_170517_0f3f', 
        'PSScene4Band:20170813_170429_100e', 
        'PSScene4Band:20170827_175918_0f24', 
        'PSScene3Band:20170827_175918_0f24', 
        'PSScene4Band:20170829_170647_0f34', 
        'PSScene4Band:20170829_170646_0f34', 
        'PSScene3Band:20170829_170647_0f34', 
        'PSScene3Band:20170829_170646_0f34', 
        'PSScene4Band:20170903_170641_0f17', 
        'PSScene4Band:20170903_170640_0f17', 
        'PSScene3Band:20170903_170641_0f17', 
        'PSScene3Band:20170903_170640_0f17', 
        'PSScene4Band:20170906_170654_1012', 
        'PSScene4Band:20170906_170653_1012', 
        'PSScene3Band:20170906_170654_1012', 
        'PSScene3Band:20170906_170653_1012', 
        'PSScene4Band:20170915_170529_100a', 
        'PSScene3Band:20170915_170529_100a', 
        'PSScene4Band:20170918_175620_0f21', 
        'PSScene3Band:20170918_175620_0f21', 
        'PSScene4Band:20170920_175654_1054', 
        'PSScene4Band:20170920_175653_1054', 
        'PSScene3Band:20170920_175654_1054', 
        'PSScene3Band:20170920_175653_1054', 
        'PSScene4Band:20170921_175651_101c', 
        'PSScene3Band:20170921_175651_101c', 
        'PSScene4Band:20171005_175435_0f4d', 
        'PSScene3Band:20171005_175435_0f4d', 
        'PSScene3Band:20171008_175414_104c', 
        'PSScene3Band:20171008_175413_104c', 
        'PSScene4Band:20171008_175414_104c', 
        'PSScene4Band:20171008_175413_104c', 
        'PSScene3Band:20171011_191452_0c44', 
        'PSScene3Band:20171012_170559_100a', 
        'PSScene4Band:20171012_170559_100a', 
        'PSScene4Band:20171012_175409_1053', 
        'PSScene4Band:20171012_175408_1053', 
        'PSScene3Band:20171012_175409_1053', 
        'PSScene3Band:20171012_175408_1053', 
        'PSScene3Band:20171013_164837_1_0c65', 
        'PSScene3Band:20171013_164837_0c65', 
        'PSScene4Band:20171013_164837_1_0c65', 
        'PSScene4Band:20171013_164837_0c65', 
        'PSScene3Band:20171017_170859_1034', 
        'PSScene3Band:20171017_170858_1034', 
        'PSScene4Band:20171017_170859_1034', 
        'PSScene4Band:20171017_170858_1034', 
        'PSScene3Band:20171019_170725_1032', 
        'PSScene4Band:20171019_170725_1032', 
        'PSScene4Band:20171024_170922_1003', 
        'PSScene3Band:20171024_170922_1003', 
        'PSScene4Band:20171024_175217_1050', 
        'PSScene4Band:20171024_175216_1050', 
        'PSScene3Band:20171024_175217_1050', 
        'PSScene3Band:20171024_175216_1050', 
        'PSScene3Band:20171029_170746_1025', 
        'PSScene4Band:20171029_170746_1025', 
        'PSScene4Band:20171029_170956_1005', 
        'PSScene3Band:20171029_170956_1005', 
        'PSScene4Band:20171029_175157_0f46', 
        'PSScene3Band:20171029_175157_0f46', 
        'Sentinel2L1C:S2A_MSIL1C_20171029T180251_N0206_R141_T13TDE_20171029T193005', 
        'PSScene3Band:20171029_183018_0c44']

ids_3band = [i.split(":")[1] for i in ids if '3Band' in i]
ids_4band = [i.split(":")[1] for i in ids if '4Band' in i]
ids_sentinel = [i.split(":")[1] for i in ids if 'Sentinel' in i]

item_id = ids_4band[0]
# item_id =  "20160707_195147_1057916_RapidEye-1"
item_type = "PSScene4Band"
asset_type = "analytic"

# setup auth
session = requests.Session()
#session.auth = (os.environ['PL_API_KEY'], '')
session.auth = ('9990896fbd9e46dfac87ad0119a3cd33', '')

# request an item
item = \
  session.get(
    ("https://api.planet.com/data/v1/item-types/" +
    "{}/items/{}/assets/").format(item_type, item_id))

# extract the activation url from the item for the desired asset
item_activation_url = item.json()[asset_type]["_links"]["activate"]

# request activation
response = session.post(item_activation_url)

print response.status_code

curl_statement = '''curl -L -H "Authorization: api-key {}" https://api.planet.com/data/v1/item-types/{}/items/{}/assets/ | jq .analytic.status'''.format('9990896fbd9e46dfac87ad0119a3cd33', item_type, item_id)