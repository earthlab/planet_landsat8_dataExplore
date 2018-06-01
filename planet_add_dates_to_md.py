import arcpy
import os
from datetime import datetime
import glob
from xml.dom import minidom

# function to return collection dates from Planet MultiSpectral image metadata
def getDates_planet(xml_file):

    xmldoc = minidom.parse(xml_file)
    date_node = xmldoc.getElementsByTagName("eop:acquisitionDate")
    date_str = str(date_node[0].firstChild.nodeValue)
    #acq_date_yyyymmdd = date_str.split('T')[0]        
    test = '2017-10-12T17:54:08+00:00'
    part1 = date_str.split('T')[0]
    part2 = date_str.split('T')[1].split('+')[0]
    rfm = ' '.join([i for i in [part1,part2]])
    
    return rfm #acq_date_yyyymmdd
    
# extract the data. it will be a bunch of lists
home_dir = r"C:\Projects\RD\planet\sample_order"

all_ims = glob.glob('{}/*/*.tif'.format(home_dir))
all_sr_ims = [im for im in all_ims if "SR" in im]
all_xml = [s.replace('MS_SR.tif', 'MS_metadata.xml') for s in all_sr_ims]

dates = [getDates_planet(xml) for xml in all_xml]

# apply the dates as a field to the mosaic dataset
filegdb = r"C:\Projects\RD\planet\planet_sample.gdb"
md_name = 'planet_sample_SR'
md = os.path.join(filegdb, md_name)

# arcpy.AddField_management(md, 'AcqDate', 'DATE')

# update the field with the dates array
with arcpy.da.UpdateCursor(fp, ['AcqDate']) as uc:
    for i,row in enumerate(uc):
        imtime = dates[i]
        # dt = datetime.strptime(imtime, '%Y-%m-%d')
        dt = datetime.strptime(imtime, '%Y-%m-%d %H:%M:%S')
        row[0] = dt
        uc.updateRow(row)
        
