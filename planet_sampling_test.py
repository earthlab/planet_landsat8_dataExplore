import os, sys
import rasterio
import numpy as np
from matplotlib import pyplot as plt
import glob
from xml.dom import minidom


def getCorrCoefs_planet(xml_file)

    xmldoc = minidom.parse(xml_file)
    nodes = xmldoc.getElementsByTagName("ps:bandSpecificMetadata")

    # XML parser refers to bands by numbers 1-4
    coeffs = {}
    for node in nodes:
        bn = node.getElementsByTagName("ps:bandNumber")[0].firstChild.data
        if bn in ['1', '2', '3', '4']:
            i = int(bn)
            value = node.getElementsByTagName("ps:reflectanceCoefficient")[0].firstChild.data
            coeffs[i] = float(value)
            


home_dir = r"C:\projects\rd\planet\planet_sample"

all_ims = glob.glob(home_dir, '*.tif')
all_xml = glob.glob(home_dir, '*.xml')

