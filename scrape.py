from __future__ import print_function
import numpy
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import pylab
import astropy as astropy
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import os.path


#######################################################
#                                                     #
#  SOURCE: https://ps1images.stsci.edu/ps1image.html  #
#                                                     #
#######################################################


def parsebcgs(file):
    raw_data = open(file)
    raw_data = raw_data.read()
    lines = raw_data.split("\n")

    data = []
    for i in range(len(lines)):
        elem = lines[i].split()
        data.append(elem)

    bcgcoords = []
    for i in range(len(data)):
        bcgcoords.append([])
        bcgcoords[i].append(data[i][0].lower())
        if bcgcoords[i][0][-1] == '.':
            bcgcoords[i][0] = bcgcoords[i][0][:-1]

        if data[i][3][-1] == '.':
            data[i][3] = data[i][3][:-1]

        if data[i][6][-1] == '.':
            data[i][6] = data[i][6][:-1]


        bcgcoords[i].append(str(data[i][1])+ " " +str(data[i][2]) + " " +str(data[i][3])+ " " + str(data[i][4]) + " " + str(data[i][5]) + " " + str(data[i][6]))

    for i in range(len(bcgcoords)):
        c = SkyCoord(bcgcoords[i][1], unit=(u.hourangle, u.deg))
        bcgcoords[i][1] = c.ra.degree
        bcgcoords[i].append(c.dec.degree)

    return bcgcoords

def getimages(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


###########################################################################################################################


bcgcoords = parsebcgs("/home/elenarom/galaxy_alignment/bcgs.txt")[159:] ### Change this number to change which are included

### For future use: set bcgcoords to a list of [name, ra, dec], where ra and dec are in decimal degrees


for i in range(len(bcgcoords)):
    ra = bcgcoords[i][1]
    dec = bcgcoords[i][2]
    fitsurl = geturl(ra, dec, size=1000, filters="r", format="fits") ### Edit this to change image parameters
    
    if len(fitsurl) == 1:
        hdul = fits.open(fitsurl[0])
        name = "/home/elenarom/galaxy_alignment/fitsoutput/" + bcgcoords[i][0] + ".fits" ### Change path to wherever you want the fits files to be saved
        if not os.path.exists(name):
            fits.writeto(name, hdul[0].data, header=hdul[0].header)
        else:
            print("already have "+ bcgcoords[i][0]) ### if file with this name already exists
    else:
        print("missed " + str(bcgcoords[i][0])) ### if no image found in the PanSTARRS catalog



