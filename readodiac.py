import numpy as np
import pylab as plt
import sys
from matplotlib import ticker, cm


#2) Coordinate definition
#
#The emission fields are defined on a geographical (longitude and latitude) coordinate system.  Longitude and latitude can be calculated as follows: 
#
#Lon (i) = -180 + dx/2 + dx x i (i=0,43199)
#Lat (j) =  90 + dy/2 + dy x j (j=0,21599)
#
#where dx = dy = 30 arc-seconds (0.008333333333333 degrees or approximately 1km). 




def readodiac(extent,monthchar):
    
    filename = 'odiac2018/1km_tif/odiac2018_1km_excl_intl_17'+monthchar+'.tif'

    import tifffile as tf

    image_stack = tf.imread(filename)

    #print(image_stack.shape)
    #print(image_stack.dtype)
    
    #image_stack = np.flip(image_stack,axis=1)
    
    print(image_stack.shape)
    print(image_stack.dtype)

    nrows = len(image_stack)
    ncols = len(image_stack[1])

    Unit  = 'Tonne Carbon/cell/month'

    dx = dy = 0.008333333333333

    ODIAClons = np.arange(ncols)*dx + dx/2. - 180.
    ODIAClats = 90 - np.arange(nrows)*dy + dy/2. 

    ind = np.where(image_stack == 0.)
    image_stack[ind] = np.nan

    #Focus on target region

    indlon = np.where((ODIAClons >= extent[0]) & (ODIAClons <= extent[1]))
    
    indlat = np.where((ODIAClats >= extent[2]) & (ODIAClats <= extent[3]))

    image_stack = np.squeeze(image_stack[:,indlon])
    image_stack = np.squeeze(image_stack[indlat,:])

    return image_stack, ODIAClons[indlon], ODIAClats[indlat]





extent = [0,50,30,60]

co2flux, odiaclons, odiaclats = readodiac(extent,'01')

X, Y = np.meshgrid(odiaclons,odiaclats)

print(np.max(co2flux),np.min(co2flux))

cflux = plt.contourf(X,Y,co2flux,locator=ticker.LogLocator(),cmap=cm.coolwarm)

plt.show()
sys.exit()
