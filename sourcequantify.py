import numpy as np
import pylab as plt
import cartopy.crs as ccrs
import sys
from xco2funcs import *
from haversine import haversine

test_outfile = 'sourceestimationtestfile.npz'


#np.savez(test_outfile,oco3co2=newgrid,spressoco3=spress_oco3_grid,\
#         odiacco2emission=odiacco2flux,odiacxco2=odiac_xco2,\
#         uwind=UWIND,vwind=VWIND,\
#         lons=ODIAClonsDOMAIN,lats=ODIAClatsDOMAIN)


npzfile       = np.load(test_outfile)
oco3xco2      = npzfile['oco3co2']
spress        = npzfile['spressoco3']
odiacco2flux  = npzfile['odiacco2emission']
odiacxco2     = npzfile['odiacxco2']
lons          = npzfile['lons']
lats          = npzfile['lats']
uwind         = npzfile['uwind']
vwind         = npzfile['vwind']



targets = [[51.6660,32.6539],[49.7013,34.0954]]

X, Y = np.meshgrid(lons, lats)




def plotdata(invar,lat,lon,label,vmin,vmax):


    surfvar = plt.pcolormesh(lat, lon,invar,\
                             vmin=v_min,vmax=v_max,\
                             cmap=plt.cm.get_cmap('rainbow'))

    cb3 = plt.colorbar(surfvar,extend='both',orientation = 'horizontal',\
                       label = label)
    cb3.set_clim(vmin,vmax)


print(np.nanmedian(odiacco2flux))
print(np.nanpercentile(odiacco2flux,10))

#
# Further define the city signal
#
ind = np.where(odiacco2flux <= np.nanpercentile(odiacco2flux,75))
oco3xco2[ind]     = np.nan
odiacco2flux[ind] = np.nan
odiacxco2[ind]    = np.nan
uwind[ind]        = np.nan
vwind[ind]        = np.nan


useextent = [50,52,35,36.5]





plt.subplot(221)
v_min = np.nanmedian([oco3xco2]) - np.nanstd([oco3xco2])
v_max = np.nanmedian([oco3xco2]) + np.nanstd([oco3xco2])
plotdata(oco3xco2,X,Y,'OCO-3 $\Delta$XCO$_2$ (ppm)', v_min,v_max)
plt.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')

plt.subplot(222)
v_min = 10; v_max = 1e3
plotdata(odiacco2flux,X,Y,'ODIAC CO2 emissions [gC/km2/s]', v_min,v_max)
plt.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')

plt.subplot(223)
v_min = 0; v_max = 3
plotdata(odiacxco2,X,Y,'ODIAC XCO2 [ppm]', v_min,v_max)
plt.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')

plt.subplot(224)
Q = plt.quiver(X,Y,uwind, vwind, \
           scale=5,scale_units='inches',color='black',alpha=0.5)
qk = plt.quiverkey(Q, 51.75, 32.1, 1, r'$1 \frac{m}{s}$', labelpos='E',
                   coordinates='data')
plt.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')


#
# Source pixel method
#

Qsum = 0.
for ii in np.arange(len(uwind)):
    for jj in np.arange(len(vwind)):

        uwinduse      = uwind[ii,jj]
        vwinduse      = vwind[ii,jj]
        distance      = 1000 # m
        odiacxco2use  = odiacxco2[ii,jj]
        spressuse     = spress[ii,jj]

        print(uwinduse,vwinduse,distance,odiacxco2use,spressuse)
        
        if ~np.isnan(odiacxco2use):
            Q = sourcepixelmethod(uwinduse,vwinduse,distance,odiacxco2use,spressuse)
            #print(Q)
            Qsum += Q

Qsum = Qsum * 3600 * 24 * 365 / 1e12
print(Qsum)

            
#
# Preliminary IME
# An emission estimate for Isfahan is 19.8 +/- 6.9 Mt CO2
# (for the year 2013, http://citycarbonfootprints.info/).

meanwindspeed = 5    # m/s
Ldimension    = 10 # km

sdist = np.zeros([len(lons),len(lats)]) # sample dist
dist = np.zeros([len(lons),len(lats)])
dist[:,:] = np.nan
sdist[:,:] = np.nan
targetcentre=[32.6539,51.6660]
for ii in np.arange(len(lats)):
    for jj in np.arange(len(lons)):
        dist[ii,jj] = haversine(targetcentre,[lats[ii],lons[jj]])

ind = np.where(dist <= 20)
print(np.shape(ind))
print(ind)
sdist[ind] = 1


wind_mask = np.zeros([len(lons),len(lats)])
wind_mask[:,:] = np.nan
wind_mask[ind] = np.sqrt(uwind[ind]**2 + vwind[ind]**2)
print(np.shape(wind_mask))

print(np.nanmedian(wind_mask))

plt.figure(20)
img = plt.pcolor(X,Y,wind_mask)
plt.colorbar(img)
plt.plot(targetcentre[1],targetcentre[0],'ko',markersize=80,markerfacecolor='none')    
plt.show()
sys.exit()

difflon = np.subtract(lons,targetcentre[0])
difflat = np.subtract(lats,targetcentre[1])



ODIACEMISSION, OCO3_IME = calc_ime(odiacco2flux,oco3xco2,meanwindspeed,Ldimension)

imestring   = writescreenstring(OCO3_IME,0)
odiacstring = writescreenstring(ODIACEMISSION,0)
print('IME/OCO-3: '+imestring)
print('ODIAC: '+odiacstring)

plt.show()
