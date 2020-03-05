import numpy  as np
import h5py
import pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
from matplotlib import ticker, cm
import pandas as pd
import os

def get_xco2_from_emission(INemission,INmonth,INboundarylayerheight):

    monthday = [31,28,31,30,31,30,\
                31,31,30,31,30,31]

    Av = 6.023e23 # Avagrado's number molecules/mole

    # number density of air
    rho_air = 2.69e25                       # molec/m^3 NOTE FIXED NEEDS ADJUSTING RE SURFACE PRESSURE
    rho_air = np.multiply(rho_air,1e9)      # molec/km^3    
    
    # NOTE: ODIAC described as tonnes C/1km box/month
    
    INemission = np.multiply(INemission,44/12.)              # tonnes CO2/1 km box/month

    dayhourconversion = monthday[INmonth-1]*24               # days/month * hours/day = hours/month

    INemission = np.divide(INemission,dayhourconversion)     # tonnes CO2/1 km box/hour
    
    INemission = np.multiply(INemission,1e6)                 # grams CO2/1 km box/hour

    INemission = np.divide(INemission,44)                    # mole CO2/1 km box/hour

    INemission = np.multiply(INemission,Av)                  # molec CO2/1 km box/hour

    INemission = np.divide(INemission,INboundarylayerheight) # molec CO2/km^3

    return np.divide(INemission,rho_air)
    


def readodiac(extent,monthchar):
    
    filename = 'odiac2018/1km_tif/odiac2018_1km_excl_intl_17'+monthchar+'.tif'

    import tifffile as tf

    image_stack = tf.imread(filename)
    #print(image_stack.shape)
    #print(image_stack.dtype)

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


def readecostress(GEOfilename,WUEfilename):

    data_geo = h5py.File(GEOfilename,'r')
    data_wue = h5py.File(WUEfilename,'r')

    lat   = np.array(data_geo['Geolocation/latitude'])
    lon   = np.array(data_geo['Geolocation/longitude'])

    wue   = np.array(data_wue['Water Use Efficiency/WUEavg'])
    
    return lat, lon, wue




def plotdata(ax,fig,lon,lat,extent,invar,colors,v_min,v_max,usecontour):

    ax.set_extent(extent)

    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    #ax.add_feature(land_50m, edgecolor='gray')


    if usecontour:

        im = ax.contourf(lon, lat, invar, levels=colors, vmin=v_min,vmax=v_max,\
                          cmap=plt.cm.get_cmap('rainbow'))

        out=fig.colorbar(im,ax=ax, extend='both',orientation='vertical',\
                         label='XXX')
        
    else:
    
        im = ax.scatter(lon, lat, c=colors, s=5, alpha=0.8, \
                        transform=ccrs.PlateCarree(),\
                        cmap=plt.cm.get_cmap('rainbow'),vmin=v_min,vmax=v_max)

        out=fig.colorbar(im,ax=ax, extend='both',orientation='vertical',\
                         label='XCO$_2$ (ppm)')
    
    #ax.add_feature(cfeature.BORDERS, edgecolor='gray')
    #ax.add_feature(cfeature.OCEAN)
    #ax.add_feature(cfeature.LAND, edgecolor='black')
    #ax.add_feature(cfeature.LAKES, edgecolor='black')
    #ax.add_feature(cfeature.RIVERS)
    #ax.add_feature(cfeature.COASTLINE)

    ax.coastlines(resolution='50m', color='black', linewidth=1)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False





#-------------------------
# Read ECOSTRESS data
#-------------------------

GEOfilename = '../ECOSTRESS/data/ECOSTRESS_L1B_GEO_06109_001_20190803T220442_0600_01.h5'
WUEfilename = '../ECOSTRESS/data/ECOSTRESS_L4_WUE_06109_001_20190803T220442_0600_01.h5'

es_lat, es_lon, es_wue = readecostress(GEOfilename,WUEfilename)

#-------------------------
# Read OCO-3 data
#-------------------------

samfile = 'SAM_TG.csv'
datadir = 'SAMS-250919/'
cvsdata = pd.read_csv(samfile,header=0)
cvsdata = np.array(cvsdata)
nsams = len(cvsdata[:,0])

for ii in np.arange(nsams):
    tmp           = cvsdata[ii,3].split('/')
    saminfilename = datadir+tmp[-1]
    sam_or_tg     = cvsdata[ii,7]
    orbitnumber   = cvsdata[ii,8]
    targetname    = cvsdata[ii,0]+'_'+orbitnumber

    print('**READ IN FILE: '+saminfilename+' Output name: '+targetname)

    # h5dump -n [filename]
    data_l2 = h5py.File(saminfilename,'r')
 
    lat   = np.array(data_l2['RetrievalGeometry/retrieval_latitude'])
    lon   = np.array(data_l2['RetrievalGeometry/retrieval_longitude'])
    flag  = np.array(data_l2['RetrievalResults/outcome_flag'])
    xco2  = np.array(data_l2['/RetrievalResults/xco2']) * 1.e6
    smodeB= np.array(data_l2['/RetrievalHeader/sounding_operation_mode'])
    pma   = np.array(data_l2['/RetrievalGeometry/retrieval_pma_motion_flag'])
    sid   = np.array(data_l2['RetrievalHeader/sounding_id'])



    #-------------------------
    # Identify and plot SAM in OCO-3 data
    #-------------------------

    # AM = area mode; TG = target

    ind = -999
    
    if sam_or_tg == 'TG': ind = np.where((pma == 0) & (smodeB ==  b'TG') & (flag <= 2) & (sid != 0))
    if sam_or_tg == 'AM': ind = np.where((pma == 0) & (smodeB ==  b'AM') & (flag <= 2) & (sid != 0))


    #
    # Error checking - ignore files that do not contain TG or AM flags
    #
    if (np.array(ind).size) != 0:

        central_lat = 0
        central_lon = 0
        extent      = [np.min(lon[ind]), np.max(lon[ind]), \
                       np.min(lat[ind]), np.max(lat[ind])]

        proj = ccrs.PlateCarree(central_longitude=central_lon)


        #
        # Remove bottom 10th percentile to identify city signal
        #
        bttenpc  = np.percentile(xco2[ind],10)
        xco2  = np.subtract(xco2,bttenpc)


        plt.figure(20,figsize=[10,10])
        
        plt.subplot(211)

        minval = np.min(xco2[ind])
        maxval = np.max(xco2[ind])
        nbins  = 100
        dbin   = (maxval-minval)/np.float(nbins)
        bins = np.arange(nbins)*dbin + minval
    
        plt.hist(xco2[ind],bins)

        

        from mpl_toolkits.mplot3d import Axes3D 
    
        plt.subplot(212,projection='3d')

        p = plt.scatter(lon[ind],lat[ind],xco2[ind],\
                        vmin = -10, vmax = 10,\
                        c=xco2[ind],marker='o',alpha=0.5,cmap=cm.coolwarm)
        plt.colorbar(p)
        plt.xlabel('Longitude [deg]')
        plt.ylabel('Latitude [deg]')
        plt.title(targetname)


        print('******************')
        print(saminfilename)
        print(extent)

        plt.savefig('pngfigures/hist.'+targetname+'.png')

        plt.close()

        fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, sharex=True, \
                                         sharey=True, figsize=(15,10),\
                                          subplot_kw={'projection': proj})

        v_min = np.median([xco2[ind]]) - np.std([xco2[ind]])
        v_max = np.median([xco2[ind]]) + np.std([xco2[ind]])

        colors = xco2[ind]
    
        plotdata(ax1,fig1,lon[ind],lat[ind],extent,xco2[ind],colors,v_min,v_max,usecontour=0)

        tokyo    = [35.6804,139.7690]
        tsukuba  = [36.0835,140.0764]
        kawasaki = [35.5308,139.7029]

        ax1.plot(tokyo[1],tokyo[0],'ko')
        ax1.text(tokyo[1],tokyo[0],'Tokyo')
        
        ax1.plot(tsukuba[1],tsukuba[0],'ko')
        ax1.text(tsukuba[1],tsukuba[0],'Tsukuba')
    
        ax1.plot(kawasaki[1],kawasaki[0],'ko')
        ax1.text(kawasaki[1],kawasaki[0],'Kawasaki')
    
        ax1.set_title(saminfilename)






        #-------------------------
        # Read ODIAC data
        #-------------------------

        co2flux, odiaclons, odiaclats = readodiac(extent,'01')

        #-------------------------
        # Plot ODIAC data
        #-------------------------

        X, Y = np.meshgrid(odiaclons,odiaclats)

        odiacdatapoints = np.count_nonzero(~np.isnan(co2flux))
        if odiacdatapoints > 0:

            surf = ax3.contourf(X,Y,co2flux,locator=ticker.LogLocator(),\
                                transform=ccrs.PlateCarree(),cmap=cm.coolwarm)

            gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                               linewidth=2, color='gray', alpha=0.5, linestyle='--')

            cb = fig1.colorbar(surf,ax=ax3,extend='both',orientation = 'vertical',\
                               label = 'ODIAC CO$_2$ emissions [Tonne Carbon/cell/month]')

            gl.xlabels_top = False
            gl.ylabels_right = False

            ax3.coastlines(resolution='50m', color='black', linewidth=1)

            tokyo    = [35.6804,139.7690]
            tsukuba  = [36.0835,140.0764]
            kawasaki = [35.5308,139.7029]
        
            ax3.plot(tokyo[1],tokyo[0],'ko')
            ax3.text(tokyo[1],tokyo[0],'Tokyo')

            ax3.plot(tsukuba[1],tsukuba[0],'ko')
            ax3.text(tsukuba[1],tsukuba[0],'Tsukuba')

            ax3.plot(kawasaki[1],kawasaki[0],'ko')
            ax3.text(kawasaki[1],kawasaki[0],'Kawasaki')

            #
            # Plot ODIAC XCO2 after making some assumptions
            #

            usemonth               = 8 # August
            useboundarylayerheight = 2 # km

            print(np.nanmin(co2flux),np.nanmax(co2flux),np.nanmedian(co2flux))
            
            odiac_xco2 = get_xco2_from_emission(co2flux,usemonth,useboundarylayerheight)

            odiac_xco2 = np.divide(odiac_xco2,1e-6)

            v_min = 0
            v_max = 5
            npts  = 100
            uselevels = np.arange(npts)*(v_max-v_min)/npts
            
            surfxco2 = ax4.contourf(X,Y,odiac_xco2,levels=uselevels,\
                                    transform=ccrs.PlateCarree(),\
                                    vmin=v_min,vmax=v_max,\
                                    cmap=cm.coolwarm)

            gl2 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                               linewidth=2, color='gray', alpha=0.5, linestyle='--')

            cb2 = fig1.colorbar(surfxco2,ax=ax4,extend='both',orientation = 'vertical',\
                                label = 'ODIAC boundary layer XCO2 [ppm]')

            gl2.xlabels_top = False
            gl2.ylabels_right = False

            ax4.coastlines(resolution='50m', color='black', linewidth=1)


            plt.savefig('pngfigures/maps.'+targetname+'.png')

            #plt.show(); sys.exit()

        plt.close()
    
    
    
plt.show()
sys.exit()





#
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#
#X, Y = np.meshgrid(lon[ind],lat[ind])
#
#ax2 = fig1.gca(projection='3d')
#
#Z = np.meshgrid(xco2[ind],xco2[ind])
#
#surf = ax2.plot_surface(X,Y,Z,cmap=cm.coolwarm)
#
#




##-------------------------
## Crudely plot ECOSTRESS data
##-------------------------
#
#v_min = np.nanmedian(es_wue) - 2*np.nanstd(es_wue)
#v_max = np.nanmedian(es_wue) + 2*np.nanstd(es_wue)
#
#npts = 10
#colors = np.arange(npts)*(v_max-v_min)/npts + v_min
#
#
#plotdata(ax2,fig1,es_lon,es_lat,\
#         extent,es_wue,colors,v_min,v_max,usecontour=1)
#
#ax2.set_title(GEOfilename)
#    
plt.show()






