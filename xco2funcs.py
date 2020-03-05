import numpy  as np
import h5py
import pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
from matplotlib import ticker, cm
import pandas as pd
import os


def readxco2(InFileName):
    
    datestring    = InFileName[46:52]
    obsyear       = datestring[0:2]
    obsmonth      = datestring[2:4]
    obsdom        = datestring[4:6]

    # h5dump -n [filename]
    data_l2 = h5py.File(InFileName,'r')

    lat   = np.array(data_l2['RetrievalGeometry/retrieval_latitude'])
    lon   = np.array(data_l2['RetrievalGeometry/retrieval_longitude'])
    flag  = np.array(data_l2['RetrievalResults/outcome_flag'])
    xco2  = np.array(data_l2['/RetrievalResults/xco2']) * 1.e6
    smodeB= np.array(data_l2['/RetrievalHeader/sounding_operation_mode'])
    pma   = np.array(data_l2['/RetrievalGeometry/retrieval_pma_motion_flag'])
    sid   = np.array(data_l2['RetrievalHeader/sounding_id'])
    co2_ratio_idp                = np.array(data_l2['/PreprocessingResults/co2_ratio_idp'])
    h2o_ratio_idp                = np.array(data_l2['/PreprocessingResults/h2o_ratio_idp'])
    brdf_reflectance_o2          = np.array(data_l2['/BRDFResults/brdf_reflectance_o2'])
    brdf_reflectance_slope_o2    = np.array(data_l2['/BRDFResults/brdf_reflectance_slope_o2'])
    brdf_reflectance_weak_co2    = np.array(data_l2['/BRDFResults/brdf_reflectance_weak_co2'])
    brdf_reflectance_strong_co2  = np.array(data_l2['/BRDFResults/brdf_reflectance_strong_co2'])    
    aerosol_total_aod            = np.array(data_l2['/AerosolResults/aerosol_total_aod'])
    surface_pressure_fph         = np.array(data_l2['/RetrievalResults/surface_pressure_fph'])
    surface_pressure_apriori_fph = np.array(data_l2['/RetrievalResults/surface_pressure_apriori_fph'])
    co2_vertical_gradient_delta  = np.array(data_l2['/RetrievalResults/co2_vertical_gradient_delta'])
    aerosol_aod                  = np.array(data_l2['/AerosolResults/aerosol_aod'])
    post_o2_column               = np.array(data_l2['/RetrievalResults/retrieved_o2_column'])

    air_column = np.zeros(len(post_o2_column))
    air_column = np.divide(post_o2_column,0.20947)

    dP_fph = np.subtract(surface_pressure_fph, surface_pressure_apriori_fph)

    # NOTES from Rob:
    # 
    # The AODs are a pain to read. You’ll notice how
    # /AerosolResults/aerosol_aod is n_soundings x 8 x 4. The 8 are the 8
    # aerosol types and the 4 are 4 parameters describing the Gaussian
    # aerosol layer. The first of those 4 is the AOD for each type. The 8
    # possible types are given by /Metadata/AllAerosolTypes (in the same
    # order), so you’ll want the 7th (“ST” = “strat”) and 8th (“Water” =
    # “water”) values of /AerosolResults/aerosol_aod's second dimension, and
    # the 1st value of the third dimension. 

    total_aod_strat = np.array(np.squeeze(aerosol_aod[:,6,0])) # strat
    total_aod_water = np.array(np.squeeze(aerosol_aod[:,7,0])) # water

    # NOTE
    # Alternative to SAM is target
    # if sam_or_tg == 'TG': ind = np.where((pma == 0) & (smodeB ==  b'TG') & (flag <= 2) & (sid != 0))

    #
    # Use filters from Mattheus: r_01 and r_02
    #
    #co2_ratio_idp.between(1.005,1.058)
    #h2o_ratio_idp.between(0.8,1.025)
    #brdf_reflectance_o2.between(0,0.5)
    #brdf_reflectance_slope_o2.between(-0.00006,0)
    #brdf_reflectance_slope_strong_co2.between(-0.00014,0.00035)
    #aerosol_total_aod.between(0.0,0.15)
    #dP_fph.between(-12,12)
    #total_aod_strat.between(0,0.025)
    #total_aod_water.between(0,0.035)
    #
    #co2_ratio_idp.between(1.0075,1.04)
    #h2o_ratio_idp.between(0.87,1.04)
    #brdf_reflectance_o2.between(0.13,0.375)
    #co2_vertical_gradient_delta.between(-75,60)
    #brdf_reflectance_weak_co2.between(0.16,0.60)
    #brdf_reflectance_strong_co2.between(0.16,0.45)
    #dP_fph.between(-15,3)
    
    ind = np.where((pma == 0) & (smodeB ==  b'AM') & (flag <= 2) & (sid != 0) & \
                   ((co2_ratio_idp >= 1.005) & (co2_ratio_idp <= 1.058)) & \
                   ((h2o_ratio_idp >= 0.8) & (h2o_ratio_idp <=1.025))    & \
                   ((brdf_reflectance_o2 >= 0) & (brdf_reflectance_o2 <= 0.5)) & \
                   ((brdf_reflectance_slope_o2 >= -0.00006) & (brdf_reflectance_slope_o2 <= 0)) & \
                   #((brdf_reflectance_strong_co2 >= -0.00014) & (brdf_reflectance_strong_co2 <= 0.00035)) & \
                   #((aerosol_total_aod >= 0.0) & (aerosol_total_aod <= 0.15 )) & \
                   #((dP_fph >= -12) & (dP_fph <= 12)) & \
                   ((total_aod_strat >= 0) & (total_aod_strat <= 0.025)) & \
                   ((total_aod_water >= 0) & (total_aod_water <= 0.035)))


    #ind = np.where((pma == 0) & (smodeB ==  b'AM') & (flag <= 2) & (sid != 0) & \
    #               ((co2_ratio_idp >= 1.0075) & (co2_ratio_idp <= 1.04)) & \
    #               ((h2o_ratio_idp >= 0.87) & (h2o_ratio_idp <=1.04))    & \
    #               ((brdf_reflectance_o2 >= 0.13) & (brdf_reflectance_o2 <= 0.375)) & \
    #               ((co2_vertical_gradient_delta >= -75) & (co2_vertical_gradient_delta <= 0.60)) & \
    #               ((brdf_reflectance_weak_co2 >= 0.16) & (brdf_reflectance_weak_co2 <= 0.60)) & \
    #               ((brdf_reflectance_strong_co2 >= 0.16) & (brdf_reflectance_strong_co2 <= 0.45)) & \
    #               ((dP_fph >= -15) & (dP_fph <= 3)))
    


    xco2use = xco2[ind]
    latuse  = lat[ind]
    lonuse  = lon[ind]
    aircolumn  = air_column[ind]
    surface_pressure_fph = surface_pressure_fph[ind]

    return xco2use, latuse, lonuse, obsmonth, surface_pressure_fph, aircolumn
    

def convertodiacemissionunits(inQ,INmonth):
    #input units are from ODIAC - tonnes C/1km box/month
    #for the Gaussian plume model we want g/s/grid square

    # Need to account for area

    monthday = [31,28,31,30,31,30,\
                31,31,30,31,30,31]
     
    inQ = np.multiply(inQ,44/12.)                         # tonnes CO2/1 km box/month

    dayhourconversion = monthday[INmonth-1]*24*3600       # days/month * secs/day = sec/month

    inQ = np.divide(inQ,dayhourconversion)                # tonnes CO2/1 km box/sec
    
    inQ = np.multiply(inQ,1e6)                            # grams CO2/1 km box/sec

    return inQ
    

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
    
    #INemission = np.divide(INemission,INboundarylayerheight) # molec CO2/km^3

    return np.divide(INemission,rho_air)
    

def readodiac(monthchar):

    extent = [-180,180,-90,90]
    
    filename = 'odiac2018/odiac2018_1km_excl_intl_17'+monthchar+'.tif'

    import tifffile as tf

    image_stack = tf.imread(filename)
    #print(image_stack.shape)
    #print(image_stack.dtype)

    nrows = len(image_stack)
    ncols = len(image_stack[1])

    #print(nrows,ncols)
    #sys.exit()

    Unit  = 'Tonne Carbon/cell/month'

    dx = dy = 0.008333333333333

    ODIAClons = np.arange(ncols)*dx + dx/2. - 180.
    ODIAClats = 90 - np.arange(nrows)*dy + dy/2. 

    ind = np.where(image_stack == 0.)
    image_stack[ind] = np.nan

    #Focus on target region

    #indlon = np.where((ODIAClons >= extent[0]) & (ODIAClons <= extent[1]))
    #indlat = np.where((ODIAClats >= extent[2]) & (ODIAClats <= extent[3]))

    #image_stack = np.squeeze(image_stack[:,indlon])
    #image_stack = np.squeeze(image_stack[indlat,:])

    return image_stack, ODIAClons, ODIAClats


def plotdata(ax,fig,lon,lat,extent,invar,colors,v_min,v_max,uselabel,usecontour,no2flag):

    ax.set_extent(extent)

    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    
    if no2flag:

        im = ax.scatter(lon, lat, c=colors, s=10, alpha=0.8, \
                        transform=ccrs.PlateCarree(),\
                        cmap=plt.cm.get_cmap('rainbow'),vmin=v_min,vmax=v_max)

    else:

        X, Y = np.meshgrid(lon,lat)
        
        im = ax.pcolormesh(X,Y,invar,\
                       transform=ccrs.PlateCarree(),\
                       vmin=v_min,vmax=v_max,\
                       cmap=plt.cm.get_cmap('rainbow'))


    out=fig.colorbar(im,ax=ax, extend='both',orientation='horizontal',\
                     label=uselabel)        

    ax.coastlines(resolution='50m', color='black', linewidth=1)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False



def definecommongrid(extent):
    
    nrows = 21600; ncols = 43200
    dx = dy = 0.008333333333333

    ODIAClons = np.arange(ncols)*dx + dx/2. - 180.
    ODIAClats = 90 - np.arange(nrows)*dy + dy/2.
    
    indlon = np.where((ODIAClons >= extent[0]) & (ODIAClons <= extent[1]))
    indlat = np.where((ODIAClats >= extent[2]) & (ODIAClats <= extent[3]))
    ODIAClons = ODIAClons[indlon]; ODIAClats = ODIAClats[indlat]
    
    return ODIAClons, ODIAClats


def regridoco3xco2(ODIAClonsDOMAIN, ODIAClatsDOMAIN,xco2use,latxco2,lonxco2):

    newgrid = np.empty([len(ODIAClonsDOMAIN),len(ODIAClatsDOMAIN)]); newgrid[:,:] = np.nan
    nobs    = np.zeros([len(ODIAClonsDOMAIN),len(ODIAClatsDOMAIN)])

    for ii in np.arange(0,len(xco2use),1):
        nearlat = np.subtract(ODIAClatsDOMAIN,latxco2[ii])
        nearlon = np.subtract(ODIAClonsDOMAIN,lonxco2[ii])
        p = np.squeeze(np.min(abs(nearlat)))
        q = np.squeeze(np.min(abs(nearlon)))
        pind = np.squeeze(np.where(abs(nearlat) == p))
        qind = np.squeeze(np.where(abs(nearlon) == q))
        newgrid[qind-1,pind-1] = xco2use[ii]
        nobs[qind-1,pind-1] += 1

    return newgrid, nobs


def regridmet(ODIAClonlatDOMAIN, wind,latlon):

    newgrid = np.empty([len(ODIAClonlatDOMAIN)]); newgrid[:] = np.nan
    nobs    = np.zeros([len(ODIAClonlatDOMAIN)])

    for ii in np.arange(0,len(latlon),30):

        nearlatlon = np.subtract(ODIAClonlatDOMAIN,latlon[ii])
        p = np.squeeze(np.min(abs(nearlatlon)))
        pind = np.squeeze(np.where(abs(nearlatlon) == p))
            
        newgrid[pind-1] = wind[ii]
        nobs[pind-1] += 1
    
    return newgrid, nobs

def filterandconvertno2(no2,latno2,lonno2,qa):
    #
    # NO2
    #
    no2    = np.reshape(no2,-1)
    latno2 = np.reshape(latno2,-1)
    lonno2 = np.reshape(lonno2,-1)
    qa     = np.reshape(qa,-1)

    no2ind = np.where((qa >= 0.5) & (no2 < 1e35))

    # convert mol/m2 to molec/cm2
    #
    no2 = np.multiply(no2,6.023e23)
    no2 = np.divide(no2,1e4)
    no2 = np.divide(no2,1e16) # 10^16 molec/cm2

    return no2[no2ind], latno2[no2ind], lonno2[no2ind]

def writescreenstring(a, b): return '{:6.1f}'.format(a)+'$\pm$'+'{:6.1f}'.format(b)+' TgCO2/yr'    


def calc_ime(odiacco2flux,xco2grid,aircolumngrid,MEANWINDSPEED,Ldimension):

    # ODIAC CO2 fluxes have units of gCO2/km2/s
    # Return TgCO2/yr

    Mw_air = 28.9628
    Mw_co2 = 44.01
    Mw_ratio = Mw_co2/Mw_air
    Av = 6.023e23                           # Avagrado's number molecules/mole        
        
    xco2grid      = np.reshape(xco2grid,-1)
    aircolumngrid = np.reshape(aircolumngrid,-1)

    ind = np.where(~np.isnan(xco2grid) & (xco2grid>0)) # remove the second criterion when L2 filter is correctly applied
    xco2grid = xco2grid[ind]
    aircolumngrid = aircolumngrid[ind]

    
    # Convert DXCO2 (ppb) to mole fraction
    tmp = np.divide(xco2grid,1e6)

    # Convert molec/cm2 column air density to kg/m2
    airrho = np.divide(aircolumngrid,Av) # moles/cm2
    airrho = np.multiply(airrho,Mw_air)   # g/cm2
    airrho = np.divide(airrho,1e3)       # kg/cm2
    airrho = np.multiply(airrho,1e4)     # kg/m2

    co2rho = np.multiply(airrho,tmp)     # kg/m2
    co2rho = np.multiply(Mw_ratio,airrho)# kg CO2/m2
    
    # MEANWINDSPEED in m/s - used 1e3 to convert to km/s
    alpha1 = 1.0; alpha2 = 0.6
    U_eff = alpha1 * np.log10(MEANWINDSPEED) + alpha2

    #U_eff = MEANWINDSPEED

    #print(np.nansum(tmp)*3600*24*365/1e12,U_eff,Ldimension*1000)
    
    # L = dimension of domain
    IME = np.nansum(co2rho)*3600*24*365/1e12 * U_eff/(Ldimension*1000)

    return np.nansum(odiacco2flux)*3600*24*365/1e12, IME


def sourcepixelmethod(uwind,vwind,Ldimension,dxco2,surfacepressure):
    #convert xco2 to Q
    U = np.sqrt(uwind**2 + vwind**2)
    return U*Ldimension*dxco2/9.81
