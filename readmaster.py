import numpy  as np
import h5py
import pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
from matplotlib import ticker, cm
import pandas as pd
import os
from gauss_func import gauss_func
from xco2funcs import *
from haversine import haversine
from readfunc import *
import glob, os

# Isfahan and Arak, Iran
targets = [[51.6660,32.6539],[49.7013,34.0954]]
ntargets = len(targets)

#-------------------------
# ____  _____    _    ___`_    _____ ___ _     _____ ____  
#|  _ \| ____|  / \  |  _ \  |  ___|_ _| |   | ____/ ___| 
#| |_) |  _|   / _ \ | | | | | |_   | || |   |  _| \___ \ 
#|  _ <| |___ / ___ \| |_| | |  _|  | || |___| |___ ___) |
#|_| \_\_____/_/   \_\____/  |_|   |___|_____|_____|____/ 
#
#-------------------------

#InNO2FileName  = 'Isfahan_Iran/S5P_OFFL_L2__NO2____20191224T082923_20191224T101053_11378_01_010302_20191226T011059.SUB.nc'
#InGEOSFileName = 'Isfahan_Iran/oco3_L2MetSC_03630a_191224_B10106_200125115513.h5'
#InOCO3FileName = 'Isfahan_Iran/oco3_L2DiaSC_03630a_191224_B10106_200205125150s.h5'

# Isfahan Iran
City = 'oco3_L2_03630a_fossil0000/' ; useextent = [51,52,32,33]
# Buenos Aires, Argentina
City = 'oco3_L2_02194a_fossil0035/' ; useextent = [-60,-58,-36,-34]
# Buenos Aires, Argentina
City = 'oco3_L2_03962a_fossil0035/' ; useextent = [-60,-58,-36,-34]
# Delhi
City = 'oco3_L2_02450a_fossil0029/' ; useextent = [75,78,27.5,29.5]
## Tehran Iran 1
#City = 'oco3_L2_02482a_fossil0014/' ; useextent = [50,52,35,36.5]
## Tehran Iran 2
#City = 'oco3_L2_02543a_fossil0014/' ; useextent = [50,52,35,36.5]




InNO2FileName  = glob.glob(City+'S5P_OFFL*.nc')[0]
InGEOSFileName = glob.glob(City+'oco3_L2MetSC_*.h5')[0]
InOCO3FileName = glob.glob(City+'oco3_L2DiaSC_*.h5')[0]

#InNO2FileName  = 'Isfahan_Iran/S5P_OFFL_L2__NO2____20191224T082923_20191224T101053_11378_01_010302_20191226T011059.SUB.nc'
#InGEOSFileName = 'Isfahan_Iran/oco3_L2MetSC_03630a_191224_B10106_200125115513.h5'
#InOCO3FileName = 'Isfahan_Iran/oco3_L2DiaSC_03630a_191224_B10106_200205125150s.h5'

#-------------------------
# Read in coincident TROPOMI NO2
#-------------------------

qa, no2, latno2, lonno2 = readno2(InNO2FileName)

#-------------------------
# Read in coincident GEOS meteorology
#-------------------------


uwind, vwind, metlat, metlon = readmet(InGEOSFileName)

#-------------------------
# Read OCO-3 data
#-------------------------

xco2use, latxco2, lonxco2, obsmonth, spress_oco3 = readxco2(InOCO3FileName)

#------------------------- 
#  ____  _____ _____ ___ _   _ _____   ____   ___  __  __    _    ___ _   _ 
# |  _ \| ____|  ___|_ _| \ | | ____| |  _ \ / _ \|  \/  |  / \  |_ _| \ | |
# | | | |  _| | |_   | ||  \| |  _|   | | | | | | | |\/| | / _ \  | ||  \| |
# | |_| | |___|  _|  | || |\  | |___  | |_| | |_| | |  | |/ ___ \ | || |\  |
# |____/|_____|_|   |___|_| \_|_____| |____/ \___/|_|  |_/_/   \_\___|_| \_|
#                                                                           
#------------------------- 

extent      = [np.min(lonxco2), np.max(lonxco2), \
               np.min(latxco2), np.max(latxco2)]
extent      = useextent

#-------------------------
# Define common 1km ODIAC grid within defined domain
#-------------------------

ODIAClonsDOMAIN, ODIAClatsDOMAIN = definecommongrid(extent)

#-------------------------
# Move OCO-3 data on a regular 1km x 1km grid
#-------------------------

newgrid, nobs = regridoco3xco2(ODIAClonsDOMAIN, ODIAClatsDOMAIN,xco2use,latxco2,lonxco2)

spress_oco3_grid, nobs = regridoco3xco2(ODIAClonsDOMAIN, ODIAClatsDOMAIN,spress_oco3,latxco2,lonxco2)

#nobsind = np.where(nobs > 1)
#newgrid = np.divide(newgrid[nobsind],nobs[nobsind])

#-------------------------
#  ____  _____ _____ ___ _   _ _____   ____  _     ___ _____ 
# |  _ \| ____|  ___|_ _| \ | | ____| |  _ \| |   / _ \_   _|
# | | | |  _| | |_   | ||  \| |  _|   | |_) | |  | | | || |  
# | |_| | |___|  _|  | || |\  | |___  |  __/| |__| |_| || |  
# |____/|_____|_|   |___|_| \_|_____| |_|   |_____\___/ |_|  
# 
#-------------------------

central_lat = 0
central_lon = 0
    
proj = ccrs.PlateCarree(central_longitude=central_lon)

fig1, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3, sharex=True, \
                                           sharey=True, figsize=(12,8),\
                                           subplot_kw={'projection': proj})

#-------------------------
# PLOT OCO-3 data
#-------------------------

# NOTE: see readxco2.py
# AM = area mode; TG = target
# sam_or_tg = 'AM'

#
# Remove bottom 10th percentile of gridded OCO-3 data to identify city signal
#
tmpind   = np.where(newgrid>0); nonzeroxco2 = newgrid[tmpind]
bttenpc  = np.nanpercentile(nonzeroxco2,10)
newgrid  = np.subtract(newgrid,bttenpc)

v_min_oco3 = np.nanmedian([newgrid]) - np.nanstd([newgrid])
v_max_oco3 = np.nanmedian([newgrid]) + np.nanstd([newgrid])

colors = newgrid

newgrid = np.transpose(newgrid)

uselabel = 'OCO-3 $\Delta$XCO$_2$ (ppm)'
plotdata(ax1,fig1,ODIAClonsDOMAIN,ODIAClatsDOMAIN,extent,newgrid,\
         colors,v_min_oco3,v_max_oco3,uselabel,usecontour=1,no2flag=0)    
ax1.set_title(InOCO3FileName)

ax1.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
ax1.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none')

#-------------------------
# PLOT CORRESPONDING GC winds
#-------------------------
        
uwind = np.squeeze(uwind[:,:,-1]); vwind = np.squeeze(vwind[:,:,-1])

#Reshape the n_soundings x footprint, for convenience
metlat = np.reshape(metlat,-1)
metlon = np.reshape(metlon,-1)    
uwind  = np.reshape(uwind,-1)
vwind  = np.reshape(vwind,-1)

indlon = np.where((metlon >= extent[0]) & (metlon <= extent[1]))
indlat = np.where((metlat >= extent[2]) & (metlat <= extent[3]))
uwind = uwind[indlon]; vwind = vwind[indlat]; metlat = metlat[indlat]; metlon = metlon[indlon]

newuwindgrid, windnobs = regridmet(ODIAClonsDOMAIN,uwind,metlon)
newvwindgrid, windnobs = regridmet(ODIAClatsDOMAIN,vwind,metlat)

X, Y = np.meshgrid(ODIAClonsDOMAIN,ODIAClatsDOMAIN)
UWIND, VWIND = np.meshgrid(newvwindgrid,newuwindgrid)

UWIND = np.transpose(UWIND); VWIND = np.transpose(VWIND)

print(np.shape(newuwindgrid),np.shape(newvwindgrid))
print(np.shape(newgrid))
print(np.shape(UWIND))
#sys.exit()

for ii in np.arange(len(ODIAClatsDOMAIN)):
    for jj in np.arange(len(ODIAClonsDOMAIN)):
        if np.isnan(newgrid[ii,jj]):
            UWIND[ii,jj] = np.nan            
            VWIND[ii,jj] = np.nan                        


#for ii in np.arange(len(newuwindgrid)):
#    for jj in np.arange(len(newvwindgrid)):
#        if (np.isnan(newgrid[jj,ii])):
#            UWIND[ii,jj] = np.nan
#            VWIND[ii,jj] = np.nan


ax3.quiver(X,Y,UWIND,VWIND,\
           scale=5,scale_units='inches',color='black',alpha=0.5)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False

ax1.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
ax1.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none')



winddir_target = []
windspd_target = []


# This is for the plume model UNCOMMENT
#for ii in np.arange(ntargets):
#    xdist = abs(np.subtract(metlon,targets[ii][0]))
#    ydist = abs(np.subtract(metlat,targets[ii][1]))
#    xminval = np.min(xdist); xminind = np.where(xdist==xminval); uselon = metlon[xminind]
#    yminval = np.min(ydist); yminind = np.where(ydist==yminval); uselat = metlat[yminind]        
#    uwind_use = uwind[xminind]; vwind_use = vwind[yminind]
#    wdir      = (270-np.rad2deg(np.arctan2(vwind_use,uwind_use)))%360
#    wspd      = np.sqrt(uwind_use**2 + vwind_use**2)
#    winddir_target.append(wdir)
#    windspd_target.append(wspd)

#------------------------- 
# PLOT TROPOMI NO2
#-------------------------    

no2, latno2, lonno2 = filterandconvertno2(no2,latno2,lonno2,qa)

v_min = np.min(no2)
v_max = np.max(no2)
colors = no2

uselabel = 'TROPOMI NO$_2$ (10$^{16}$ molec/cm$^2$)'
plotdata(ax2,fig1,lonno2,latno2,extent,no2,colors,v_min,v_max,uselabel,usecontour=0,no2flag=1)
ax2.set_title(InNO2FileName)

ax2.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
ax2.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none')

#-------------------------
# Read ODIAC data
#-------------------------

odiac_outfile = 'odiac.'+'{}'.format(obsmonth)+'.npz'
#
#tmp = ['01','02','03','04','05','06',\
#       '07','08','09','10','11','12']
#
## Saves out global field - does using npz file save time?
## Unit = Tonnes of carbon/1km2/month
#for ii in np.arange(len(tmp)):
#    obsmonth = tmp[ii]
#    odiac_outfile = 'odiac.'+obsmonth+'.npz'
#    print(odiac_outfile)
#    co2flux, odiaclons, odiaclats = readodiac(obsmonth)    
#    np.savez(odiac_outfile,co2flux=co2flux,odiaclons=odiaclons,odiaclats=odiaclats)
#sys.exit()

npzfile       = np.load(odiac_outfile)
co2flux       = npzfile['co2flux']
odiaclons     = npzfile['odiaclons']
odiaclats     = npzfile['odiaclats']

#print(np.shape(co2flux))
#co2flux = np.flip(co2flux,axis=1) # flip lats
print(np.shape(co2flux))
print(np.shape(odiaclons),np.shape(odiaclats))
#sys.exit()

# emission units arre no gCO2/km2/s
odiacco2flux = convertodiacemissionunits(co2flux,int(obsmonth))

indodiaclon  = np.where((odiaclons >= extent[0]) & (odiaclons <= extent[1]))
indodiaclat  = np.where((odiaclats >= extent[2]) & (odiaclats <= extent[3]))

odiaclats    = np.squeeze(odiaclats[indodiaclat])
odiaclons    = np.squeeze(odiaclons[indodiaclon])

odiacco2flux = np.squeeze(odiacco2flux[indodiaclat,:])
odiacco2flux = np.squeeze(odiacco2flux[:,indodiaclon])


#print(odiaclats[0:9])
#print(ODIAClatsDOMAIN[0:9])
#print(odiaclons[0:9])
#print(ODIAClonsDOMAIN[0:9])
#print(np.shape(odiaclats),np.shape(ODIAClatsDOMAIN))
#print(np.shape(odiaclons),np.shape(ODIAClonsDOMAIN))

#-------------------------
# Plot ODIAC emission data
#-------------------------

print(np.shape(odiaclons),np.shape(odiaclats))
print(np.shape(newgrid),np.shape(odiacco2flux))
#sys.exit()


for ii in np.arange(len(odiaclats)):
    for jj in np.arange(len(odiaclons)):
        if np.isnan(newgrid[ii,jj]): odiacco2flux[ii,jj] = np.nan





        
X,Y = np.meshgrid(odiaclons,odiaclats)

odiacdatapoints = np.count_nonzero(~np.isnan(co2flux))
if odiacdatapoints > 0:

    #ind = np.where((np.isnan(np.transpose(newgrid))))
    #odiacco2flux[ind] = np.nan

    from matplotlib.colors import LogNorm
    
    surfemiss = ax4.pcolor(X,Y,odiacco2flux,\
                               transform=ccrs.PlateCarree(),\
                               norm=LogNorm(vmin=np.nanmin(co2flux),vmax=np.nanmax(co2flux)),\
                               vmin=10,vmax=1e3,\
                               cmap=plt.cm.get_cmap('rainbow'))
    
    gl2 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')
    
    cb2 = fig1.colorbar(surfemiss,ax=ax4,extend='both',orientation = 'horizontal',\
                        label = 'ODIAC CO2 emissions [gC/km2/s]')

    gl2.xlabels_top = False
    gl2.ylabels_right = False
        
    ax4.coastlines(resolution='50m', color='black', linewidth=1)

    ax4.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
    ax4.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none')

#-------------------------
# Plot ODIAC XCO2 after making some assumptions
#-------------------------

    useboundarylayerheight = 1 # km

    odiac_xco2 = get_xco2_from_emission(odiacco2flux,int(obsmonth),useboundarylayerheight)
    
    odiac_xco2 = np.divide(odiac_xco2,1e-6)

    print(np.nanmin(odiac_xco2),np.nanmax(odiac_xco2))
    print(np.nanmedian(odiac_xco2))
    #sys.exit()
    
    odiacQ_target = []

    for ii in np.arange(ntargets):
        xdist = abs(np.subtract(odiaclons,targets[ii][0]))
        ydist = abs(np.subtract(odiaclats,targets[ii][1]))
        xminval = np.min(xdist); xminind = np.where(xdist==xminval); uselon = odiaclons[xminind]
        yminval = np.min(ydist); yminind = np.where(ydist==yminval); uselat = odiaclats[yminind]
        xminind = np.squeeze(xminind); yminind = np.squeeze(yminind)
        Q_use   = np.nansum(co2flux[xminind-2:xminind+2,yminind-2:yminind+2])
        Q_use_unit = convertodiacemissionunits(Q_use,int(obsmonth))
        odiacQ_target.append(Q_use_unit)

    print('----ODIAC Q values (g CO2/km box/sec) ---')
    print(odiacQ_target)
    print('----')
        
    v_min = 0
    v_max = 3#np.nanmax(odiac_xco2)
    npts  = 100
    uselevels = np.arange(npts)*(v_max-v_min)/npts        


    surfxco2 = ax5.pcolormesh(X,Y,odiac_xco2,\
                              transform=ccrs.PlateCarree(),\
                              vmin=v_min_oco3,vmax=v_max_oco3,\
                              cmap=plt.cm.get_cmap('rainbow'))
    
    gl2 = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')
    
    cb3 = fig1.colorbar(surfxco2,ax=ax5,extend='both',orientation = 'horizontal',\
                        label = 'ODIAC XCO2 [ppm]')
    cb3.set_clim(v_min_oco3,v_max_oco3)
    
    gl2.xlabels_top = False
    gl2.ylabels_right = False
    
    ax5.coastlines(resolution='50m', color='black', linewidth=1)

    ax5.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
    ax5.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none') 



# Saves out field to test source quantification methods
#test_outfile = 'sourceestimationtestfile.npz'
#np.savez(test_outfile,oco3co2=newgrid,spressoco3=spress_oco3_grid,\
#         odiacco2emission=odiacco2flux,odiacxco2=odiac_xco2,\
#         uwind=UWIND,vwind=VWIND,\
#         lons=ODIAClonsDOMAIN,lats=ODIAClatsDOMAIN)
#sys.exit()
    

        
#-------------------------             
#  ___ __  __ _____                  _   _               _ 
# |_ _|  \/  | ____|  _ __ ___   ___| |_| |__   ___   __| |
#  | || |\/| |  _|   | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` |
#  | || |  | | |___  | | | | | |  __/ |_| | | | (_) | (_| |
# |___|_|  |_|_____| |_| |_| |_|\___|\__|_| |_|\___/ \__,_|
#
# Integrated mass enhancement
#-------------------------             

meanwindspeed = 5    # m/s
Ldimension    = 10 # km

print(np.shape(odiacco2flux),np.shape(newgrid),np.shape(xco2use))

ODIACEMISSION, OCO3_IME = calc_ime(odiacco2flux,newgrid,meanwindspeed,Ldimension)

imestring   = writescreenstring(OCO3_IME,0)
odiacstring = writescreenstring(ODIACEMISSION,0)
ypos = np.mean(useextent[2:3]); xpos = np.mean(useextent[0:1])
dy = np.abs(useextent[2]-useextent[3]); dx = np.abs(useextent[1]-useextent[0])
ax1.annotate('OCO-3 IME: '+imestring,(xpos,ypos+0.25*dy),horizontalalignment='left',fontsize=8, color='black')
ax1.annotate('ODIAC: '+odiacstring,(xpos,ypos+0.15*dy),horizontalalignment='left',fontsize=8, color='black')

#-------------------------             
# Ad hoc emission estimate for urban core
#-------------------------             

model = np.reshape(odiacco2flux,-1)
oco3  = np.reshape(newgrid,-1)
ind = np.where(~np.isnan(oco3) & ~np.isnan(model))
oco3 = oco3[ind]; model = model[ind]

r = np.corrcoef(model,oco3)
print(r)
print(np.shape(model),np.shape(oco3))

#print(model)
#print(oco3)
#plt.figure(20)
#plt.plot(oco3,model,'ko')




#-------------------------             
#  ____  _   _ _   _   ____  _    _   _ __  __ _____ 
# |  _ \| | | | \ | | |  _ \| |  | | | |  \/  | ____|
# | |_) | | | |  \| | | |_) | |  | | | | |\/| |  _|  
# |  _ <| |_| | |\  | |  __/| |__| |_| | |  | | |___ 
# |_| \_\\___/|_| \_| |_|   |_____\___/|_|  |_|_____|
#                                                    
#  __  __  ___  ____  _____ _     
# |  \/  |/ _ \|  _ \| ____| |    
# | |\/| | | | | | | |  _| | |    
# | |  | | |_| | |_| | |___| |___ 
# |_|  |_|\___/|____/|_____|_____|
#                                 
#-------------------------      

#-------------------------
# Run Gaussian plume model
#-------------------------        

# Stability
CONSTANT_STABILITY=1
stab1 = 1
stability = stab1
stability_str=['Very unstable','Moderately unstable','Slightly unstable', \
               'Neutral','Moderately stable','Very stable'];    
# Define model grid
#latres = lonres = 0.05
#mlon = np.mgrid[extent[0]:extent[1]:lonres]
#mlat = np.mgrid[extent[2]:extent[3]:latres]
#x = mlon; y = mlat

#
# Use lat/lon from 1-km ODIAC
#
x = ODIAClonsDOMAIN; y = ODIAClatsDOMAIN

# Number of stacks
stacks = 2
# Define stack locations
stack_x=[targets[0][0], targets[1][0], -500.];
stack_y=[targets[0][1], targets[1][1], -200.];
# Define stacks emissions
Q=[6., 8., 40.]
# Define wind directions and speeds
#wind_dir   = [0,23,0]
#wind_speed = [10,5,12]    
#
#**Eventually need to include injection height
#


#
# Initialise model field
#
# UNCOMMENT
#C1=np.zeros((len(x),len(y)))
#for j in range(0,stacks):
#    C=np.ones((len(x),len(y)))
#    
#    #C=gauss_func(odiacQ_target[j],windspd_target[j],winddir_target[j],x,y,
#    C=gauss_func(Q[j],windspd_target[j],winddir_target[j],x,y,                     
#                 stack_x[j],stack_y[j],stability);
#    C1[:,:]=C1[:,:]+C;

#
# Plot output
#

[y,x]=np.meshgrid(y,x); # x and y defined at all positions on the grid

# use levels from ODIAC
#uselevels = np.arange(25)*0.01

v_min = 0
v_max = 0.5#np.nanmax(odiac_xco2)
npts  = 50
uselevels = np.arange(npts)*(v_max-v_min)/npts        

# C1 units should be g/m2 [g/s / m/s * m]
# We need XCO2
Av      = 6.023e23         # Avagrado's number molecules/mole
rho_air = 2.69e25          # molec/m3
#C1 = np.divide(C1,44)      # mole  CO2/m3
#C1 = np.multiply(C1,Av)    # molec CO2/m3
#C1 = np.divide(C1,rho_air) # mole fraction
#C1 = np.divide(C1,1e-6)    # ppm
# UNCOMMENT
#print(np.min(C1),np.max(C1))


#C1ind = np.where(np.isnan(newgrid))
#C1[C1ind] = np.nan


#plotdata(ax4,fig1,x,y,extent,C1,uselevels,v_min,v_max,'tmp',usecontour=0,no2flag=0)    

#levels=uselevels,\

# UNCOMMENT
#modelxco2 = ax6.pcolor(x,y,C1,\
#                        transform=ccrs.PlateCarree(),\
#                        cmap=plt.cm.get_cmap('rainbow'))
#
#gl3 = ax6.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                    linewidth=2, color='gray', alpha=0.5, linestyle='--')
#
#cb3 = fig1.colorbar(modelxco2,ax=ax6,extend='both',orientation = 'horizontal',\
#                    label = 'Gaussian plume model XCO2 [ppm]')
#cb3.set_clim(v_min,v_max)
#
#gl3.xlabels_top = False
#gl3.ylabels_right = False
#

ax6.plot(targets[0][0],targets[0][1],'ko',markersize=10,markerfacecolor='none')
ax6.plot(targets[1][0],targets[1][1],'ko',markersize=10,markerfacecolor='none')

    
    
plt.show()
sys.exit()




#    
plt.show()






