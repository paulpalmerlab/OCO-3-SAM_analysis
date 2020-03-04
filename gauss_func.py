import numpy as np
import sys
from scipy.special import erfcinv as erfcinv
from calc_sigmas import calc_sigmas 
from haversine import haversine
from pylab import plt

def gauss_func(Q,u,dir1,x,y,xs,ys,STABILITY):

   #
   # Calculate wind speed and direction
   #
   # components of u in x and y directions
   wx=u*np.sin((dir1)*np.pi/180.);
   wy=u*np.cos((dir1)*np.pi/180.);

   #
   # Distance between grid points (expressed in degrees) and the stack in km
   # Also x and y components
   #
   hypotenuse = np.zeros((len(x),len(y)))
   xdist = np.zeros((len(x),len(y)))
   ydist = np.zeros((len(x),len(y)))
   for ii in np.arange(len(x)):
      for jj in np.arange(len(y)):
         hypotenuse[ii,jj] = haversine([xs,ys],[x[ii],y[jj]])*1e3
         xdist[ii,jj]      = haversine([xs,ys],[x[ii],ys])*1e3
         ydist[ii,jj]      = haversine([xs,ys],[xs,y[jj]])*1e3
         if (xs < x[ii]): xdist[ii,jj] = -xdist[ii,jj]
         if (ys < y[jj]): ydist[ii,jj] = -ydist[ii,jj]         
         
   # Need angle between point x, y and the wind direction, so use scalar product:
   dot_product=wx*xdist+wy*ydist;
   # product of magnitude of vectors:
   magnitudes=u*np.sqrt(xdist**2.+ydist**2.); 
   # angle between wind and point (x,y)
   subtended=np.arccos(dot_product/(magnitudes+1e-15));

   distfromsource = np.sqrt(xdist**2.+ydist**2.)
   
   # SANITY CHECK
   # distance to point x,y from stack
   #hypotenuse=np.sqrt(xdist**2.+ydist**2.);
   #
   #img = plt.pcolormesh(magnitudes)
   #plt.colorbar(img)
   #plt.show()
   #sys.exit()
   
   # Distance along the wind direction to perpendilcular line that intesects
   # x,y
   downwind=np.cos(subtended)*hypotenuse;

   # Distance cross wind.
   crosswind=np.sin(subtended)*hypotenuse;

   ind=np.where(downwind>0.);
   
   C=np.zeros((len(x),len(y)));

   # calculate sigmas based on stability and distance downwind
   (sig_y,sig_z)=calc_sigmas(STABILITY,downwind,distfromsource);

   C[ind]=Q/(2.*np.pi*u*sig_y[ind]) \
       * np.exp(-crosswind[ind]**2./(2.*sig_y[ind]**2.))
   
   return C

   
   
