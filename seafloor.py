from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib.mlab import griddata
import scipy.interpolate as interpolate
import sys
import random
import math
from scipy.stats import norm
from scipy.stats import truncnorm
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import griddata


from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import box

import ast
from plot_csv import read_bathy_csv
from coral_mask import process_poly, read_bathy

def findall_nearest(array,x, y, test):
	"""
	sort points by proximity to test
	"""
	min_dist = np.inf;
	idx = -1;
	dists = []
	for i,a in enumerate(array):
		v = ( ( y[a[0]]-y[test[0]])**2+(x[a[1]]-x[test[1]])**2 )
		dists.append( [ i, 1.0/math.sqrt(v) ] );
	dists.sort(reverse=True, key=lambda x:x[1]);	
	return np.array( dists );
			

def get_coral_polys(filename='coral_polygons.dat', my_buffer=25):
	"""
	##Grid of where the flipping coral is##
	add a buffer around the coral polygons to ensure we get everything inside
	"""

	polygons = []
	polygon_coords = []
	with open(filename, 'r') as infile:
		for line in infile:
			p = ast.literal_eval( line );
			poly = Polygon( p ).buffer(my_buffer);
			polygons.append( poly );
			
			pols, ipols = process_poly( poly );
			polygon_coords.append( pols );
			
	return polygons, polygon_coords

if __name__ == "__main__":

		
	my_buffer = 25;
	polygons, polygon_coords = get_coral_polys(my_buffer=25);
	x, y, ht, mk = read_bathy('PS_bathy.csv', 'PS_bathy_no_coral.csv');
		

	print( "Finding edges..." )
	edges = [ [] for p in polygons ];
	edge_map = []
	mj = len(x)
	mi = len(y)
	for i in range( mi ):
		tmp = []
		for j in range( mj ):
			inp = False;
			
			if i>0 and j>0 and i<mi-1 and j<mj-1:
				for r,p in enumerate(polygons):
					pt = Point(x[j], y[i]);
					if p.contains(pt): break;
					
					ptu = Point(x[j], y[i+1]);
					ptd = Point(x[j], y[i-1]);
					ptl = Point(x[j-1], y[i]);
					ptr = Point(x[j+1], y[i]);
					if (p.contains(ptu) or p.contains(ptd) or p.contains(ptl) or p.contains(ptr)):
						if not inp: tmp.append( 2000 );
						edges[r].append( (i,j) );
						inp = True;
						#break;
					
			if not inp:
				tmp.append( 0 );
		edge_map.append(tmp)
	#edge_map = np.array( edge_map )
	#CS = plt.contourf(x,y,edge_map, 60, cmap=plt.cm.ocean_r)
	#plt.colorbar()  
	#plt.show()
	#plt.close();


	print( "Finding floor..." )
	seafloor = [];
	nearest = [1,2,4,6,8,10,20,50,100,150,200,250,500];
	mj = len(x)
	mi = len(y)
	for i in range( mi ):
		print( mi-i )
		tmp = [];
		for j in range( mj ):
			est = []
			inp = False;
			pt = Point(x[j], y[i])

			for r,p in enumerate(polygons):
				
				if p.contains( pt ): ##in the coral patch
	 
					weights = findall_nearest(edges[r], x, y ,( i,j ) ); #all points close to i,j
					for n in nearest:
						av = 0;
						norm = 0;
						for k in range( min(n, len(weights)) ): #weighted average
							idx = int(weights[k][0]);
							av += weights[k][1] * ( ht[ edges[r][idx][0] ][ edges[r][idx][1] ] );
							norm += weights[k][1];
						av = av/norm;
						est.append( av );
					inp = True;
					break;

			if not inp:
				for n in nearest:
					est.append( ht[i][j] );
			
			tmp.append( est )
		seafloor.append(tmp)
	seafloor = np.array( seafloor );		

	for k in range(len(nearest)):
		print("plot", nearest[k])
		CS = plt.contourf(x,y,seafloor[:,:,k], 60, cmap=plt.cm.ocean_r)
		plt.colorbar()  
		plt.savefig("seafloor_extrap_" + str(my_buffer) + "_" + str(nearest[k]) + ".png");
		plt.close();

		with open("seafloor_extrap_" + str(my_buffer) + "_" + str(nearest[k]) + ".csv", 'w') as outfile:
			for i in range( mi ):
				for j in range( mj ):
					outfile.write("{},{},{}\n".format( x[j], y[i], seafloor[i][j][k]/1000.) )
