import csv
import numpy as np
import matplotlib.pyplot as plt
import copy

from shapely.geometry import Polygon
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
from matplotlib.patches import Polygon as pgn

import sys

def read_bathy(filename, mask_filename):
	"""
		## The CSV grid is not exactly uniform, but it is close enough ##
		## Record heights and coral presence
		## filename: csv file with bathymetry data
		## mask_filename: csv with masked coral positions
	"""

	coral = [];
	with open('PS_bathy_no_coral.csv', 'r') as csvfile:
		spamreader = csv.reader(csvfile, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
		for row in spamreader:	
			if row[2] == '1.70E+38':
				x = float(row[0]);
				y = float(row[1]);
				if x > 451920 and x < 454520 and y > 7885390 and y < 7888970: #in the sea
					if y < 7885670 and x < 452182: continue; #corner piece
					coral.append( (row[0],row[1]) );
				
	x = [];
	y = [];
	mk = [];
	ht = [];

	with open(filename, 'r') as csvfile:
		spamreader = csv.reader(csvfile, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
		line = [];
		cline = [];
		xlast = np.inf;
		ylast = np.inf;
		
		gotx = False;
		goty = False;
	
		for row in spamreader:
			
			xcoord = float(row[0]);
			ycoord = float(row[1]);
			
			if row[2] == '1.70E+38':
				zcoord = 0; #max height in data set
			else:
				zcoord = float(row[2])*1000; #convert to mm
			
			if (row[0],row[1]) in coral:
				iscoral = 1.0;
			else:
				iscoral = 0.0;
			
			if len(line) == 0:
				y.append( ycoord );
				
			if len(line) > 0 and xcoord < xlast:
				ht.append( line );
				mk.append( cline );
				line = [ zcoord ];
				cline = [ iscoral ];
				gotx = True;
				y.append( ycoord );
				
			else:
				line.append( zcoord );
				cline.append( iscoral );
				if not gotx: x.append( xcoord );
				
			xlast = xcoord
			ylast = ycoord
		
	ht.append( line )
	mk.append( cline )

	##Convert to np arrays
	ht = np.array( ht );
	mk = np.array( mk );
	x = np.array(x)
	y = np.array(y)
	return x, y, ht, mk;

def find_nearest(array, x, y, test):
	"""
	find the point in 'array' closest to 'test'
	"""
	min_dist = np.inf;
	idx = -1;
	for i,a in enumerate(array):
		v = ( ( y[a[0]]-y[test[0]])**2+(x[a[1]]-x[test[1]])**2 )
		if v < min_dist:
			min_dist = v;
			idx = i;		
	return idx, min_dist;
		
def process_poly(poly):
	"""
	Turn shapely polygons into things lists
	"""
	lons, lats = poly.exterior.coords.xy
	x,y=(lons, lats);
	pols = list(zip(x,y))
	ipols = [];
	for i in poly.interiors:
		lons, lats = i.coords.xy;
		x,y=(lons, lats);
		ipols.append( list(zip(x,y)) )
			
	return pols, ipols
	

if __name__ == "__main__":

	x, y, ht, mk = read_bathy('PS_bathy.csv', 'PS_bathy_no_coral.csv');


	##Find the grid positions which do not contain coral but border some coral
	edges = [];
	mj = len(x)
	mi = len(y)
	for i in range( 1,mi-1 ):
		for j in range( 1,mj-1 ):
			if mk[i][j] == 0 and ( mk[i][j+1] == 1 or mk[i][j-1] == 1 or mk[i+1][j] == 1 or mk[i-1][j] == 1): #diamond neighbourhood 
			#and ( mk[i+1][j+1] == 1 or mk[i+1][j-1] == 1 or mk[i-1][j+1] == 1 or mk[i-1][j-1] == 1): #moore neighbourhood
				edges.append( (i,j) );
	edge_copy = copy.deepcopy(edges)


	poly = []; #list of coral polygons
	pt = edge_copy[0]; 
	poly.append( ( x[pt[1]], y[pt[0]]) ); #insert the first edge point
	del edge_copy[0]

	poly_group = []
	poly_list_group = [];
	while len(edge_copy) > 0: ##look at every border point in the coral mask 
		idx, dist = find_nearest(edge_copy, x, y , pt ); ##find the closest border point to the current position

		fin_poly = False;
		if abs(pt[0] -  edge_copy[idx][0]) > 1 or abs(pt[1] -  edge_copy[idx][1]) > 1: ##check if the polygon is closed
			fin_poly = True;

		if fin_poly:
			if len(poly) > 2: ##we at least have a triangle!
				poly_list_group.append( poly ) ##add to co-ord list
				pl = Polygon( poly )			##convert to shapely polygon
				poly_group.append( pl );		##add to polygon list
			poly = [];
					
		pt = edge_copy[idx]; ##next point
		poly.append( ( x[pt[1]], y[pt[0]]) );
		del edge_copy[idx]
		
	if len(poly) > 2: ##finish off last polygon
		poly_list_group.append( poly )
		pl = Polygon( poly )
		poly_group.append( pl );
						
	print(len(edge_copy))# should be 1;

	##Plot the polygons
	fig, ax = plt.subplots()
	for pl in poly_group:
		pols, ipols = process_poly( pl );
		county = pgn(pols, fc='w', ec=(0,0,0,1), lw=1, alpha=1, fill=None, zorder=100)
		ax.add_artist(county)
	##drawing in the edges
	edge_map = [];
	for i in range( mi ):
		tmp = [];
		for j in range( mj ):
			if (i,j) in edges:
				tmp.append(1);
			else:
				tmp.append(0)
		edge_map.append(tmp)
	##contour plot
	CS = plt.contourf(x,y,edge_map,1,zorder=0);

	plt.colorbar()  
	plt.xlim([451920,454520 ])
	plt.ylim([7885390,7888970])
	plt.show()


	with open( "coral_polygons.dat", 'w' ) as outfile:
		for p in poly_list_group:
			outfile.write( str(p) + "\n" )
			
			
