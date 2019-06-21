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


from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import box

import ast
from plot_csv import read_bathy_csv
from coral_mask import process_poly, read_bathy
from seafloor import get_coral_polys


if len(sys.argv) < 8:
	print("usage: python3 extrap_coral.py seed err sea_rise nn_floor reverse Tmax Ti")
	print("seed = random seed (integer)")
	print("err = Central, Gaussian or BadYear.")
	print("sea_rise = sea_level rise in mm/yr")
	print("nn_floor. extrapolation of seafloor.") 
	print("if reverse == reverse: extrapolate backwards. default forwards");
	print("Tmax = how far to extrapolate");
	print("Ti = how often to write out data");
	
##Parameters
my_buffer = 25;
lateral_growth = 0.25;
dist_size = 1000; #for the BadYear setting

seed = int(sys.argv[1]);	
np.random.seed(seed=(seed+123456789))
err = sys.argv[2]
sea_rise_py = float(sys.argv[3]);
nn_floor = sys.argv[4]
reverse = False;
if sys.argv[5].lower() == "reverse":
	reverse = True;
Tmax = int(sys.argv[6]);
Ti = int(sys.argv[7]);


growth = [ 
[ [-4000,-3500], 1.4, 0, 'k' ],
[ [-3500,-3000], 5.5,  5.2, 'r'],
[ [-3000,-2500], 6.4 , 7.8, 'g' ],
[ [-2500,-2000], 6.9 , 9.4, 'b' ],
[ [-2000,-1500], 5.1 , 4.0, 'y' ],
[ [-1500,-1000], 2.5 , 1.3, 'c' ],
[ [-1000,-500], 3.8 , 3.8, 'm' ],
[ [-500,0], 3.4 , 5.1, 'seagreen'],
[ [0,1500], 1.4 , 1.1, 'peru']
]

##read bathy grid
polygons, polygon_coords = get_coral_polys(my_buffer=my_buffer);
x, y, ht = read_bathy_csv('PS_bathy.csv')

##coral true/false array
mk = np.empty( ht.shape );
mj = len(x)
mi = len(y)
for i in range( mi ):
	for j in range( mj ):
		got = False;
		pt = Point( x[j], y[i] );
		for r,p in enumerate(polygons):
			if p.contains( pt ):
				mk[i][j] = 1
				got = True;
				break;
		if not got:
			mk[i][j] = 0;
sizex, sizey = mk.shape
				
## Get the sea floor
tx, ty, seafloor = read_bathy_csv('seafloor_extrap_' + str(my_buffer) + '_' + str(nn_floor) + '.csv');

##average grid spacing (x_i - x_i-1) should be same for all i anyway.
step = 0;
for i in range(1,len(x)):
	step += x[i] - x[i-1];
step /= (len(x) - 1);

lg_step = step/lateral_growth;
lg_add = lateral_growth/step;


if err == "BadYear":
	rands = [];
	for g in growth:
		r = norm.rvs(loc=g[1], scale=g[2], size=dist_size)
		r.sort();
		rands.append(r);


if reverse:
	filename=err +str(seed) + "_r" + str(sea_rise_py)  + "_" + str(my_buffer) + '_' + str(nn_floor) + "_backward"; 	
	tag = "-"
	tdir = -1
else:
	filename=err +str(seed) + "_r" + str(sea_rise_py)  + "_" + str(my_buffer) + '_' + str(nn_floor) + "_forward"; 
	tag = "+"
	tdir = 1;

	
numc = None;
num = 0
cl = 0.999;
##dealing with a few edge cases where seafloor extrapolation introduces anomalies.
if reverse:
	for i in range(sizex):
		for j in range(sizey):	
			if mk[i][j] > cl: #is coral
				if seafloor[i][j] > ht[i][j]: ##seafloor is higher than current coral height! 
					ht[i][j] = seafloor[i][j]; #can only go down to the seafloor
					mk[i][j] = 0; #land not coral
					


with open(filename + "_t" + tag + str(0) + ".csv", 'w') as ofile:
	for i in range(sizex):
		for j in range(sizey):
			ofile.write("{}, {}, {}, {}\n".format( x[j], y[i], ht[i][j]/1000., mk[i][j]) )



T = Tmax;
skip=Ti

for t in range(T):
	if err == "BadYear":
		year_idx = np.random.random_integers(0, high=(dist_size-1))

	num = 0;
	##vertical growth
	for i in range(sizex):
		for j in range(sizey):	
			if mk[i][j] >= cl:
				if (not reverse) or (ht[i][j] > seafloor[i][j] and ht[i][j] > -4000): ##stop going back in time when you hit the floor
										
					for k,g in enumerate(growth):
						if ht[i][j] >= g[0][0] and ht[i][j] <= g[0][1]: ##find appropriate growth rate

							if err == "Central":
								ht[i][j] += g[1]*tdir;	
							elif err == "Gaussian":
								ht[i][j] += norm.rvs(loc=g[1], scale=g[2])*tdir 
							else:
								ht[i][j] += rands[k][year_idx]*tdir; #BadYear
						
							num += 1;						
							break;
			if reverse:
				if ht[i][j] < seafloor[i][j]: #clip overshooting
					ht[i][j] = seafloor[i][j];
					mk[i][j] = 0;
				if ht[i][j] < -4000:
					mk[i][j] = 0;
						
	if numc is None: numc = float( num ); ##how many grid cells changed
	
	##horizontal growth	
	if not reverse:
		##lateral growth
		new_coral = np.zeros( mk.shape );
		for i in range(sizex):
			for j in range(sizey):	
				#not coral
				if mk[i][j] < cl:
					# neighbour is coral
					if i+1 < sizex and mk[i+1][j] >= cl: new_coral[i][j] += lg_add;
					if i-1 >= 0    and mk[i-1][j] >= cl: new_coral[i][j] += lg_add;
					if j+1 < sizey and mk[i][j+1] >= cl: new_coral[i][j] += lg_add;
					if j-1 >= 0    and mk[i][j-1] >= cl: new_coral[i][j] += lg_add;

		for i in range(sizex):
			for j in range(sizey):	
				mk[i][j] += new_coral[i][j];
		
		##sea_level rise		
		for i in range(sizex):
			for j in range(sizey):	
				ht[i][j] -= sea_rise_py*tdir;

						
	if (t+1) % skip == 0:
		print(t+1);
		with open(filename + "_t" + tag + str(t+1) + ".csv", 'w') as ofile:
			for i in range(sizex):
				for j in range(sizey):			
					ofile.write("{}, {}, {}, {}\n".format( x[j], y[i], ht[i][j]/1000., mk[i][j]) )

