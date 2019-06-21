import matplotlib.pyplot as plt
import csv
import sys
import numpy as np

def read_bathy_csv(filename):

	x = [];
	y = [];
	z = [];
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
				zcoord = float(row[2])*1000; #convert to millimeters
		
			if len(line) == 0:	y.append( ycoord );
				
			if len(line) > 0 and xcoord < xlast: #we finished reading left to right and got to a new line in the grid
				ht.append( line );
				line = [ zcoord ];
				y.append( ycoord );
				gotx = True;
			else: #only get in here on first line where we read the x-coordinates
				line.append( zcoord );
				if not gotx: x.append( xcoord );
				
			xlast = xcoord
			ylast = ycoord
	ht.append( line )

	ht = np.array( ht );
	x = np.array(x)
	y = np.array(y)

	return x, y, ht;

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("usage: python3 plot_csv.py bathymetry.csv")
		print("output: plot of bathymetry grid");

	n_contour = 60;
	x, y, ht = read_bathy_csv(sys.argv[1]);
	print(len(x), len(y), len(ht))

	CS = plt.contourf(x,y,ht, n_contour, cmap=plt.cm.ocean_r)
	plt.colorbar()  
	plt.savefig(sys.argv[1].replace("csv", "png"))
	plt.show()


