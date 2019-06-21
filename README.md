# Coral
Modelling coral growth

The file `PS_bathy.csv` contains the bathymetric dataset. `PS_bathy_no_coral.csv` has the 9 coral areas masked out.

The script `coral_mask.py` pulls out the grid coordinates that define the coral areas into `coral_polygons.dat`.

The script `seafloor.py` estimates the sea floor depth using different nearest neighbour weightings.

The script `extrap_coral.py` extrapolates into the future or the past.

