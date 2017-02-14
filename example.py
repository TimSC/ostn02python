from __future__ import print_function
from OSGB import *
from OSTN02 import *
import transform

#OSTN02 for Python
#=================

#See README.md for usage.
#See COPYING for redistribution terms

if __name__=="__main__":

	print("OSGB36_to_ETRS89")
	print("================")
	print("Take OS map reference: TR143599")

	xin, yin = parse_grid("TR", 14300, 59900)

	print("OS X (Eastings) "+str(xin))
	print("OS Y (Northings) "+str(yin))

	(x,y,h) = OSGB36_to_ETRS89 (xin, yin)

	print("Using the OSGB36_to_ETRS89 conversion gives us the grid position:")
	print(str(x) + "," + str(y) + "," + str(h))

	(gla, glo) = grid_to_ll(x, y)

	print("The grid position converts to ETRS89 lat,lon (using grid_to_ll) of:")
	print(str((gla, glo)))

	print("Actual answer: 51.297880, 1.072628")

	print("ETRS89 is within a metre of WGS84 (as used by GPS receivers), at time of writing (2011).")

	print("\nGridRef to ETRS89")
	print("==================")

	ref = "TR143599"
	print("Starting point: "+ref)
	print(transform.OSGB36GridRefToETRS89(ref))
	print("Expected:", 51.297880, 1.072628)

	print("\nETRS89_to_OSGB36")
	print("=================")

	gla = 51.297880
	glo = 1.072628
	h = 44.621

	print("To ETRS89 grid (using ll_to_grid):")
	(x2,y2) = ll_to_grid(gla, glo)
	print(str((x2,y2)))

	print("To OS Eastings/Northings (using ETRS89_to_OSGB36):")
	print(ETRS89_to_OSGB36(x2,y2,h))

	print("Actual Answer: 614300, 159900, 0")

	print("\nExceptions")
	print("==========")
	print("Some areas do not have OSTN02 coverage for example OS grid 622129,185038")
	print("This causes an exception to be raised, but it can be caught:")
	try:
		(x,y,h) = OSGB36_to_ETRS89 (622129,185038)
	except Exception as e:
		print('Exception occurred, value:', e)

	print("\nAn approximate transform is available where OSTN02 is not defined.")
	lat, lon = 51.520557, 1.200446
	print("Convert {}, {}".format(lat, lon))
	(lat2, lon2, alt2) = shift_ll_from_WGS84(lat, lon, 0.0)
	print ("Result:", ll_to_grid(lat2, lon2))
	print ("Expected: 622128, 185038 (according to streetmap.co.uk)")

	print("All done!")

