
#OSTN02 for Python
#=================

#This is a port of the perl module Geo::Coordinates::OSTN02 by Toby Thurston 
#Toby kindly allowed his code to be used for any purpose.
#The python port is (c) 2010 Tim Sheerman-Chase
#The OSTN02 transform is Crown Copyright (C) 2002
#See COPYING for redistribution terms

from OSGB import *
from OSTN02 import *
import transform, six

six.print_("OSGB36_to_ETRS89")
six.print_("================")
six.print_("Take OS map reference: TR143599")

xin, yin = parse_grid("TR", 14300, 59900)

six.print_("OS X (Eastings) "+str(xin))
six.print_("OS Y (Northings) "+str(yin))

(x,y,h) = OSGB36_to_ETRS89 (xin, yin)

six.print_("Using the OSGB36_to_ETRS89 conversion gives us the grid position:")
six.print_(str(x) + "," + str(y) + "," + str(h))

(gla, glo) = grid_to_ll(x, y)

six.print_("The grid position converts to ETRS89 lat,lon (using grid_to_ll) of:")
six.print_(str((gla, glo)))

six.print_("Actual answer: 51.297880, 1.072628")

six.print_("ETRS89 is within a metre of WGS84 (as used by GPS receivers), at time of writing (2011).")

six.print_("\nGridRef to ETRS89")
six.print_("==================")

ref = "TR143599"
six.print_("Starting point: "+ref)
six.print_(transform.OSGB36GridRefToETRS89(ref))
six.print_("Expected:", 51.297880, 1.072628)

six.print_("\nETRS89_to_OSGB36")
six.print_("=================")

gla = 51.297880
glo = 1.072628
h = 44.621

six.print_("To ETRS89 grid (using ll_to_grid):")
(x2,y2) = ll_to_grid(gla, glo)
six.print_(str((x2,y2)))

six.print_("To OS Eastings/Northings (using ETRS89_to_OSGB36):")
six.print_(ETRS89_to_OSGB36(x2,y2,h))

six.print_("Actual Answer: 614300, 159900, 0")

six.print_("\nExceptions")
six.print_("==========")
six.print_("Some areas do not have OSTN02 coverage for example OS grid 622129,185038")
six.print_("This causes an exception to be raised, but it can be caught:")
try:
	(x,y,h) = OSGB36_to_ETRS89 (622129,185038)
except Exception as e:
	six.print_('Exception occurred, value:', e)

six.print_("All done!")

