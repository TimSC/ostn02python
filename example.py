
#OSTN02 for Python
#=================

#This is a port of the perl module Geo::Coordinates::OSTN02 by Toby Thurston 
#Toby kindly allowed his code to be used for any purpose.
#The python port is (c) 2010 Tim Sheerman-Chase
#The OSTN02 transform is Crown Copyright (C) 2002
#See COPYING for redistribution terms

from OSGB import *
from OSTN02 import *

print "OSGB36_to_ETRS89"
print "================"
print "Take OS map reference: TR143599"
xin = 614350
yin = 159950
print "OS X (Eastings) "+str(xin)
print "OS Y (Northings) "+str(yin)

(x,y,h) = OSGB36_to_ETRS89 (xin, yin)

print "Using the OSGB36_to_ETRS89 conversion gives us the grid position:"
print str(x) + "," + str(y) + "," + str(h)

(gla, glo) = grid_to_ll(x, y)

print "The grid position converts to ETRS89 lat,lon (using grid_to_ll) of:"
print str((gla, glo))

print "Actual answer: 51.29831006, 1.07337394, 44.621"

print "ETRS89 is within a metre of WGS84 (as used by GPS receivers), at time of writing (2011)."

print "\nETRS89_to_OSGB36"
print "================="

gla = 51.29831006
glo = 1.07337394
h = 44.621


print "To ETRS89 grid (using ll_to_grid):"
(x2,y2) = ll_to_grid(gla, glo)
print str((x2,y2))

print "To OS Eastings/Northings (using ETRS89_to_OSGB36):"
print ETRS89_to_OSGB36(x2,y2,h)

print "Actual Answer: 614350, 159950, 0"

print "\nExceptions"
print "=========="
print "Some areas do not have OSTN02 coverage for example OS grid 622129,185038"
print "This causes an exception to be raised, but it can be caught:"
try:
	(x,y,h) = OSGB36_to_ETRS89 (622129,185038)
except Exception as e:
	print 'Exception occurred, value:', e

print "All done!"

