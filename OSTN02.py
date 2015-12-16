import os
import bz2, six

#OSTN02 for Python
#=================

#This is a port of the perl module Geo::Coordinates::OSTN02 by Toby Thurston (c) 2008
#Toby kindly allowed his code to be used for any purpose.
#The python port is (c) 2010-2011 Tim Sheerman-Chase
#The OSTN02 transform is Crown Copyright (C) 2002
#See COPYING for redistribution terms

MAX_EASTING  =  700000
MAX_NORTHING = 1250000

MIN_X_SHIFT =  86.275
MIN_Y_SHIFT = -81.603
MIN_Z_SHIFT =  43.982

def ostn():
    filename = os.path.join(os.path.dirname(__file__), 'ostn02data.txt.bz2')
    fi = bz2.BZ2File(filename)
    lines = fi.readlines()
    out = {}
    for line in lines:
        line = line.rstrip(b'\r\n').decode("ascii")
        #six.print_(line)
        ne = line[:6]
        offset = (int(line[6:10],16),int(line[10:14],16),int(line[14:18],16))
        out[ne] = offset
    return out

ostn_data = ostn() # load all the data from below
ostn_shift_for= {}


def ETRS89_to_OSGB36(x,y,z=0.0):

    if ( 0 <= x and x <= MAX_EASTING and 0 <= y and y <= MAX_NORTHING ):
        (dx, dy, dz) = _find_OSTN02_shifts_at(x,y)
        (x, y, z) = _round_to_nearest_mm(x+dx, y+dy, z-dz) # note $z sign differs
    
    else:
        raise Exception('OSTN02 is not defined at '+str(x)+', '+str(y)+')')


    return (x, y, z)

def OSGB36_to_ETRS89 (x0, y0, z0 = 0.0):
    epsilon = 0.00001
    (dx, dy, dz) = _find_OSTN02_shifts_at(x0,y0)
    (x,  y,  z ) = (x0-dx, y0-dy, z0+dz)
    (last_dx, last_dy) = (dx, dy)
    #APPROX:
    while 1:
        (dx, dy, dz) = _find_OSTN02_shifts_at(x,y)
        (x, y) = (x0-dx, y0-dy)
        if abs(dx-last_dx)<epsilon and abs(dy-last_dy)<epsilon: break #last APPROX 
        (last_dx, last_dy) = (dx, dy)

    (x, y, z) = _round_to_nearest_mm(x0-dx, y0-dy, z0+dz)

    return (x, y, z)


def _round_to_nearest_mm(x,  y,  z):

    x = int(round(x*1000.))/1000.
    y = int(round(y*1000.))/1000.
    z = int(round(z*1000.))/1000.
    return (x, y, z)


def _find_OSTN02_shifts_at(x,y):

    e_index = int(x/1000.)
    n_index = int(y/1000.)

    s0_ref = _get_ostn_ref(e_index+0, n_index+0)
    #print s0_ref
    s1_ref = _get_ostn_ref(e_index+1, n_index+0)
    s2_ref = _get_ostn_ref(e_index+0, n_index+1)
    s3_ref = _get_ostn_ref(e_index+1, n_index+1)

    if s0_ref is None or s1_ref is None or s2_ref is None or s3_ref is None:
        raise Exception("[OSTN02 not defined at ("+str(x)+","+str(y)+")]")

    x0 = e_index * 1000
    y0 = n_index * 1000

    dx = x - x0 # offset within square
    dy = y - y0

    t = dx/1000
    u = dy/1000

    f0 = (1-t)*(1-u)
    f1 =    t *(1-u)
    f2 = (1-t)*   u
    f3 =    t *   u

    se = f0*s0_ref[0] + f1*s1_ref[0] + f2*s2_ref[0] + f3*s3_ref[0]
    sn = f0*s0_ref[1] + f1*s1_ref[1] + f2*s2_ref[1] + f3*s3_ref[1]
    sg = f0*s0_ref[2] + f1*s1_ref[2] + f2*s2_ref[2] + f3*s3_ref[2]

    #if se*sn*sg==0.:
    #    print("[OSTN02 defined as zeros at ($x, $y), coordinates unchanged]")

    return (se, sn, sg)

def _get_ostn_ref(x,y):

    key = "%03x%03x" % (y, x)
    if key in ostn_shift_for:
        return ostn_shift_for[key]

    if key in ostn_data:
        data = ostn_data[key]
        data2 = (data[0]/1000.0 + MIN_X_SHIFT,data[1]/1000.0 +MIN_Y_SHIFT,data[2]/1000.0 + MIN_Z_SHIFT)
        ostn_shift_for[key] = data2
        return data2



