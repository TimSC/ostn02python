import transform, csv, gpxutils

if __name__=="__main__":
	fi = open("../pos.csv")
	out = gpxutils.GpxWriter("out.gpx")

	centrePos = transform.OSGB36GridRefToOSFB36EastNorth("TM390914")

	for li in csv.reader(fi):
		name = li[0]

		#Isolate the map reference
		pos = li[1].split(" ")[0]

		#Remove unwanted characters
		pos = pos.replace(';',"")

		print name, pos

		if len(pos) == 0:
			#Skip zero length map references
			continue

		lat, lon = transform.OSGB36GridRefToETRS89(pos)

		#Get easting and northing on OS grid (removes letter code of grid square)
		e, n = transform.OSGB36GridRefToOSFB36EastNorth(pos)
		
		#Distance in km using pythagoras
		dist = round(((e - centrePos[0]) ** 2. + (n - centrePos[1]) ** 2.) ** 0.5 / 1000., 1)

		print dist

		formattedName = "{0} {1}".format(name,dist)

		out.Waypoint(lat, lon, formattedName)

	del out

