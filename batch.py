import transform, csv, gpxutils

if __name__=="__main__":
	fi = open("../pos.csv")
	out = gpxutils.GpxWriter("out.gpx")

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
		
		out.Waypoint(lat, lon, name)

	del out

