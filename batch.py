import transform, csv, gpxutils, kmlutils, codecs

#Python gets rather confused with utf-8 encoded csv. These are helper classes
#to manage the problem
def unicode_csv_reader(unicode_csv_data, dialect=csv.excel, **kwargs):
    # csv.py doesn't do Unicode; encode temporarily as UTF-8:
    csv_reader = csv.reader(utf_8_encoder(unicode_csv_data),
                            dialect=dialect, **kwargs)
    for row in csv_reader:
        # decode UTF-8 back to Unicode, cell by cell:
        yield [unicode(cell, 'utf-8') for cell in row]

def utf_8_encoder(unicode_csv_data):
    for line in unicode_csv_data:
        yield line.encode('utf-8')

#Main program
if __name__=="__main__":
	#Read CSV file
	fiRows = list(unicode_csv_reader(codecs.open("pos_essex2.csv", "r", "utf-8")))

	#Output file
	outKml = kmlutils.KmlWriter("out.kml")
	outGpx = gpxutils.GpxWriter("out.gpx")

	#Use first row as dictionary keys
	fi = [dict(zip(fiRows[0], tmp)) for tmp in fiRows[1:]]
	
	centrePos = transform.OSGB36GridRefToOSFB36EastNorth("TM390914")

	for li in fi:

		name = li["NAME,C,40"]
		condition = li["CONDITION,N,9,0"]
		rating = li["Rating"]
		url = li["URL"]
		typ = li["TYPE,C,60"]
		ambience = li["AMBIENCE,N,9,0"]	
		access = li["ACCESS,N,9,0"]	
		colour = li["COLOUR,C,9"]
		info = li["INFO"]

		#Isolate the map reference
		pos = li["MAP_REF,C,9"].split(" ")[0]

		#Remove unwanted characters
		pos = pos.replace(';',"")

		print name, pos

		if len(pos) == 0:
			#Skip zero length map references
			continue

		lat, lon = transform.OSGB36GridRefToETRS89(str(pos))

		#Get easting and northing on OS grid (removes letter code of grid square)
		e, n = transform.OSGB36GridRefToOSFB36EastNorth(str(pos))
		
		#Distance in km using pythagoras
		dist = round(((e - centrePos[0]) ** 2. + (n - centrePos[1]) ** 2.) ** 0.5 / 1000., 1)

		print dist

		formattedName = u"{0}".format(name)
		description = []
		description.append(u"Type: {0}<br/>".format(typ))
		description.append(u"Info: {0}<br/>".format(info))
		description.append(u"rating: {0}<br/>".format(rating))
		description.append(u"Condition: {0}<br/>".format(condition))
		description.append(u"Ambience: {0}<br/>".format(ambience))
		description.append(u"Access: {0}<br/>".format(access))
		description.append(u"\n")

		descriptionStr = "".join(description)
		print descriptionStr

		outKml.Waypoint(lat, lon, formattedName, description=descriptionStr)
		outGpx.Waypoint(lat, lon, formattedName, description=descriptionStr)

	#Probably not necessary but just in case
	del outKml
	del outGpx

