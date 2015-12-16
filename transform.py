import OSTN02, OSGB, six

def OSGB36GridRefToETRS89(mapRef):

	if len(mapRef) < 4:
		raise ValueError("Map ref too short")
	if len(mapRef) % 2 == 1:
		six.print_(ValueError("Unexpected input length"))

	coordLen = int((len(mapRef) - 2) / 2)
	code = mapRef[:2]
	east = int(mapRef[2:2+coordLen])*pow(10,5-coordLen)
	nrth = int(mapRef[2+coordLen:])*pow(10,5-coordLen)

	x1, y1 = OSGB.parse_grid(code, east, nrth)

	(x,y,h) = OSTN02.OSGB36_to_ETRS89 (x1, y1)
	(gla, glo) = OSGB.grid_to_ll(x, y)
	return (gla, glo)

