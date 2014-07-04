import OSTN02, OSGB, transform

if __name__=="__main__":
	fi = open("../pos.csv")
	for li in fi:
		vals = li.strip().split(",")
		#print vals
		name = vals[0]
		pos = vals[1].split(" ")[0]
		print name, pos
		print transform.OSGB36GridRefToETRS89(pos)
		
