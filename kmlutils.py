import codecs
from xml.sax.saxutils import escape

class KmlWriter(object):

	def __init__(self, filename):
		self.outFile = codecs.open(filename, "w", "utf-8")

		self.outFile.write('<?xml version="1.0" encoding="UTF-8"?>\n')

		self.outFile.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
		self.outFile.write('<Document>\n')

	def __del__(self):
		self.outFile.flush()

		self.outFile.write('</Document>\n')
		self.outFile.write('</kml>\n')

	def Waypoint(self, lat, lon, name, description = None, alt = 0.):

		self.outFile.write('<Placemark>\n')

		self.outFile.write('<name>{0}</name>\n'.format(escape(name)))
		if description is not None:
			self.outFile.write('<description>\n')
			self.outFile.write('<![CDATA[\n')
			self.outFile.write(description)
			self.outFile.write(']]>\n')
			self.outFile.write('</description>\n')
		self.outFile.write('<Point>\n')
		self.outFile.write('<coordinates>{0},{1},{2}</coordinates>\n'.format(float(lon), float(lat), float(alt)))
		self.outFile.write('</Point>\n')
		self.outFile.write('</Placemark>\n')

