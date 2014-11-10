
#OSTN02 for Python
#=================

#This is a port of the perl module Geo::Coordinates::OSTN02 by Toby Thurston (c) 2008
#Toby kindly allowed his code to be used for any purpose.
#The python port is (c) 2010-2011 Tim Sheerman-Chase
#The OSTN02 transform is Crown Copyright (C) 2002
#See COPYING for redistribution terms

import math

RAD = math.pi / 180
DAR = 180 / math.pi

WGS84_MAJOR_AXIS = 6378137.000
WGS84_FLATTENING = 1.0 / 298.257223563

# set defaults for Britain
ellipsoid_shapes = {
	'WGS84': [ 6378137.0000, 6356752.31425 ],
	'ETRS89' : [ 6378137.0000, 6356752.31425 ],
	'ETRN89' : [ 6378137.0000, 6356752.31425 ],
	'GRS80'  : [ 6378137.0000, 6356752.31425 ],
	'OSGB36' : [ 6377563.396,  6356256.910  ],
	'OSGM02' : [ 6377563.396,  6356256.910  ] }
# yes lots of synonyms

# constants for OSGB mercator projection
LAM0 = RAD * -2.0  # lon of grid origin
PHI0 = RAD * 49.0  # lat of grid origin
E0 =  400000   # Easting for origin
N0 = -100000   # Northing for origin
F0 = 0.9996012717 # Convergence factor

# useful patterns
#our $Real_Pattern	= qr /^[-+][\.\d]+$/;
#our $GSq_Pattern	 = qr /[GHJMNORST][A-Z]/i;
#our $LR_Pattern	  = qr /^(\d{1,3})\D+(\d{3})\D?(\d{3})$/;
#our $GR_Pattern	  = qr /^($GSq_Pattern)\s?(\d{3})\D?(\d{3})$/;
#our $Long_GR_Pattern = qr /^($GSq_Pattern)\s?(\d{5})\D?(\d{5})$/;
#our $ISO_LL_Pattern = qr{\A
#						([-+])(	\d{2,6})(?:\.(\d+))?
#						([-+])([01]\d{2,6})(?:\.(\d+))?
#						([-+][\.\d]+)?
#						\/
#						\Z}xo;

def ll_to_grid(lat,lon,alt=0.0,shape = 'WGS84'):

	#my $shape = defined $ellipsoid_shapes{$_[-1]} ? pop : 'OSGB36'; # last argument (or omitted)
	#if ($lat =~ $ISO_LL_Pattern ) {
	#	($lat, $lon, $alt) = parse_ISO_ll($lat);
	#}

	(a,b) = ellipsoid_shapes[shape]
	
	e2 = (a**2.-b**2.)/a**2.
	n = (a-b)/(a+b)
	
	phi = RAD * lat
	lam = RAD * lon

	sp2  = math.sin(phi)**2.
	nu   = a * F0 * (1. - e2 * sp2 ) ** -0.5
	rho  = a * F0 * (1. - e2) * (1. - e2 * sp2 ) ** -1.5
	eta2 = nu/rho - 1.
	
	M = _compute_M(phi, b, n)

	cp = math.cos(phi)
	sp = math.sin(phi)
	tp = math.tan(phi)
	tp2 = tp*tp
	tp4 = tp2*tp2

	I	= M + N0
	II   = nu/2.  * sp * cp
	III  = nu/24. * sp * cp**3. * (5.-tp2+9.*eta2)
	IIIA = nu/720.* sp * cp**5. *(61.-58.*tp2+tp4)

	IV   = nu*cp
	V	= nu/6.   * cp**3. * (nu/rho-tp2)
	VI   = nu/120. * cp**5. * (5.-18.*tp2+tp4+14.*eta2-58.*tp2*eta2)

	l = lam - LAM0
	north = I  + II*l**2. + III*l**4. + IIIA*l**6.;
	east  = E0 + IV*l	+   V*l**3. +   VI*l**5.;

	# round to 3dp (mm)
	#($east, $north) = map { sprintf "%.3f", $_ } ($east, $north);

	return (east,north)

def grid_to_ll(E,N,shape='WGS84'):

	#if ( $E =~ $GR_Pattern || $E =~ $Long_GR_Pattern || $E =~ $LR_Pattern ) {
	#	($E, $N) = parse_grid($E);
	#}

	(a,b) = ellipsoid_shapes[shape]

	e2 = (a**2.-b**2.)/a**2.
	n = (a-b)/(a+b)

	dN = N - N0

	phi = PHI0 + dN/(a * F0)

	M = _compute_M(phi, b, n);
	while (dN-M >= 0.001):
	   phi = phi + (dN-M)/(a * F0)
	   M = _compute_M(phi, b, n)

	sp2  = math.sin(phi)**2.;
	nu   = a * F0 *			 (1. - e2 * sp2 ) ** -0.5
	rho  = a * F0 * (1. - e2) * (1. - e2 * sp2 ) ** -1.5
	eta2 = nu/rho - 1.

	tp = math.tan(phi)
	tp2 = tp*tp
	tp4 = tp2*tp2

	VII  = tp /   (2.*rho*nu)
	VIII = tp /  (24.*rho*nu**3.) *  (5. +  3.*tp2 + eta2 - 9.*tp2*eta2)
	IX   = tp / (720.*rho*nu**5.) * (61. + 90.*tp2 + 45.*tp4)

	sp = 1.0 / math.cos(phi) 
	tp6 = tp4*tp2

	X	= sp/nu
	XI   = sp/(   6.*nu**3.)*(nu/rho + 2.*tp2)
	XII  = sp/( 120.*nu**5.)*(	  5. + 28.*tp2 +   24.*tp4)
	XIIA = sp/(5040.*nu**7.)*(	61. + 662.*tp2 + 1320.*tp4 + 720.*tp6)

	e = E - E0

	phi = phi		- VII*e**2. + VIII*e**4. -   IX*e**6.
	lam = LAM0 + X*e -  XI*e**3. +  XII*e**5. - XIIA*e**7.

	phi = phi * DAR
	lam = lam * DAR

	return (phi, lam)
	#return format_ll_ISO($phi,$lam);


def _compute_M(phi, b, n):
	p_plus  = phi + PHI0
	p_minus = phi - PHI0
	return b * F0 * (
		   (1 + n * (1 + 5./4*n*(1 + n)))*p_minus
		 - 3*n*(1+n*(1+7./8*n))  * math.sin(p_minus) * math.cos(p_plus)
		 + (15./8*n * (n*(1+n))) * math.sin(2*p_minus) * math.cos(2*p_plus)
		 - 35./24*n**3			 * math.sin(3*p_minus) * math.cos(3*p_plus)
		   )

'''

our @Grid = ( [ qw( V W X Y Z ) ],
			  [ qw( Q R S T U ) ],
			  [ qw( L M N O P ) ],
			  [ qw( F G H J K ) ],
			  [ qw( A B C D E ) ] );
'''
Big_off = {#East then north
				 'G' : ( -1, 2 ),
				 'H' : ( 0, 2 ),
				 'J' : ( 1, 2 ),
				 'M' : ( -1, 1 ),
				 'N' : ( 0, 1 ),
				 'O' : ( 1, 1 ),
				 'R' : ( -1, 0 ),
				 'S' : ( 0, 0 ),
				 'T' : ( 1, 0 ),
		   }

Small_off = {#East then north
				 'A' : ( 0, 4 ),
				 'B' : ( 1, 4 ),
				 'C' : ( 2, 4 ),
				 'D' : ( 3, 4 ),
				 'E' : ( 4, 4 ),

				 'F' : ( 0, 3 ),
				 'G' : ( 1, 3 ),
				 'H' : ( 2, 3 ),
				 'J' : ( 3, 3 ),
				 'K' : ( 4, 3 ),

				 'L' : ( 0, 2 ),
				 'M' : ( 1, 2 ),
				 'N' : ( 2, 2 ),
				 'O' : ( 3, 2 ),
				 'P' : ( 4, 2 ),

				 'Q' : ( 0, 1 ),
				 'R' : ( 1, 1 ),
				 'S' : ( 2, 1 ),
				 'T' : ( 3, 1 ),
				 'U' : ( 4, 1 ),

				 'V' : ( 0, 0 ),
				 'W' : ( 1, 0 ),
				 'X' : ( 2, 0 ),
				 'Y' : ( 3, 0 ),
				 'Z' : ( 4, 0 ),
		   }

BIG_SQUARE = 500000
SQUARE	 = 100000
'''
# Landranger sheet data
# These are the full GRs (as metres from Newlyn) of the SW corner of each sheet.
our %LR = (
1   => [ 429000 ,1179000 ] ,
2   => [ 433000 ,1156000 ] ,
3   => [ 414000 ,1147000 ] ,
4   => [ 420000 ,1107000 ] ,
5   => [ 340000 ,1020000 ] ,
6   => [ 321000 , 996000 ] ,
7   => [ 315000 , 970000 ] ,
8   => [ 117000 , 926000 ] ,
9   => [ 212000 , 940000 ] ,
10  => [ 252000 , 940000 ] ,
11  => [ 292000 , 929000 ] ,
12  => [ 300000 , 939000 ] ,
13  => [  95000 , 903000 ] ,
14  => [ 105000 , 886000 ] ,
15  => [ 196000 , 900000 ] ,
16  => [ 236000 , 900000 ] ,
17  => [ 276000 , 900000 ] ,
18  => [  69000 , 863000 ] ,
19  => [ 174000 , 860000 ] ,
20  => [ 214000 , 860000 ] ,
21  => [ 254000 , 860000 ] ,
22  => [  57000 , 823000 ] ,
23  => [ 113000 , 836000 ] ,
24  => [ 150000 , 830000 ] ,
25  => [ 190000 , 820000 ] ,
26  => [ 230000 , 820000 ] ,
27  => [ 270000 , 830000 ] ,
28  => [ 310000 , 833000 ] ,
29  => [ 345000 , 830000 ] ,
30  => [ 377000 , 830000 ] ,
31  => [  50000 , 783000 ] ,
32  => [ 130000 , 800000 ] ,
33  => [ 170000 , 790000 ] ,
34  => [ 210000 , 780000 ] ,
35  => [ 250000 , 790000 ] ,
36  => [ 285000 , 793000 ] ,
37  => [ 325000 , 793000 ] ,
38  => [ 365000 , 790000 ] ,
39  => [ 120000 , 770000 ] ,
40  => [ 160000 , 760000 ] ,
41  => [ 200000 , 750000 ] ,
42  => [ 240000 , 750000 ] ,
43  => [ 280000 , 760000 ] ,
44  => [ 320000 , 760000 ] ,
45  => [ 360000 , 760000 ] ,
46  => [  92000 , 733000 ] ,
47  => [ 120000 , 732000 ] ,
48  => [ 120000 , 710000 ] ,
49  => [ 160000 , 720000 ] ,
50  => [ 200000 , 710000 ] ,
51  => [ 240000 , 720000 ] ,
52  => [ 270000 , 720000 ] ,
53  => [ 294000 , 720000 ] ,
54  => [ 334000 , 720000 ] ,
55  => [ 164000 , 680000 ] ,
56  => [ 204000 , 682000 ] ,
57  => [ 244000 , 682000 ] ,
58  => [ 284000 , 690000 ] ,
59  => [ 324000 , 690000 ] ,
60  => [ 110000 , 640000 ] ,
61  => [ 131000 , 662000 ] ,
62  => [ 160000 , 640000 ] ,
63  => [ 200000 , 642000 ] ,
64  => [ 240000 , 645000 ] ,
65  => [ 280000 , 650000 ] ,
66  => [ 316000 , 650000 ] ,
67  => [ 356000 , 650000 ] ,
68  => [ 157000 , 600000 ] ,
69  => [ 175000 , 613000 ] ,
70  => [ 215000 , 605000 ] ,
71  => [ 255000 , 605000 ] ,
72  => [ 280000 , 620000 ] ,
73  => [ 320000 , 620000 ] ,
74  => [ 357000 , 620000 ] ,
75  => [ 390000 , 620000 ] ,
76  => [ 195000 , 570000 ] ,
77  => [ 235000 , 570000 ] ,
78  => [ 275000 , 580000 ] ,
79  => [ 315000 , 580000 ] ,
80  => [ 355000 , 580000 ] ,
81  => [ 395000 , 580000 ] ,
82  => [ 195000 , 530000 ] ,
83  => [ 235000 , 530000 ] ,
84  => [ 265000 , 540000 ] ,
85  => [ 305000 , 540000 ] ,
86  => [ 345000 , 540000 ] ,
87  => [ 367000 , 540000 ] ,
88  => [ 407000 , 540000 ] ,
89  => [ 290000 , 500000 ] ,
90  => [ 317000 , 500000 ] ,
91  => [ 357000 , 500000 ] ,
92  => [ 380000 , 500000 ] ,
93  => [ 420000 , 500000 ] ,
94  => [ 460000 , 485000 ] ,
95  => [ 213000 , 465000 ] ,
96  => [ 303000 , 460000 ] ,
97  => [ 326000 , 460000 ] ,
98  => [ 366000 , 460000 ] ,
99  => [ 406000 , 460000 ] ,
100 => [ 446000 , 460000 ] ,
101 => [ 486000 , 460000 ] ,
102 => [ 326000 , 420000 ] ,
103 => [ 360000 , 420000 ] ,
104 => [ 400000 , 420000 ] ,
105 => [ 440000 , 420000 ] ,
106 => [ 463000 , 420000 ] ,
107 => [ 500000 , 420000 ] ,
108 => [ 320000 , 380000 ] ,
109 => [ 360000 , 380000 ] ,
110 => [ 400000 , 380000 ] ,
111 => [ 430000 , 380000 ] ,
112 => [ 470000 , 385000 ] ,
113 => [ 510000 , 386000 ] ,
114 => [ 220000 , 360000 ] ,
115 => [ 240000 , 345000 ] ,
116 => [ 280000 , 345000 ] ,
117 => [ 320000 , 340000 ] ,
118 => [ 360000 , 340000 ] ,
119 => [ 400000 , 340000 ] ,
120 => [ 440000 , 350000 ] ,
121 => [ 478000 , 350000 ] ,
122 => [ 518000 , 350000 ] ,
123 => [ 210000 , 320000 ] ,
124 => [ 250000 , 305000 ] ,
125 => [ 280000 , 305000 ] ,
126 => [ 320000 , 300000 ] ,
127 => [ 360000 , 300000 ] ,
128 => [ 400000 , 308000 ] ,
129 => [ 440000 , 310000 ] ,
130 => [ 480000 , 310000 ] ,
131 => [ 520000 , 310000 ] ,
132 => [ 560000 , 310000 ] ,
133 => [ 600000 , 310000 ] ,
134 => [ 617000 , 290000 ] ,
135 => [ 250000 , 265000 ] ,
136 => [ 280000 , 265000 ] ,
137 => [ 320000 , 260000 ] ,
138 => [ 345000 , 260000 ] ,
139 => [ 385000 , 268000 ] ,
140 => [ 425000 , 270000 ] ,
141 => [ 465000 , 270000 ] ,
142 => [ 504000 , 274000 ] ,
143 => [ 537000 , 274000 ] ,
144 => [ 577000 , 270000 ] ,
145 => [ 200000 , 220000 ] ,
146 => [ 240000 , 225000 ] ,
147 => [ 270000 , 240000 ] ,
148 => [ 310000 , 240000 ] ,
149 => [ 333000 , 228000 ] ,
150 => [ 373000 , 228000 ] ,
151 => [ 413000 , 230000 ] ,
152 => [ 453000 , 230000 ] ,
153 => [ 493000 , 234000 ] ,
154 => [ 533000 , 234000 ] ,
155 => [ 573000 , 234000 ] ,
156 => [ 613000 , 250000 ] ,
157 => [ 165000 , 201000 ] ,
158 => [ 189000 , 190000 ] ,
159 => [ 229000 , 185000 ] ,
160 => [ 269000 , 205000 ] ,
161 => [ 309000 , 205000 ] ,
162 => [ 349000 , 188000 ] ,
163 => [ 389000 , 190000 ] ,
164 => [ 429000 , 190000 ] ,
165 => [ 460000 , 195000 ] ,
166 => [ 500000 , 194000 ] ,
167 => [ 540000 , 194000 ] ,
168 => [ 580000 , 194000 ] ,
169 => [ 607000 , 210000 ] ,
170 => [ 269000 , 165000 ] ,
171 => [ 309000 , 165000 ] ,
172 => [ 340000 , 155000 ] ,
173 => [ 380000 , 155000 ] ,
174 => [ 420000 , 155000 ] ,
175 => [ 460000 , 155000 ] ,
176 => [ 495000 , 160000 ] ,
177 => [ 530000 , 160000 ] ,
178 => [ 565000 , 155000 ] ,
179 => [ 603000 , 133000 ] ,
180 => [ 240000 , 112000 ] ,
181 => [ 280000 , 112000 ] ,
182 => [ 320000 , 130000 ] ,
183 => [ 349000 , 115000 ] ,
184 => [ 389000 , 115000 ] ,
185 => [ 426000 , 116000 ] ,
186 => [ 465000 , 125000 ] ,
187 => [ 505000 , 125000 ] ,
188 => [ 545000 , 125000 ] ,
189 => [ 585000 , 115000 ] ,
190 => [ 207000 ,  87000 ] ,
191 => [ 247000 ,  72000 ] ,
192 => [ 287000 ,  72000 ] ,
193 => [ 310000 ,  90000 ] ,
194 => [ 349000 ,  75000 ] ,
195 => [ 389000 ,  75000 ] ,
196 => [ 429000 ,  76000 ] ,
197 => [ 469000 ,  90000 ] ,
198 => [ 509000 ,  97000 ] ,
199 => [ 549000 ,  94000 ] ,
200 => [ 175000 ,  50000 ] ,
201 => [ 215000 ,  47000 ] ,
202 => [ 255000 ,  32000 ] ,
203 => [ 132000 ,  11000 ] ,
204 => [ 172000 ,  14000 ] ,
);

sub format_grid_trad {
	use integer;
	my ($sq, $e, $n) = format_grid_GPS(@_);
	($e,$n) = ($e/100,$n/100);
	return ($sq, $e, $n) if wantarray;
	return sprintf "%s %03d %03d", $sq, $e, $n;
}

sub format_grid_GPS {
	my $e = shift;
	my $n = shift;

	croak "Easting must not be negative" if $e<0;
	croak "Northing must not be negative" if $n<0;

	# round to nearest metre
	($e,$n) = map { $_+0.5 } ($e, $n);
	my $sq;

	{
		use integer;
		$sq = sprintf "%s%s", _letter( 2 + $e/BIG_SQUARE		 , 1+$n/BIG_SQUARE		),
							  _letter(($e % BIG_SQUARE ) / SQUARE, ( $n % BIG_SQUARE )/SQUARE );

		($e,$n) = map { $_ % SQUARE } ($e, $n);
	}

	return ($sq, $e, $n) if wantarray;
	return sprintf "%s %05d %05d", $sq, $e, $n;
}

sub format_grid_landranger {
	use integer;
	my ($e,$n) = @_;
	my @sheets = ();
	for my $sheet (1..204) {
		my $de = $e-$LR{$sheet}->[0];
		my $dn = $n-$LR{$sheet}->[1];
		push @sheets, $sheet if $de>=0 && $de < 40000
							 && $dn>=0 && $dn < 40000;
	}
	my $sq;
	($sq, $e, $n) = format_grid_trad($e,$n);

	return ($sq, $e, $n, @sheets) if wantarray;

	return sprintf("%s %03d %03d is not on any OS Sheet", $sq, $e, $n) unless @sheets;
	return sprintf("%s %03d %03d on OS Sheet %d"		, $sq, $e, $n, $sheets[0]) if 1==@sheets;
	return sprintf("%s %03d %03d on OS Sheets %d and %d", $sq, $e, $n, @sheets)	if 2==@sheets;
	return sprintf("%s %03d %03d on OS Sheets %s", $sq, $e, $n, join(', ', @sheets[0..($#sheets-1)], "and $sheets[-1]"));

}

sub _letter {
	my $x = shift;
	my $y = shift;
	die "Argument out of range in _letter\n"
		unless defined $x && $x=~/^\d+$/ && $x>=0 && $x<5
			&& defined $y && $y=~/^\d+$/ && $y>=0 && $y<5;

	return $Grid[$y][$x];
}


sub parse_grid {
	my $s = shift;
	return parse_trad_grid($s) if $s =~ $GR_Pattern;
	return parse_GPS_grid($s)  if $s =~ $Long_GR_Pattern;
	return parse_landranger_grid($1, $2, $3) if $s =~ $LR_Pattern;
	return parse_landranger_grid($s) if $s =~ /^\d{1,3}$/ && $s < 205;
	croak "$s <-- this does not match my grid ref patterns";
	return
}

sub parse_trad_grid {
	my ($letters, $e, $n);
	if	( @_ == 1 && $_[0] =~ $GR_Pattern ) {
		($letters, $e, $n) = ($1,$2,$3)
	}
	elsif ( @_ == 2 && $_[0] =~ $GSq_Pattern && $_[1] =~ /^(\d{3})(\d{3})$/ ) {
		$letters = $_[0]; ($e, $n) = ($1,$2)
	}
	elsif ( @_ == 3 && $_[0] =~ $GSq_Pattern && $_[1] =~ /^\d{1,3}$/
											 && $_[2] =~ /^\d{1,3}$/ ) {
		($letters, $e, $n) = @_
	}
	else { croak "Cannot parse @_ as a traditional grid reference"; }

	return _parse_grid($letters, $e*100, $n*100)
}

sub parse_GPS_grid {
	my ($letters, $e, $n);
	if	( @_ == 1 && $_[0] =~ $Long_GR_Pattern ) {
		($letters, $e, $n) = ($1,$2,$3)
	}
	elsif ( @_ == 2 && $_[0] =~ $GSq_Pattern && $_[1] =~ /^(\d{5})(\d{5})$/ ) {
		$letters = $_[0]; ($e, $n) = ($1,$2)
	}
	elsif ( @_ == 3 && $_[0] =~ $GSq_Pattern && $_[1] =~ /^\d{5}$/ && $_[2] =~ /^\d{5}$/ ) {
		($letters, $e, $n) = @_
	}
	else { croak "Cannot parse @_ as a GPS-style grid reference"; }

	return _parse_grid($letters, $e, $n)
}
'''

def parse_grid (letters,e=0.0,n=0.0):
	"""Convert from OS grid references with letter prefix to the national grid.
	e and n are in metres.
	e.g. OSGB.parse_grid("ST", 00000, 00000) converts to (300000, 100000)
	NY462754 is OSGB.parse_grid("NY", 46200, 75400) and converts to (346200, 575400)
	"""

	letters = str.upper(letters)

	c = letters[0:1].upper()
	e += Big_off[c][0]*BIG_SQUARE
	n += Big_off[c][1]*BIG_SQUARE
	
	d = letters[1:2].upper()
	e += Small_off[d][0]*SQUARE
	n += Small_off[d][1]*SQUARE

	return (e, n)

def grid_to_small_code(e, n):
	er = int(e % BIG_SQUARE) / SQUARE
	nr = int(n % BIG_SQUARE) / SQUARE
	found = None
	for cd in Small_off:
		if er == Small_off[cd][0] and nr == Small_off[cd][1]:
			found = cd
	ebig = e - er *  SQUARE
	nbig = n - nr *  SQUARE
	return found, ebig, nbig

def grid_to_big_code(e, n):
	found = None
	for cd in Big_off:
		if e == Big_off[cd][0] and n == Big_off[cd][1]:
			found = cd
	return found

def os_streetview_tile_to_grid(tile_name):
	#Convert OS Street View tile name (e.g. SO02NW) to grid (e.g. 
	e,n = parse_grid(tile_name)
	e += 10000 * int(tile_name[2:3])
	n += 10000 * int(tile_name[3:4])
	
	if tile_name[4:6].upper() == "NW": 
		n += 5000
	if tile_name[4:6].upper() == "NE": 
		n += 5000
		e += 5000
	if tile_name[4:6].upper() == "SE": 
		e += 5000
	return e, n

def grid_to_os_streetview_tile(grid):

	e = grid[0]
	n = grid[1]

	#Deal with small offsets
	e, eSmall = divmod(e, 5000)
	n, nSmall = divmod(n, 5000)
	e = int(e * 5000)
	n = int(n * 5000)

	#Find which corner of the tile we are in
	eRem = e % 10000
	nRem = n % 10000
	if eRem >= 5000.: 
		etile = True
		e -= 5000
	else: etile = False
	if nRem >= 5000.: 
		ntile = True
		n -= 5000
	else: ntile = False
	corner = None
	if etile and ntile: corner = "NE"
	if etile and not ntile: corner = "SE"
	if not etile and ntile: corner = "NW"
	if not etile and not ntile: corner = "SW"

	#Mid level offset
	eMidOffset = int(e % 100000) / 10000
	nMidOffset = int(n % 100000) / 10000
	e -= eMidOffset * 10000
	n -= nMidOffset * 10000

	smallCode, e, n = grid_to_small_code(e, n)
	bigCode = grid_to_big_code(e / BIG_SQUARE, n / BIG_SQUARE)

	codeOut = "{0}{1}{2}{3}{4}".format(bigCode, smallCode, eMidOffset, nMidOffset, corner)
	return codeOut, eSmall, nSmall

'''
sub _get_en {
	my $e = shift || croak "You need to supply a grid reference";
	my $n = shift;
	if ( $e =~ /^(\d{3})(\d{3})$/ && not defined $n  ) { return ($1*100, $2*100) }
	if ( $e =~ /^\d{3}$/		  && $n =~ /^\d{3}$/ ) { return ($e*100, $n*100) }
	if ( $e =~ /^\d{4}$/		  && $n =~ /^\d{4}$/ ) { return ($e*10,  $n*10 ) }
	if ( $e =~ /^\d{5}$/		  && $n =~ /^\d{5}$/ ) { return ($e*1,   $n*1  ) }
	croak "I was expecting a grid reference, not this: @_";
}


sub parse_landranger_grid {

	return unless defined wantarray;

	my $sheet = shift;

	croak "You need to supply an OS Sheet number" unless defined $sheet;
	croak "I do not know the OS Sheet number ($sheet) that you have given"
	unless defined $LR{$sheet};

	unless (@_) {
		return wantarray ? @{$LR{$sheet}} : format_grid_trad(@{$LR{$sheet}});
	}

	use integer;

	my ($e,$n) = &_get_en; # convert grid refs to metres

	my ($lle,$lln) = @{$LR{$sheet}};

	# offset from start, corrected if we are in the next 100km sq
	my $offset = $e - $lle%100_000 ; $offset += 100_000 if $offset < 0;

	if ( $offset >= 40_000 ) {
		croak sprintf "The easting you have given is %.1f km east of Sheet %s", $offset/1000-40, $sheet;
	}
	elsif ( $offset < 0 ) {
		croak sprintf "The easting you have given is %.1f km west of Sheet %s", abs($offset/1000), $sheet;
	}
	else {
		$e = $lle + $offset;
	}

	# now the same for the northing
	$offset = $n - $lln%100_000 ; $offset += 100_000 if $offset < 0;
	if ( $offset >= 40_000 ) {
		croak sprintf "The northing you have given is %.1f km north of Sheet %s", $offset/1000-40, $sheet;
	}
	elsif ( $offset < 0 ) {
		croak sprintf "The northing you have given is %.1f km south of Sheet %s", abs($offset/1000), $sheet;
	}
	else {
		$n = $lln + $offset;
	}

	return ($e, $n);

}

sub parse_ISO_ll {
	return unless defined wantarray;
	my $ISO_string = shift;

	croak "I can't parse  an ISO 6709 lat/lon string from your input ($ISO_string)" unless
	my ($lat_sign, $lat_ip, $lat_fp,
		$lon_sign, $lon_ip, $lon_fp, $alt ) = $ISO_string =~ m{$ISO_LL_Pattern};

	my $l_lat = length($lat_ip);
	my $l_lon = length($lon_ip);
	croak "Not ISO 6709 string: $ISO_string" unless $l_lat%2==0 && $l_lon%2==1; # (2,4,6) + (3,5,7)
	croak "Bad ISO 6709 string: $ISO_string" unless $l_lon-$l_lat == 1; # (2,3) (4,5) or (6,7)

	my ($lat, $lon) = (0,0);
	$lat_fp = (defined $lat_fp) ? ".$lat_fp" : '';
	$lon_fp = (defined $lon_fp) ? ".$lon_fp" : '';
	if	( $l_lat == 2 ) {
		$lat = $lat_ip.$lat_fp;
		$lon = $lon_ip.$lon_fp;
	}
	elsif ( $l_lat == 4 ) {
		$lat = substr($lat_ip,0,2) + ( substr($lat_ip,2,2) . $lat_fp ) / 60;
		$lon = substr($lon_ip,0,3) + ( substr($lon_ip,3,2) . $lon_fp ) / 60;
	}
	else {
		$lat = substr($lat_ip,0,2) + substr($lat_ip,2,2)/60 + ( substr($lat_ip,4,2) . $lat_fp )/3600;
		$lon = substr($lon_ip,0,3) + substr($lon_ip,3,2)/60 + ( substr($lon_ip,5,2) . $lon_fp )/3600;
	}

	croak "Latitude cannot exceed 90 degrees"   if $lat > 90;
	croak "Longitude cannot exceed 180 degrees" if $lon > 180;

	$lat = $lat_sign . $lat;
	$lon = $lon_sign . $lon;

	return ($lat, $lon, $alt) if wantarray;
	return format_ll_ISO($lat,$lon);
}

#	 Latitude and Longitude in Degrees:
#		 +-DD.DDDD+-DDD.DDDD/		 (eg +12.345-098.765/)
#	  Latitude and Longitude in Degrees and Minutes:
#		 +-DDMM.MMMM+-DDDMM.MMMM/	 (eg +1234.56-09854.321/)
#	  Latitude and Longitude in Degrees, Minutes and Seconds:
#		 +-DDMMSS.SSSS+-DDDMMSS.SSSS/ (eg +123456.7-0985432.1/)
#
#   where:
#
#		+-DD   = three-digit integer degrees part of latitude (through -90 ~ -00 ~ +90)
#		+-DDD  = four-digit integer degrees part of longitude (through -180 ~ -000 ~ +180)
#		MM	= two-digit integer minutes part (00 through 59)
#		SS	= two-digit integer seconds part (00 through 59)
#		.DDDD = variable-length fraction part in degrees
#		.MMMM = variable-length fraction part in minutes
#		.SSSS = variable-length fraction part in seconds
#
#		* Latitude is written in the first, and longitude is second.
#		* The sign is always necessary for each value.
#		  Latitude : North="+" South="-"
#		  Longitude: East ="+" West ="-"
#		* The integer part is a fixed length respectively.
#		  And padding character is "0".
#		  (Note: Therefor, it is shown explicitly that the first is latitude and the second is
#				 longitude, from the number of figures of the integer part.)
#		* It is variable-length below the decimal point.
#		* "/"is a terminator.
#
#   Altitude can be added optionally.
#	  Latitude, Longitude (in Degrees) and Altitude:
#		 +-DD.DDDD+-DDD.DDDD+-AAA.AAA/		 (eg +12.345-098.765+15.9/)
#	  Latitude, Longitude (in Degrees and Minutes) and Altitude:
#		 +-DDMM.MMMM+-DDDMM.MMMM+-AAA.AAA/	 (eg +1234.56-09854.321+15.9/)
#	  Latitude, Longitude (in Degrees, Minutes and Seconds) and Altitude:
#		 +-DDMMSS.SSSS+-DDDMMSS.SSSS+-AAA.AAA/ (eg +123456.7-0985432.1+15.9/)
#
#   where:
#
#		+-AAA.AAA = variable-length altitude in meters [m].
#
#		* The unit of altitude is meter [m].
#		* The integer part and the fraction part of altitude are both variable-length.
#


sub format_ll_trad {
	return unless defined wantarray;
	my ($lat, $lon) = @_;
	my ($lad, $lam, $las, $is_north) = _dms($lat); my $lah = $is_north ? 'N' : 'S';
	my ($lod, $lom, $los, $is_east ) = _dms($lon); my $loh = $is_east  ? 'E' : 'W';
	return ($lah, $lad, $lam, $las, $loh, $lod, $lom, $los) if wantarray;
	return sprintf("%s%d:%02d:%02d %s%d:%02d:%02d", $lah, $lad, $lam, $las, $loh, $lod, $lom, $los);
}

sub _dms {
	my $dd = shift;
	my $is_positive = ($dd>=0);
	$dd = abs($dd);
	my $d = int($dd);	 $dd = $dd-$d;
	my $m = int($dd*60);  $dd = $dd-$m/60;
	my $s = $dd*3600;
	return $d, $m, $s, $is_positive;
}

sub format_ll_ISO {
	return unless defined wantarray;
	my ($lat, $lon, $option) = @_;

	my ($lasign, $lad, $lam, $las) = _get_sdms($lat);
	my ($losign, $lod, $lom, $los) = _get_sdms($lon);

	# return rounded to nearest minute unless we ask for more accuracy
	if (defined $option && (uc($option) eq 'SECONDS')) {
		return ($lasign, $lad, $lam, $las, $losign, $lod, $lom, $los) if wantarray;
		return sprintf("%s%02d%02d%02d%s%03d%02d%02d/", $lasign, $lad, $lam, $las, $losign, $lod, $lom, $los);
	}

	($lad, $lam) = _round_up($lad, $lam, $las);
	($lod, $lom) = _round_up($lod, $lom, $los);



	return ($lasign ,$lad, $lam, $losign, $lod, $lom) if wantarray;
	return sprintf("%s%02d%02d%s%03d%02d/", $lasign ,$lad, $lam, $losign, $lod, $lom);
}

sub _round_up {
	my ($d, $m, $s) = @_;
	return unless defined wantarray;
	return ($d, $m) if $s<30;

	$m++;
	if ($m==60) {
		$m = 0;
		$d = $d+1;
	}
	return ($d, $m);
}

sub _get_sdms {
	my $r = shift;
	return unless defined wantarray;

	my $sign = $r>=0 ? '+' : '-';
	$r = abs($r);
	my $deg = int($r);
	my $exact_minutes = 60*($r-$deg);
	my $whole_minutes = int($exact_minutes);
	my $exact_seconds = 60 * ($exact_minutes-$whole_minutes);
	my $whole_seconds = int(0.5+$exact_seconds);
	if ( $whole_seconds > 59) {
		$whole_minutes++;
		$whole_seconds=0;
		if ($whole_minutes > 59 ) {
			$deg++;
			$whole_minutes = 0;
		}
	}
	return ($sign, $deg, $whole_minutes, $whole_seconds);
}

my %parameters_for_datum = (

	"OSGB36" => [ 573.604, 0.119600236/10000, 375, -111, 431 ],
	"OSGM02" => [ 573.604, 0.119600236/10000, 375, -111, 431 ],

	);

sub shift_ll_from_WGS84 {

	my ($lat, $lon, $elevation) = (@_, 0);

	my $parameter_ref = $parameters_for_datum{'OSGM02'};
	my $target_da = -1 * $parameter_ref->[0];
	my $target_df = -1 * $parameter_ref->[1];
	my $target_dx = -1 * $parameter_ref->[2];
	my $target_dy = -1 * $parameter_ref->[3];
	my $target_dz = -1 * $parameter_ref->[4];

	my $reference_major_axis = WGS84_MAJOR_AXIS;
	my $reference_flattening = WGS84_FLATTENING;

	return _transform($lat, $lon, $elevation,
					  $reference_major_axis, $reference_flattening,
					  $target_da, $target_df,
					  $target_dx, $target_dy, $target_dz);
}

sub shift_ll_into_WGS84 {
	my ($lat, $lon, $elevation) = (@_, 0);

	my $parameter_ref = $parameters_for_datum{'OSGM02'};
	my $target_da = $parameter_ref->[0];
	my $target_df = $parameter_ref->[1];
	my $target_dx = $parameter_ref->[2];
	my $target_dy = $parameter_ref->[3];
	my $target_dz = $parameter_ref->[4];

	my $reference_major_axis = WGS84_MAJOR_AXIS - $target_da;
	my $reference_flattening = WGS84_FLATTENING - $target_df;

	return _transform($lat, $lon, $elevation,
					  $reference_major_axis, $reference_flattening,
					  $target_da, $target_df,
					  $target_dx, $target_dy, $target_dz);
}

sub _transform {
	return unless defined wantarray;

	my $lat = shift;
	my $lon = shift;
	my $elev = shift || 0; # in case $elevation was passed as undef

	my $from_a = shift;
	my $from_f = shift;

	my $da = shift;
	my $df = shift;
	my $dx = shift;
	my $dy = shift;
	my $dz = shift;

	my $sin_lat = sin( $lat * RAD );
	my $cos_lat = cos( $lat * RAD );
	my $sin_lon = sin( $lon * RAD );
	my $cos_lon = cos( $lon * RAD );

	my $b_a	  = 1 - $from_f;
	my $e_sq	 = $from_f*(2-$from_f);
	my $ecc	  = 1 - $e_sq*$sin_lat*$sin_lat;

	my $Rn	   = $from_a / sqrt($ecc);
	my $Rm	   = $from_a * (1-$e_sq) / ($ecc*sqrt($ecc));

	my $d_lat = ( - $dx*$sin_lat*$cos_lon
				  - $dy*$sin_lat*$sin_lon
				  + $dz*$cos_lat
				  + $da*($Rn*$e_sq*$sin_lat*$cos_lat)/$from_a
				  + $df*($Rm/$b_a + $Rn*$b_a)*$sin_lat*$cos_lat
				) / ($Rm + $elev);


	my $d_lon = ( - $dx*$sin_lon
				  + $dy*$cos_lon
				) / (($Rn+$elev)*$cos_lat);

	my $d_elev = + $dx*$cos_lat*$cos_lon
				 + $dy*$cos_lat*$sin_lon
				 + $dz*$sin_lat
				 - $da*$from_a/$Rn
				 + $df*$b_a*$Rn*$sin_lat*$sin_lat;

	my ($new_lat, $new_lon, $new_elev) = (
		 $lat + $d_lat * DAR,
		 $lon + $d_lon * DAR,
		 $elev + $d_elev,
	   );

	return ($new_lat, $new_lon, $new_elev) if wantarray;
	return sprintf "%s, (%s m)", format_ll_ISO($new_lat, $new_lon), $new_elev;

}

1;
__END__

=head1 NAME

Geo::Coordinates::OSGB - Convert coordinates between Lat/Lon and the British National Grid

An implementation of co-ordinate conversion for England, Wales, and Scotland
based on formulae published by the Ordnance Survey of Great Britain.

These modules will convert accurately between an OSGB national grid
reference and lat/lon coordinates based on the OSGB geoid model.  (For an
explanation of what a geoid model is and why you should care, read the
L<Theory> section below.) The OSGB geoid model fits mainland Britain very
well, but is rather different from the international WGS84 model that has
rapidly become the de facto universal standard model thanks to the
popularity of GPS devices and maps on the Internet.  So, if you are trying
to translate from an OSGB grid reference to lat/lon coordinates that can be
used in Google Earth, Wikipedia, or some other Internet based tool, you will
need to do two transformations:  first translate your grid ref into OSGB
lat/lon; then nudge the result into WGS84.  Routines are provided to do both
of these operations, but they are only approximate.  The inaccuracy of the
approximation varies according to where you are in the country but may be as
much as several metres in some areas.

To get more accurate results you need to combine this module with its
companion L<Geo::Coordinates::OSTN02> which implements the transformation
that now defines the relationship between GPS survey data based on WGS84 and
the British National Grid.  Using this module you should be able to get
results that are accurate to within a few centimetres, but it is slightly
slower and requires more memory to run.

Note that the OSGB (and therefore this module) does not cover the whole of
the British Isles, nor even the whole of the UK, in particular it covers
neither the Channel Islands nor Northern Ireland.  The coverage that is
included is essentially the same as the coverage provided by the OSGB
"Landranger" 1:50000 series maps.

=head1 SYNOPSIS

  use Geo::Coordinates::OSGB qw(ll_to_grid grid_to_ll);

  # Basic conversion routines
  ($easting,$northing) = ll_to_grid($lat,$lon);
  ($lat,$lon) = grid_to_ll($easting,$northing);

=head1 DESCRIPTION

These modules provide a collection of routines to convert between coordinates expressed as latitude & longtitude and
map grid references, using the formulae given in the British Ordnance Survey's excellent information leaflet, referenced
below in L<Theory>.  There are some key concepts explained in that section that you need to know in order to use these
modules successfully, so you are recommended to at least skim through it now.

The module is implemented purely in Perl, and should run on any Perl platform.

In this description `OS' means `the Ordnance Survey of Great Britain': the British
government agency that produces the standard maps of England, Wales, and
Scotland.  Any mention of `sheets' or `maps' refers to one or more of the 204
sheets in the 1:50,000 scale `Landranger' series of OS maps.

This code is fine tuned to the British national grid system.  You could use it
elsewhere but you would need to adapt it.  Some starting points for doing this
are explained in the L<Theory> section below.


=head1 FUNCTIONS

The following functions can be exported from the C<Geo::Coordinates::OSGB>
module:

	grid_to_ll				  ll_to_grid

	shift_ll_into_WGS84		 shift_ll_from_WGS84

	parse_grid
	parse_trad_grid			 format_grid_trad
	parse_GPS_grid			  format_grid_GPS
	parse_landranger_grid	   format_grid_landranger

	parse_ISO_ll				format_ll_trad
								format_ll_ISO

None of these is exported by default, so pick the ones you want or use an C<:all> tag to import them all at once.

  use Geo::Coordinates::OSGB ':all';

=over 4

=item ll_to_grid(lat,lon)

When called in a void context, or with no arguments C<ll_to_grid> does nothing.

When called in a list context, C<ll_to_grid> returns two numbers that represent
the easting and the northing corresponding to the latitude and longitude
supplied.

The parameters can be supplied as real numbers representing decimal degrees, like this

	my ($e,$n) = ll_to_grid(51.5, 2.1);

Following the normal convention, positive numbers mean North or East, negative South or West.
If you have data with degrees, minutes and seconds, you can convert them to decimals like this:

	my ($e,$n) = ll_to_grid(51+25/60, 0-5/60-2/3600);

Or you can use a single string in ISO 6709 form, like this:

	my ($e,$n) = ll_to_grid('+5130-00005/');

To learn exactly what is matched by this last option, read the source of the module and look for the
definition of C<$ISO_LL_Pattern>.  Note that the neither the C<+> or C<-> signs at the
beginning and in the middle, nor the trailing C</> may be omitted.

If you have trouble remembering the order of the arguments, note that
latitude comes before longitude in the alphabet too.

The easting and northing will be returned as a whole number of metres from
the point of origin of the British Grid (which is a point a little way to the
south-west of the Scilly Isles).

If you want the result presented in a more traditional grid reference format you should pass the results to one of the
grid formatting routines, which are described below.  Like this.

	$gridref = format_grid_trad(ll_to_grid(51.5,-0.0833));
	$gridref = format_grid_GPS(ll_to_grid(51.5,-0.0833));
	$gridref = format_grid_landranger(ll_to_grid(51.5,-0.0833));

However if you call C<ll_to_grid> in a scalar context, it will
automatically call C<format_grid_trad> for you.

It is not needed for any normal work, but C<ll_to_grid()> also takes an
optional argument that sets the ellipsoid model to use.  This normally
defaults to `OSGB36', the name of the normal model for working with British
maps.  If you are working with the highly accurate OSTN02 conversions
supplied in the companion module in this distribution, then you will need to
produce pseudo-grid references as input to those routines.  For these
purposes you should call C<ll_to_grid()> like this:

	my $pseudo_gridref = ll_to_grid(51.2, -0.4, 'WGS84');

and then transform this to a real grid reference using C<ETRS89_to_OSGB36()>
from the companion module.  This is explained in more detail below.

=item format_grid_trad(e,n)

Formats an (easting, northing) pair into traditional `full national grid
reference' with two letters and two sets of three numbers, like this
`TQ 102 606'.  If you want to remove the spaces, just apply C<s/\s//g> to it.

	$gridref = format_grid_trad(533000, 180000); # TQ 330 800
	$gridref =~ s/\s//g;						 # TQ330800

If you want the individual components call it in a list context.

	($sq, $e, $n) = format_grid_trad(533000, 180000); # (TQ,330,800)

Note the easting and northing are truncated to hectometers (as the OS system
demands), so the grid reference refers to the lower left corner of the
relevant 100m square.

=item format_grid_GPS(e,n)

Users who have bought a GPS receiver may initially have been puzzled by the
unfamiliar format used to present coordinates in the British national grid format.
On my Garmin Legend C it shows this sort of thing in the display.

	TQ 23918
   bng 00972

and in the track logs the references look like this C<TQ 23918 00972>.

These are just the same as the references described on the OS sheets, except
that the units are metres rather than hectometres, so you get five digits in
each of the easting and northings instead of three.  So in a scalar context
C<format_grid_GPS()> returns a string like this:

	$gridref = format_grid_GPS(533000, 180000); # TQ 33000 80000

If you call it in a list context, you will get a list of square, easting, and northing, with the easting and northing as
metres within the grid square.

	($sq, $e, $n) = format_grid_GPS(533000, 180000); # (TQ,33000,80000)

Note that, at least until WAAS is working in Europe, the results from your
GPS are unlikely to be more accurate than plus or minus 5m even with perfect
reception.  Most GPS devices can display the accuracy of the current fix you
are getting, but you should be aware that all normal consumer-level GPS
devices can only ever produce an approximation of an OS grid reference, no
matter what level of accuracy they may display.  The reasons for this are
discussed below in the section on L<Theory>.

=item format_grid_landranger(e,n)

This routine does the same as C<format_grid_trad>, but it appends the number of the relevant OS Landranger 1:50,000
scale map to the traditional grid reference.  Note that there may be several or no sheets returned.  This is because
many (most) of the Landranger sheets overlap, and many other valid grid references are not on any of the sheets (because
they are in the sea or a remote island.  This module does not yet cope with the detached insets on some sheets.

In a list context you will get back a list like this:  (square, easting,
northing, sheet) or (square, easting, northing, sheet1, sheet2) etc.  There
are a few places where three sheets overlap, and one corner of Herefordshire
which appears on four maps (sheets 137, 138, 148, and 149).  If the GR is not
on any sheet, then the list of sheets will be empty.

In a scalar context you will get back the same information in a helpful
string form like this "NN 241 738 on OS Sheet 44".  Note that the easting and
northing will have been truncated to the normal hectometre three
digit form.  The idea is that you'll use this form for people who might actually
want to look up the grid reference on the given map sheet, and the traditional
GR form is quite enough accuracy for that purpose.

=item parse_trad_grid(grid_ref)

Turns a traditional grid reference into a full easting and northing pair in
metres from the point of origin.  The I<grid_ref> can be a string like
C<'TQ203604'> or C<'SW 452 004'>, or a list like this C<('TV', '435904')> or a list
like this C<('NN', '345', '208')>.


=item parse_GPS_grid(grid_ref)

Does the same as C<parse_trad_grid> but is looking for five digit numbers
like C<'SW 45202 00421'>, or a list like this C<('NN', '34592', '20804')>.

=item parse_landranger_grid(sheet, e, n)

This converts an OS Landranger sheet number and a local grid reference
into a full easting and northing pair in metres from the point of origin.

The OS Landranger sheet number should be between 1 and 204 inclusive (but
I may extend this when I support insets).  You can supply C<(e,n)> as 3-digit
hectometre numbers or 5-digit metre numbers.  In either case if you supply
any leading zeros you should 'quote' the numbers to stop Perl thinking that
they are octal constants.

This module will croak at you if you give it an undefined sheet number, or
if the grid reference that you supply does not exist on the sheet.

In order to get just the coordinates of the SW corner of the sheet, just call
it with the sheet number.  It is easy to work out the coordinates of the
other corners, because all OS Landranger maps cover a 40km square (if you
don't count insets or the occasional sheet that includes extra details
outside the formal margin).

=item parse_grid(grid_ref)

Attempts to match a grid reference some form or other
in the input string and will then call the appropriate grid
parsing routine from those defined above.  In particular it will parse strings in the form
C<'176-345210'> meaning grid ref 345 210 on sheet 176, as well as C<'TQ345210'> and C<'TQ 34500 21000'> etc.

=item grid_to_ll(e,n) or grid_to_ll(grid_ref)

When called in list context C<grid_to_ll()> returns a pair of numbers
representing longitude and latitude coordinates, as real numbers.  Following
convention, positive numbers are North and East, negative numbers are South
and West.  The fractional parts of the results represent fractions of
degrees.

When called in scalar context it returns a string in ISO longitude and latitude
form, such as C<'+5025-00403/'> with the result rounded to the nearest minute (the
formulae are not much more accurate than this).  In a void context it does
nothing.

The arguments must be an (easting, northing) pair representing the absolute grid reference in metres from the point of
origin.  You can get these from a grid reference string by calling C<parse_grid()> first.

An optional last argument defines the geoid model to use just as it does for C<ll_to_grid()>.  This is only necessary is
you are working with the pseudo-grid references produced by the OSTN02 routines.  See L<Theory> for more discussion.

=item format_ll_trad(lat, lon)

Takes latitude and longitude in decimal degrees as arguments and returns a string like this

	N52:12:34 W002:30:27

In a list context it returns all 8 elements (hemisphere, degrees, minutes, seconds for each of lat and lon) in a list.
In a void context it does nothing.

=item format_ll_ISO(lat, lon)

Takes latitude and longitude in decimal degrees as arguments and returns a string like this

	+5212-00230/

In a list context it returns all 6 elements (sign, degrees, minutes for each of lat and lon) in a list.
In a void context it does nothing.


=item parse_ISO_ll(ISO_string)

Reads an ISO 6709 formatted location identifier string such as '+5212-00230/'.
To learn exactly what is matched by this last option, read the source of the module and look for the
definition of C<$ISO_LL_Pattern>.  Note that the neither the C<+> or C<-> signs at the
beginning and in the middle, nor the trailing C</> may be omitted.  These strings can also include
the altitude of a point, in metres, like this: '+5212-00230+140/'.  If you omit the altitude, 0 is assumed.

In a list context it returns ($lat, $lon, $altitude).  So if you don't want or don't need the altitude, you should just
drop it, for example like this:

   my ($lat, $lon) = parse_ISO_ll('+5212-00230/')

In normal use you won't notice this.  In particular you don't need to worry about it when
passing the results on to C<ll_to_grid>, as that routine looks for an optional altitude after the lat/lon.

=item shift_ll_from_WGS84(lat, lon, altitude)

Takes latitude and longitude in decimal degrees (plus an optional altitude in metres) from a WGS84 source (such as your
GPS handset or Google Earth) and returns an approximate equivalent latitude and longitude according to the OSGM02 model.
To determine the OSGB grid reference for given WGS84 lat/lon coordinates, you should call this before
you call C<ll_to_grid>.  Like so:

  ($lat, $lon, $alt) = shift_ll_from_WGS84($lat, $lon, $alt);
  ($e, $n) = ll_to_grid($lat,$lon);

You don't need to call this to determine a grid reference from lat/lon coordinates printed on OSGB maps (the so called
"graticule intersections" marked in pale blue on the Landranger series).

This routine provide a fast approximation; for a slower, more accurate approximation use the companion
L<Geo::Coordinates::OSTN02> modules.

=item shift_ll_into_WGS84(lat, lon, altitude)

Takes latitude and longitude in decimal degrees (plus an optional altitude in metres) from an OSGB source (such as
coordinates you read from a Landranger map, or more likely coordinates returned from C<grid_to_ll()>) and adjusts them
to fit the WGS84 model.

To determine WGS84 lat/lon coordinates (for use in Wikipedia, or Google Earth etc) for a given OSGB grid
reference, you should call this after you call C<grid_to_ll()>.  Like so:

  ($lat, $lon) = grid_to_ll($e, $n);
  ($lat, $lon, $alt) = shift_ll_into_WGS84($lat, $lon, $alt);

This routine provide a fast approximation; for a slower, more accurate approximation use the companion
L<Geo::Coordinates::OSTN02> modules.

=back

=head1 THEORY

The algorithms and theory for these conversion routines are all from
I<A Guide to Coordinate Systems in Great Britain>
published by the OSGB, April 1999 and available at
http://www.ordnancesurvey.co.uk/oswebsite/gps/information/index.html

You may also like to read some of the other introductory material there.
Should you be hoping to adapt this code to your own custom Mercator
projection, you will find the paper called I<Surveying with the
National GPS Network>, especially useful.

The routines are intended for use in Britain with the Ordnance Survey's
National Grid, however they are written in an entirely generic way, so that
you could adapt them to any other ellipsoid model that is suitable for your
local area of the earth.   There are other modules that already do this that
may be more suitable (which are referenced in the L<See Also> section), but
the key parameters are all defined at the top of the module.

	$ellipsoid_shapes{OSGB36} = [ 6377563.396,  6356256.910  ];
	use constant LAM0 => RAD * -2;  # lon of grid origin
	use constant PHI0 => RAD * 49;  # lat of grid origin
	use constant E0   =>  400000;   # Easting for origin
	use constant N0   => -100000;   # Northing for origin
	use constant F0   => 0.9996012717; # Convergence factor

The ellipsoid model is defined by two numbers that represent the major and
minor radius measured in metres.  The Mercator grid projection is then
defined by the other five parameters.  The general idea is that you pick a
suitable point to start the grid that minimizes the inevitable distortion
that is involved in a Mercator projection from spherical to Euclidean
coordinates.  Such a point should be on a meridian that bisects the area of
interest and is nearer to the equator than the whole area.  So for Britain
the point of origin is 2degW and 49degN (in the OSGB geoid model) which is near
the Channel Islands.  This point should be set as the C<LAM0> and C<PHI0>
parameters (as above) measured in radians.  Having this True Point of Origin
in the middle and below (or above if you are antipodean) minimizes
distortion but means that some of the grid values would be negative unless
you then also adjust the grid to make sure you do not get any negative
values in normal use.  This is done by defining the grid coordinates of the
True Point of Origin to be such that all the coordinates in the area of
interest will be positive.  These are the parameters C<E0> and C<N0>.
For Britain the coordinates are set as 400000 and -100000, so the that point
(0,0) in the grid is just to the south west of the Scilly Isles.  This (0,0)
point is called the False Point of Origin.  The fifth parameter affects the
convergence of the Mercator projection as you get nearer the pole; this is
another feature designed to minimize distortion, and if in doubt set it to
1 (which means it has no effect).  For Britain, being so northerly it is set
to slightly less than 1.

=head2 The British National Grid

One consequence of the True Point of Origin of the British Grid being set to
C<+4900-00200/> is that all the vertical grid lines are parallel to the 2degW
meridian; you can see this on the appropriate OS maps (for example
Landranger sheet 184), or on the C<plotmaps.pdf> picture supplied with this
package.  The effect of moving the False Point of Origin to the far south
west is that all grid references always positive.

Strictly grid references are given as whole numbers of metres from this
point, with the easting always given before the northing.  For everyday use
however, the OSGB suggest that grid references need only to be given within
the local 100km square as this makes the numbers smaller.  For this purpose
they divide Britain into a series of 100km squares identified in pair of
letters:  TQ, SU, ND, etc.  The grid of the big squares actually used is
something like this:

							   HP
							   HU
							HY
				   NA NB NC ND
				   NF NG NH NJ NK
				   NL NM NN NO NP
					  NR NS NT NU
					  NW NX NY NZ
						 SC SD SE TA
						 SH SJ SK TF TG
					  SM SN SO SP TL TM
					  SR SS ST SU TQ TR
				   SV SW SX SY SZ TV

SW covers most of Cornwall, TQ London, and HU the Shetlands.  Note that it
has the neat feature that N and S are directly above each other, so that
most Sx squares are in the south and most Nx squares are in the north.

Within each of these large squares, we only need five digit coordinates ---
from (0,0) to (99999,99999) --- to refer to a given square metre.  For daily
use however we don't generally need such precision, so the normal
recommended usage is to use units of 100m (hectometres) so that we only need
three digits for each easting and northing --- 000,000 to 999,999.  If we
combine the easting and northing we get the familiar traditional six figure
grid reference.  Each of these grid references is repeated in each of the
large 100km squares but for local use with a particular map, this does not
usually matter.  Where it does matter, the OS suggest that the six figure
reference is prefixed with the identifier of the large grid square to give a
`full national grid reference', such as TQ330800.  This system is described
in the notes of in the corner of every Landranger 1:50,000 scale map.

Modern GPS receivers can all display coordinates in the OS grid system.  You
just need to set the display units to be `British National Grid' or whatever
similar name is used on your unit.  Most units display the coordinates as two
groups of five digits and a grid square identifier.  The units are metres within
the grid square (although beware that the GPS fix is unlikely to be accurate down
to the last metre).


=head2 Geoid models

This section explains the fundamental problems of mapping a spherical earth
onto a flat piece of paper (or computer screen).  A basic understanding of
this material will help you use these routines more effectively.  It will
also provide you with a good store of ammunition if you ever get into an
argument with someone from the Flat Earth Society.

It is a direct consequence of Newton's law of universal gravitation (and in
particular the bit that states that the gravitational attraction between two
objects varies inversely as the square of the distance between them) that all
planets are roughly spherical.  (If they were any other shape gravity would
tend to pull them into a sphere).  On the other hand, most useful surfaces
for displaying large scale maps (such as pieces of paper or screens) are
flat.  There is therefore a fundamental problem in making any maps of the
earth that its curved surface being mapped must be distorted at least
slightly in order to get it to fit onto the flat map.

This module sets out to solve the corresponding problem of converting
latitude and longitude coordinates (designed for a spherical surface) to and
from a rectangular grid (for a flat surface).  This projection is in itself
is a fairly lengthy bit of maths, but what makes it extra complicated is
that the earth is not quite a sphere.  Because our planet spins about a
vertical axis, it tends to bulge out slightly in the middle, so it is more
of an oblate spheroid than a sphere.  This makes the maths even longer, but
the real problem is that the earth is not a regular oblate spheroid either,
but an irregular lump that closely resembles an oblate spheroid and which is
constantly (if slowly) being rearranged by plate tectonics.  So the best we
can do is to pick an imaginary regular oblate spheroid that provides a good
fit for the region of the earth that we are interested in mapping.  The
British Ordnance Survey did this back in 1830 and have used it ever since as
the base on which the National Grid for Great Britain is constructed.  You
can also call an oblate spheroid an ellipsoid if you like.  The general term
for an ellipsoid model of the earth is a "geoid".

The first standard OSGB geoid is known as "Airy 1830" after the year of its
first development.  It was revised in 1936, and that version, generally
known as OSGB36, is the basis of all current OSGB mapping.  In 2002 the
model was redefined (but not functionally changed) as a transformation from
the international geoid model WGS84.  This redefinition is called OSGM02.
For the purposes of these modules (and most other purposes) OSGB36 and
OSGM02 may be treated as synonyms.

The general idea is that you can establish your latitude and longitude by
careful observation of the sun, the moon, the planets, or your GPS handset,
and that you then do some clever maths to work out the corresponding grid
reference using a suitable geoid.  These modules let you do the clever
maths, and the geoid they use is the OSGM02 one.  This model provides a good
match to the local shape of the Earth in the British Isles, but is not
designed for use in the rest of the world; there are many other models in
use in other countries.

In the mid-1980s a new standard geoid model was defined to use with the
fledgling global positioning system (GPS).  This model is known as WGS84, and
is designed to be a compromise model that works equally well for all parts of
the globe (or equally poorly depending on your point of view --- for one
thing WGS84 defines the Greenwich observatory in London to be not quite on
the 0deg meridian).  Nevertheless WGS84 has grown in importance as GPS systems
have become consumer items and useful global mapping tools (such as Google
Earth) have become freely available through the Internet.  Most latitude and
longitude coordinates quoted on the Internet (for example in Wikipedia) are
WGS84 coordinates.

One thing that should be clear from the theory is that there is no such
thing as a single definitive set of coordinates for every unique spot on
earth.  There are only approximations based on one or other of the accepted
geoid models, however for most practical purposes good approximations are
all you need.  In Europe the official definition of WGS84 is sometime
referred to as ETRS89.  For all practical purposes in Western Europe the OS
advise that one can regard ETRS89 as identical to WGS84 (unless you need to
worry about tectonic plate movements).

=head2 Practical implications

If you are working exclusively with British OS maps and you merely want
to convert from the grid to the latitude and longitude coordinates printed (as
faint blue crosses) on those maps, then all you need from these modules are
the plain C<grid_to_ll()> and C<ll_to_grid()> routines.  On the other hand if
you want to produce latitude and longitude coordinates suitable for Google
Earth or Wikipedia from a British grid reference, then you need an extra
step.  Convert your grid reference using C<grid_to_ll()> and then shift it
from the OSGB model to the WGS84 model using C<shift_ll_into_WGS84()>.  To
go the other way round, shift your WGS84 lat/lon coordinated into OSGB,
using C<shift_ll_from_WGS84()>, before you convert them using
C<ll_to_grid()>.

If you have a requirement for really accurate work (say to within a
millimetre or two) then you need to use the OS's transformation matrix
called OSTN02.  This monumental work published in 2002 re-defined the
British grid in terms of offsets from WGS84 to allow really accurate grid
references to be determined from really accurate GPS readings (the sort you
get from professional fixed base stations, not from your car's sat nav or
your hand-held device).  The problem with it is that it defines the grid in
terms of a deviation in three dimensions from a pseudo-grid based on WGS84
and it does this separately for every square km of the country, so the data
set is huge and takes a second or two to load even on a fast machine.
Nevertheless a Perl version of OSTN02 is included as a separate module in
this distribution just in case you really need it (but you don't need it for
any "normal" work).  Because of the way OSTN02 is defined, the sequence of
conversion and shifting works differently from the approximate routines
described above.

Starting with a really accurate lat/lon reading in WGS84 terms, you need to
transform it into a pseudo-grid reference using C<ll_to_grid()> using an
optional argument to tell it to use the WGS84 geoid parameters instead of
the default OSGB parameters.  The L<Geo::Coordinates::OSTN02> package
provides a routine called C<ETRS89_to_OSGB36()> which will shift this pseudo-grid
reference into an accurate OSGB grid reference.  To go back the other way,
you use C<OSGB36_to_ETRS89()> to make a pseudo-grid reference, and then call
C<grid_to_ll()> with the WGS84 parameter to get WGS84 lat/long coordinates.


   ($lat, $lon, $height) = (51.5, -1, 10);
   ($x, $y) = ll_to_grid($lat, $lon, 'WGS84');
   ($e, $n, $elevation) = ETRS89_to_OSGB36($x, $y, $height);

   ($x, $y, $z) = OSGB36_to_ETRS89($e, $n, $elevation);
   ($lat, $lon) = grid_to_ll($x, $y, 'WGS84');


=head1 EXAMPLES

  # to import everything try...
  use Geo::Coordinates::OSGB ':all';

  # Get full coordinates in metres from GR
  ($e,$n) = parse_trad_grid('TQ 234 098');

  # Latitude and longitude according to the OSGB geoid (as
  # printed on OS maps), if you want them to work in Google
  # Earth or some other tool that uses WGS84 then adjust results
  ($lat, $lon) = grid_to_ll($e, $n);
  ($lat, $lon, $alt) = shift_ll_into_WGS84($lat, $lon, $alt);
  # and to go the other way
  ($lat, $lon, $alt) = shift_ll_from_WGS84($lat, $lon, $alt);
  ($e, $n) = ll_to_grid($lat,$lon);
  # In both cases the elevation is in metres (default=0m)

  # Reading and writing grid references
  # Format full easting and northing into traditional formats
  $gr1 = format_grid_trad($e, $n);	  # "TQ 234 098"
  $gr1 =~ s/\s//g;					  # "TQ234098"
  $gr2 = format_grid_GPS($e, $n);	   # "TQ 23451 09893"
  $gr3 = format_grid_landranger($e, $n);# "TQ 234 098 on Sheet 176"
  # or call in list context to get the individual parts
  ($sq, $e, $n) = format_grid_trad($e, $n); # ('TQ', 234, 98)

  # parse routines to convert from these formats to full e,n
  ($e,$n) = parse_trad_grid('TQ 234 098');
  ($e,$n) = parse_trad_grid('TQ234098'); # spaces optional
  ($e,$n) = parse_trad_grid('TQ',234,98); # or even as a list
  ($e,$n) = parse_GPS_grid('TQ 23451 09893'); # as above..

  # You can also get grid refs from individual maps.
  # Sheet between 1..204; gre & grn must be 3 or 5 digits long
  ($e,$n) = parse_landranger_grid(176,123,994);

  # With just the sheet number you get GR for SW corner
  ($e,$n) = parse_landranger_grid(184);

  # Reading and writing lat/lon coordinates
  ($lat, $lon) = parse_ISO_ll("+52-002/");
  $iso = format_ll_ISO($lat,$lon);	# "+520000-0020000/"
  $str = format_ll_trad($lat,$lon);   # "N52:00:00 W002:00:00"


=head1 BUGS

The conversions are only approximate.   So after

  ($a1,$b1) = grid_to_ll(ll_to_grid($a,$b));

neither C<$a==$a1> nor C<$b==$b1>. However C<abs($a-$a1)> and C<abs($b-$b1)>
should be less than C<0.00001> which will give you accuracy to within a
metre. In the middle of the grid 0.00001 degrees is approximately 1 metre.
Note that the error increases the further away you are from the
central meridian of the grid system.

The C<format_grid_landranger()> does not take account of inset areas on the
sheets.  So if you feed it a reference for the Scilly Isles, it will tell you
that the reference is not on any Landranger sheet, whereas in fact the
Scilly Isles are on an inset in the SW corner of Sheet 203.  There is
nothing in the design that prevents me adding the insets, they just need to
be added as extra sheets with names like "Sheet 2003 Inset 1" with their own
reference points and special sheet sizes.  Collecting the data is another
matter.

Not enough testing has been done.  I am always grateful for the feedback I
get from users, but especially for problem reports that help me to make this
a better module.


=head1 AUTHOR

Toby Thurston ---  6 Nov 2008

toby@cpan.org

=head1 SEE ALSO

The UK Ordnance Survey's theory paper referenced above in L<Theory>.

See L<Geo::Coordinates::Convert> for a general approach (not based on the above
paper).

See L<Geo::Coordinates::Lambert> for a French approach.

=cut
'''
