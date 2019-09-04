import pyvo as vo
import os.path
from os import path
import warnings

#These modules output a multiple Warnings that are currently unavoidable yet do not affect 
#the output.  Hence, I have filtered them out.

warnings.filterwarnings("ignore", module='astropy.io.votable.tree')
warnings.filterwarnings("ignore", module='astropy.io.votable.converters')
warnings.filterwarnings("ignore", module='astropy.extern.six')

#The result outputs a circular search of radius (degrees).  Note: Radii larger than 1.5 
#degrees typically time out and the DR2 website must be used directly.
radius=0.3

#Cut all targets from output that have a parallax (mas) below parallax_cut or a error ratio below parallax_error_ratio
parallax_cut=1
parallax_error_ratio=3

#The file name, which includes the radius to differentiate variable search sizes.
Filename="n2547gaiafield"+str(radius)+".csv"

#The cluster name that will be searched for in SIMBAD to find its central coordinates.
Objectname="NGC 2547"

#Function for converting RA and DEC in H:M:S (or H:M) to decimal degrees
def coordconv(RAobj,DECobj):
	RAset=0
	DECset=0
	for x in range(len(RAobj)):
		RAset+=float(RAobj[x])*15/(60**x)
	
	for y in range(len(DECobj)):
		if y == 0:
			DECset+=float(DECobj[y])
		else:
			if DECset >= 0:
				DECset+=float(DECobj[y])/60**y
			else:
				DECset-=float(DECobj[y])/60**y
	
	print("RA = "+str(round(RAset,4))+" DEC = "+str(round(DECset,4)))
	return RAset,DECset

#Checks if file already exists, if not, continues search for SIMBAD coordinates of Objectname and outputs it for conversion.
if not path.exists(Filename):
	from astroquery.simbad import Simbad
	
	result_table = Simbad.query_object(Objectname)

	RAobj=result_table["RA"][0].split()
	DECobj=result_table["DEC"][0].split()
	RAset,DECset=coordconv(RAobj,DECobj)
		
#If you want to define RA and DEC directly, uncomment them and define them here (degrees)
	#RAset=116
	#DECset=-38
	
#Initialize tap service
	tap_service = vo.dal.TAPService('https://gaia.aip.de/tap')
	
# optional: use your API token to use your account
	#vo.utils.http.session.headers['Authorization'] = 'Token input'
	
#Search Gaia database and return to dataframe (Filename) the columns of interest for a circular search of defined radius and parallax cuts around ObjectName.
	tap_result = tap_service.search("SELECT gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp FROM gdr2.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS',"+str(RAset)+","+str(DECset)+","+str(radius)+")) AND gaia_source.parallax>="+str(parallax_cut)+" AND gaia_source.parallax_over_error>="+str(parallax_error_ratio))
	    
	tap_result.to_table().to_pandas().to_csv(Filename)

