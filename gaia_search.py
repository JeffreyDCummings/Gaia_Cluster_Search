""" Script for calling Gaia DR2 API and outputs searched parameters results surrounding
 either an input cluster name or a set RA and DEC position."""
import warnings
import os.path
import pyvo as vo

# These modules output a multiple Warnings that are currently unavoidable yet do
# not affect the output.  Hence, I have filtered them out.
warnings.filterwarnings("ignore", module='astropy.io.votable.tree')
warnings.filterwarnings("ignore", module='astropy.io.votable.converters')
warnings.filterwarnings("ignore", module='astropy.extern.six')

# The result outputs a circular search of radius (degrees).  Note: Radii larger
# than 1.5 degrees typically will time out and the DR2 website must be used
# directly with ADQL and selecting all of the appropriate columns used in the
# hdbscan analysis.
FIELD_RADIUS = 0.15

# The cluster name that will be searched for in SIMBAD to find its central coordinates.
CLUSTER_NAME_SEARCH = "Collinder173"
CLUSTER_NAME = "c173"

# Cut all targets from output that have a parallax (mas) below PARALLAX_CUT or a error
# ratio below PARALLAX_ERROR_RATIO
PARALLAX_CUT = 1
PARALLAX_ERROR_RATIO = 3

# The file name, which includes the radius to differentiate variable search sizes.
FILENAME = CLUSTER_NAME+"gaiafield"+str(FIELD_RADIUS)+".csv"


def coordconv(ra_input, dec_input):
    """ Takes input RA and DEC from H:M:S (or H:M) and returns in decimal degrees. """
    ra_output = 0
    dec_output = 0
    for order, component in enumerate(ra_input):
        ra_output += float(component)*15/(60**order)

    for order, component in enumerate(dec_input):
        if order == 0:
            dec_output += float(component)
        else:
            if dec_output >= 0:
                dec_output += float(component)/60**order
            else:
                dec_output -= float(component)/60**order

    print("RA = "+str(round(ra_output, 4))+" DEC = "+str(round(dec_output, 4)))
    return ra_output, dec_output


def gaia_call(ra_obj, dec_obj):
    """ Initilize and search Gaia DR2 database with API based on input parameters
     and outputs csv file. """

# Initialize tap service.
    tap_service = vo.dal.TAPService('https://gaia.aip.de/tap')

# optional: Use your API token to use your account.
    #vo.utils.http.session.headers['Authorization'] = 'Token input'

# Search Gaia database and return to csv (FILENAME) the columns of interest for a
# circular search of defined radius and parallax cuts around CLUSTER_NAME_SEARCH.
    tap_result = tap_service.search("SELECT gaia_source.source_id, gaia_source.ra,\
    gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,\
    gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.\
    phot_g_mean_mag,gaia_source.bp_rp FROM gdr2.gaia_source WHERE 1=CONTAINS(POINT('ICRS',\
    ra, dec),CIRCLE('ICRS',"+str(ra_obj)+","+str(dec_obj)+","+str(FIELD_RADIUS)+")) AND\
    gaia_source.parallax>="+str(PARALLAX_CUT)+" AND\
    gaia_source.parallax_over_error>="+str(PARALLAX_ERROR_RATIO))

    tap_result.to_table().to_pandas().to_csv(FILENAME)


def main():
    """ Checks if file already exists.  If not, calls SIMBAD API for object coordinates
     and Gaia DR2 API for output parameters. """

# Checks if file already exists, if not, continues search for SIMBAD coordinates
# of CLUSTER_NAME_SEARCH and outputs it for conversion.
    if not os.path.exists(FILENAME):
        from astroquery.simbad import Simbad

        result_table = Simbad.query_object(CLUSTER_NAME_SEARCH)

        ra_obj, dec_obj = coordconv(result_table["RA"][0].split(),\
        result_table["DEC"][0].split())

# If you want to define RA and DEC directly, uncomment them and define them here (degrees)
        #ra_obj = 116
        #dec_obj = -38

        gaia_call(ra_obj, dec_obj)

main()
