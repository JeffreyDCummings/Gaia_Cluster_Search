# Gaia_Cluster_Search
Purpose: The massive Gaia DR2 database has a 5-parameter astrometric solution for over 1.3 billion stars, which includes RA and DEC,
parallax (distance), and 2D proper motion.  Additionally, G magnitudes and BP-RP colors are provided.  These data can be used to select 
out the coeval stellar populations, which will have the same position and velocity information, and can be differentiated from the field 
stars.  

This program uses the HDBSCAN clustering algorithm to select out multi-dimensional overdensities (e.g., star clusters), and it can handle 
clusters of variable density (unlike the similar DBSCAN) and allows the definition of minimum cluster size.  This means that the program 
does not cluster small spurious overdensities that can happen in big data sets and do not represent true clusters.  Testing this on 
multiple fields shows that it is very successful in identifying one to multiple clusters within a single field, while typically not 
clustering field stars.  However, the overdensity of field stars near proper motions of approximately 0,0 are commonly clustered, but they
produce one large cluster with a large proper motion dispersion that is very distinct from true clusters with narrow dispersions.
Additionally, these false star clusters tend to have large distance dispersions and be primarly composed of faint stars that do not create
a main sequence like structure in color and magnitude.  Hence, these false clusters can be filtered out in the final analysis.

This directory contains several files: The python programs: gaia_search.py, cluster_hdbscan.py, and cluster_plot.py, 
cluster_hdbscan_wdsearch.py, and several example Gaia field csv files: blanco1gaiafield2.5.csv, velaOB2gaiafield2.csv, and 
n2422.n2423gaiafield2.csv.

These files are described in more detail below:

gaia_search.py:
This script does a circular search of the Gaia database for a defined radius of the sky centered on either a known astronomical object or
a defined RA and DEC.  Note that large searches (r > 1.5 degrees) in dense fields are prone to time out the Gaia API call.  In these 
cases, use directly the Gaia DR2 archive (https://gea.esac.esa.int/archive/) and its ADQL search, which can output large searches as csv
files.  An example is shown below, where to input your desired coordinates and radius, change the CIRCLE call to 
"CIRCLE('ICRS',RA,DEC,radius), where RA, DEC, and radius are all in decimal degrees.  Similarly, edit the parallax and parallax error 
ratio cut as desired:

SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp
FROM gaiadr2.gaia_source 
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',58.833333333333336,-20.55138888888889,2.5))=1  AND  (gaiadr2.gaia_source.parallax>=1 AND gaiadr2.gaia_source.parallax_over_error>=3)

cluster_hdbscan.py:
With the input cluster field csv, the stellar parameters (RA, DEC, distance, and 2D proper motions) are scaled for a more appropriate 
dimension for clustering.  Additionally, RA is wrapped around if necessary.  These are input to the HDBSCAN method, and the clustering 
results and membership probabilities are created and added to the dataframe.  The selected cluster number's data is output to a csv file
for future analysis.  Lastly, this entire dataframe and whether the membership probability is desired to be displayed is input into our 
defined clusterplot method.  

cluster_plot.py:
This method plots the important output information for each cluster that passes the necessary proper motion distribution cut (e.g., Each 
dimensions PM IQR in physical space is less than 3 km/s), the distance distribution cut (distance IQR less than 300 pc), and a check that
at least not more than 80% of the clustered stars have M_G fainter than 10.  These clusters are color coded and plotted in proper motion
space, scaled RA and DEC space, distance, and M_G and bp-rp space.  If membership probabilities are requested, a color map is applied to 
each cluster's data points to represent these probabilities.  A different figure window is produced for each cluster and their members are
given larger data size relative to the other clusters/nonclustered data for clarity/comparison.  Lastly, if a 3D plot is requested, each
cluster has an additional plot window centered on its 3D spatial information, and all other nearby clusters are shown for reference.

cluster_hdbscan_wdsearch.py:
This program is a variant of cluster_hdbscan.py that performs that same analysis while also checking if clusters that are consistent with
star clusters contain any member sources with characteristics consistent with a white dwarf.  If yes, the absolute magnitude (corrected by
its own parallax-based distance rather than the median parallax of the cluster) and color are output to screen.  For nearby white dwarfs
(within 150 pc), its own parallax should be used, but for further analysis, more distant white dwarfs should consider correcting 
magnitudes based on their median cluster parallax.  Additionally, zero reddening and extinction are adopted in the output, but they may 
play an important role and should be considered before further analysis with this photometry.
