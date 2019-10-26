"""Clusters Gaia DR2 data into groups with open cluster characteristics.  Outputs to
screen objects consistent with white dwarf cluster members."""
import hdbscan
import scipy.stats as spy
import numpy as np
import pandas as pd
from cluster_plot import cluster_plot, iqr_calc_angle

# This selects which output cluster data the program will write to a csv file
# for subsequent analysis.  The star cluster of interest is commonly "0", but
# it can sometimes be a variety of numbers depending on the number of clusters
# in the field.
CLUSTER_EXTRACT_NUM = 0

# This flag determines whether or not to add a colormap representing each
# star's cluster membership probability to the final cluster plots.  Enter 1
# for yes and 0 (or anything that isn't 1) for no.  Plotting probabilities is
# recommended, but when probability information isn't critical and simpler
# easier to read figures are desired, consider not adding the information.
PLOT_MEMBERSHIP_PROB = 1

# CLUSTER_NAME identifies the input file (usually using some type of
# shorthand) and output cluster file names.  FIELD_RADIUS is the input
# GAIA DR2 search radius (in degrees), and is part of the input csv filename.
CLUSTER_NAME = "velaOB2"
FIELD_RADIUS = 2

MIN_SAMPLE = 64

def gaia_dr2_read_setup():
    """ Read the input csv file (from GaiaSearch.py) into a dataframe.  Return the
    calculated distance and further important parameters and drop stars that are
    missing these parameters.  """
    clusterfull = pd.read_csv(CLUSTER_NAME+"gaiafield"+str(FIELD_RADIUS)+".csv", delimiter=",")

# A second dataframe is created that only contains the parameters that will
# be used to base the HDBSCAN clustering on and for analyzing the results; stars
# that are missing a parameter are dropped.
    fieldparfull = clusterfull[["ra", "dec", "pmra", "pmdec", "parallax", \
        "phot_g_mean_mag", "bp_rp"]]
    fieldpar = fieldparfull.dropna().copy()

# A new distance (pc) column is created in the dataframe based on the observed
# parallax.  A 0.03 zeropoint correction is applied.
    fieldpar['distance'] = 1/((fieldpar['parallax']+0.03)*0.001)
    return fieldpar


def ra_wrapper(ra_dataframe):
    """ The input RA are checked to see if they wrap around RA=0,360 and are
    corrected.  In either case, the RA is also transformed to a uniform scale
    and added to the main dataframe and returned.  Note that this program may
    have some difficulty around the poles, but there are no clusters to analyze
    in those regions."""
    ramax = max(ra_dataframe["ra"])
    ramin = min(ra_dataframe["ra"])
    wrap_check = 0
    ra_dataframe["rawrapped"] = ra_dataframe["ra"]
    # If the maximum and minimum RA are more than 200 degrees apart, it is safe to
    # assume that this is because they span RA=0,360.  In this case, the minimum
    # of the large values (~355 deg) and the maximum of the small value (~5 deg)
    # must be calculated/corrected and then applied to properly scale and tie
    # together the RA.
    if ramax - ramin > 200:
        wrap_check = 1
        ramax = ra_dataframe.loc[ra_dataframe["ra"] < 200, ["ra"]].max().at["ra"]
        ramin = ra_dataframe.loc[ra_dataframe["ra"] > 200, ["ra"]].min().at["ra"]
        if ramin < (360 - ramax):
            ramax += 360
            ra_dataframe.loc[ra_dataframe['rawrapped'] >= 0, 'rawrapped'] = \
             ra_dataframe['rawrapped']+360
        else:
            ramin -= 360
            ra_dataframe.loc[ra_dataframe['rawrapped'] > 180, 'rawrapped'] = \
             ra_dataframe['rawrapped']-360
    racen = (ramin + ramax)/2
    ra_dataframe["ratransform"] = ra_dataframe["ra"] - racen
    if wrap_check == 1:
        ra_dataframe.loc[ra_dataframe['ratransform'] > 180, 'ratransform'] = \
         ra_dataframe['ratransform']-360
        ra_dataframe.loc[ra_dataframe['ratransform'] < -180, 'ratransform'] = \
         ra_dataframe['ratransform']+360
    return ra_dataframe


def parameter_scaler(fieldparscaled, fieldpar):
    """ Apply appropriate scales to normalize the parameters or place them in a
     uniform space.  The most appropriate scaling can be cluster dependent.  Most
     notably with increasing distance, the cluster distance and PM variations of
     true members increase due to increased errors. I recommend scaling these so
     that the output unscaled IQRs of the cluster members for each parameter are
     roughly equal by adjusted what fieldpar['distance'] is divided by and what
     fieldpar['pmra'] and ['pmdec'] are multiplied by. """

    fieldparscaled['parallax'] = fieldpar['distance']/5
    fieldparscaled['ratransform'] = fieldparscaled['ratransform']* \
        np.cos(fieldpar["dec"]*3.14159/180)
    fieldparscaled['ra'] = fieldparscaled['ratransform']*3.14159/180* \
        fieldpar['distance']
    fieldparscaled['dec'] = fieldpar['dec']*3.14159/180*fieldpar['distance']
    fieldparscaled['pmra'] = fieldpar['pmra']*10
    fieldparscaled['pmdec'] = fieldpar['pmdec']*10
    return fieldparscaled


def clustering_algorithm(fieldparscaled, fieldpar):
    """ The scaled parameters are applied to the hdbscan function. A minimum cluster
     and sample size are set.  64 is typically an appropriate min_samples for a
     5-dimensional hdbscan, but I recommend looking at potential differences in
     output when this parameter is varied by up to a factor of 2 (if not more).
     Especially adjust min_samples if a clear cluster in PM space is missed in
     the output.  A smaller min_samples increases noise but helps differentiate
     from the field the clusters with proper motions comparable to the field.
     Note that we can also use G magnitude and BP-RP color to increase clustering
     accuracy.  They are commented out in this input list, but they can help to
     create a cleaner cluster main sequence, BUT at the expense of removing many
     cluster giants/subdwarfs/white dwarfs. """

    clustering = hdbscan.HDBSCAN(min_cluster_size=50, min_samples=MIN_SAMPLE, cluster_selection_method="leaf").\
        fit(fieldparscaled[["ra", "dec", "pmra", "pmdec", "parallax"]])
        # ,"phot_g_mean_mag","bp_rp"]])

    # The cluster selection and probability results, plus transformed RA, are
    # added as new columns to the original dataframe.
    clusterselection = clustering.labels_
    clusterprobabilities = clustering.probabilities_
    fieldpar["clusternum"] = clusterselection.tolist()
    fieldpar["clusterprob"] = clusterprobabilities.tolist()
    fieldpar["ratransform"] = fieldparscaled["ratransform"]
    fieldpar["rawrapped"] = fieldparscaled["rawrapped"]

    # The stars absolute G (based on Gaia parallax) is calculated and the selected
    # CLUSTER_EXTRACT_NUM cluster is output to a csv file.
    fieldpar["M_G"] = fieldpar["phot_g_mean_mag"] - \
        np.log10(fieldpar["distance"]/10)*5
    fieldpar[fieldpar["clusternum"] == CLUSTER_EXTRACT_NUM]. \
        to_csv(CLUSTER_NAME+"clustering.csv")
    return fieldpar

def white_dwarf_identification(fieldpar):
    """ This photometrically selects out the information on cluster members consistent with white
    dwarf photometry and prints them to screen.  Note that the output photometry assumes zero
    reddening and adopts the white dwarf distance rather than the cluster distance. """
    for step in range(max(fieldpar.clusternum.unique())+1):
        distcen = np.median(fieldpar.loc[fieldpar["clusternum"] == step, ["distance"]])
        distiqr = spy.iqr(fieldpar.loc[fieldpar["clusternum"] == step, ["distance"]])
        pmraiqr = iqr_calc_angle(fieldpar.loc[fieldpar["clusternum"] == step], distcen,\
         "pmra")*0.271795
        pmdeciqr = iqr_calc_angle(fieldpar.loc[fieldpar["clusternum"] == step], distcen,\
         "pmdec")*0.271795
        true_cluster = []
        if (pmraiqr < 3 and pmdeciqr < 3 and distiqr < 500):
            true_cluster.append(step)
        for _, row in fieldpar.iterrows():
            if row["clusternum"] in true_cluster and row["bp_rp"] < 0.25 and row["M_G"] > 9 and\
             row["M_G"] > row["bp_rp"] * 5.556 + 10.111:
                print("WD in Cluster "+str(int(row["clusternum"]))+": Membership Prob. = "\
                 +str(round(row["clusterprob"], 3))+", M_G = "+str(round(row["M_G"],\
                 3))+" and BP-RP = "+str(round(row["bp_rp"], 3)))

def main():
    """ Main program series, which calls data and clustering functions and then plots
    them """
    fieldpar = gaia_dr2_read_setup()
    fieldparscaled = fieldpar.copy()
    fieldparscaled = ra_wrapper(fieldparscaled)
    fieldparscaled = parameter_scaler(fieldparscaled, fieldpar)
    fieldpar = clustering_algorithm(fieldparscaled, fieldpar)
    white_dwarf_identification(fieldpar)
    # Plots the output clustered data.  See clusterplot.py for details.
    cluster_plot(fieldpar, PLOT_MEMBERSHIP_PROB)

main()
