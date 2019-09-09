import hdbscan
import numpy as np
import pandas as pd
from clusterplot import clusterplot

#This selects which output cluster data the program will write to a csv file for subsequent analysis.  The star cluster of interest is 
#commonly "0", but it can sometimes be a variety of numbers depending on the number of clusters in the field.
extractnum=0

#This flag determines whether or not to add a colormap representing each star's cluster membership probability to the final cluster plots.
#Enter 1 for yes and 0 (or anything that isn't 1) for no.  Plotting probabilities is recommended, but when probability information
#isn't critical and simpler/easier to read figures are desired, consider not adding the information.
PlotMembershipProb=1

#The cluster name identifies the input file (usually using some type of shorthand) and output cluster file names.
clustername="velaOB2"

#Read the input csv file (from GaiaSearch.py) into a dataframe.  The format of the name is clustername+"gaiafield"+radius size of Gaia 
#search in degrees+".csv".
clusterfull=pd.read_csv(clustername+"gaiafield2.csv",delimiter=",")

#A second dataframe is created that only contains the parameters that will be used to base the HDBSCAN clustering on and for analyzing the 
#results; stars that are missing a parameter are dropped.  
fieldparfull=clusterfull[["ra","dec","pmra","pmdec","parallax","phot_g_mean_mag","bp_rp"]]
fieldpar=fieldparfull.dropna().copy()

#A new distance (pc) column is created in the dataframe based on the observed parallax.  A 0.03 zeropoint correction is applied.
fieldpar['distance']=1/((fieldpar['parallax']+0.03)*0.001)

#A deep copy of the dataframe is created, I then apply several checks/transformations to wrap around the RAs together in the case 
#where they fall around RA=0,360.  Note that this program may have some difficulty around the poles, but there are no clusters to analyze
#in those regions.
fieldparscaled = fieldpar.copy()
ramax=max(fieldpar["ra"])
ramin=min(fieldpar["ra"])
if ramax - ramin > 200:
	ramax=fieldpar.loc[fieldpar["ra"]<200,["ra"]].max().at["ra"]
	ramin=fieldpar.loc[fieldpar["ra"]>200,["ra"]].min().at["ra"]
	if ramin < 360-ramax:
		ramin+=360
	else:
		ramax-=360
racen=(ramin+ramax)/2
fieldparscaled["ratransform"]=fieldpar["ra"]-racen
fieldparscaled.loc[fieldparscaled['ratransform'] > 300, 'ratransform'] = fieldparscaled['ratransform']-360
fieldparscaled.loc[fieldparscaled['ratransform'] < -300, 'ratransform'] = fieldparscaled['ratransform']+360

#We apply appropriate scales to normalize the parameters or place them in a uniform space.
#The most appropriate scaling can be cluster dependent.  Most notably with increasing distance, the cluster distance and PM variations 
#of true members increase due to increased errors. I recommend scaling these so that the output unscaled IQRs of the cluster members 
#for each parameter are roughly equal.  
fieldparscaled['parallax'] = fieldpar['distance']/5
fieldparscaled['ratransform']=fieldparscaled['ratransform']*np.cos(fieldparscaled["dec"]*3.14159/180)
fieldparscaled['ra']=fieldparscaled['ratransform']*3.14159/180*fieldpar['distance']
fieldparscaled['dec'] = fieldpar['dec']*3.14159/180*fieldpar['distance']
fieldparscaled['pmra'] = fieldpar['pmra']*10
fieldparscaled['pmdec'] = fieldpar['pmdec']*10

#The scaled parameters are applied to the hdbscan function.  A minimum cluster and sample size are set.  64 is typically an appropriate 
#min_samples for a 5-dimensional hdbscan, but I recommend looking at potential differences in output when this parameter is varied by up 
#to a factor of 2 (if not more).  Especially adjust min_samples if a clear cluster in PM space is missed in the output.  A smaller 
#min_samples increases noise but helps differentiate from the field the clusters with proper motions comparable to the field.  Note that 
#we can also use G magnitude and BP-RP color to increase clustering accuracy.  They are commented out in this input list, but they can 
#help to create a cleaner cluster main sequence, BUT at the expense of removing many cluster giants/subdwarfs/white dwarfs.
clustering = hdbscan.HDBSCAN(min_cluster_size=100, min_samples=64).fit(fieldparscaled[["ra","dec","pmra","pmdec","parallax"]]) #,"phot_g_mean_mag","bp_rp"]])

#The cluster selection and probability results, plus transformed RA, are added as new columns to the original dataframe.
clusterselection=clustering.labels_
clusterprobabilities=clustering.probabilities_
fieldpar["clusternum"]=clusterselection.tolist()
fieldpar["clusterprob"]=clusterprobabilities.tolist()
fieldpar["ratransform"]=fieldparscaled["ratransform"]

#The stars absolute G (based on Gaia parallax) is calculated and the selected "extractnum" cluster is output to a csv file.
fieldpar["M_G"]= fieldpar["phot_g_mean_mag"]-np.log10(fieldpar["distance"]/10)*5
fieldpar[fieldpar["clusternum"]==extractnum].to_csv(clustername+"clustering.csv")

#Plots the output clustered data.  See clusterplot.py for details.
clusterplot(fieldpar,PlotMembershipProb)
