def clusterplot(ClusteredData):
	import seaborn as sns
	import matplotlib.pyplot as plt
	import scipy.stats as spy
	import numpy as np

#This loop iterates through all clustered groups found with hdbscan and plots the requested clusters that pass the checks.
	for step in range(max(ClusteredData["clusternum"])+1):

#The median distance is found, which is used to transform proper motions and RA and DEC into spatial scales.
#Because these distributions are not always Gaussian, the interquartile ranges are used to analyze parameter variations throughout a cluster.
		distcen=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["distance"]])
		pmraiqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmra"]])*3.14159/180*distcen*0.271795
		pmdeciqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmdec"]])*3.14159/180*distcen*0.271795

#Hdbscan will typically create a cluster from the overdensity in the field disk stars, producing apparent clusters but with large PM variation.
#Therefore, we require that the transformed PM variation be less than 3 km/s in each axis to be displayed as a valid cluster.  If a legitimate 
#cluster is identified but retains a lot of noise/non-members, 3 km/s may remove it from display.  In these cases, you may increase the cut to 
#4 or 5 km/s, etc., to first analyze its output, but in the end I recommend adjusting the hdbscan min_samples or the parameter scalings.
		if (pmraiqr < 3 and pmdeciqr < 3):

#To set up the plots, we set a separate Figure for each cluster.  We set up the color palette, set a color to each group, 
#and set all clustered data to large size and nonclustered data to small data size.  Lastly, set tick label font size to 12 for clarity.
			plt.figure(step,figsize=(12, 10))
			palette = sns.color_palette('bright', np.unique(ClusteredData["clusternum"]).max() + 1)
			colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in ClusteredData["clusternum"]]
			sizes = [4 if x == step else 0.005 for x in ClusteredData["clusternum"]]
			plt.rcParams['xtick.labelsize'] = 12
			plt.rcParams['ytick.labelsize'] = 12

#A subplot illustrating the cluster distance histogram.  Additionally, the median cluster distance and distance IQR are plotted. 
			plt.subplot(223)
			plt.xlabel("Distance (pc)",fontsize=15)
			plt.ylabel("Number",fontsize=15)
			plt.xlim(distcen-250,distcen+250)
			bins = np.linspace(distcen-300,distcen+300, 30)
			memberdist=ClusteredData[ClusteredData["clusternum"] == step]["distance"]
			pmdistiqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["distance"]])
			plt.hist(memberdist,bins, alpha=0.5,color=palette[step])
			plt.gcf().text(0.09, 0.43, "Distance = "+str(int(distcen))+" pc", fontsize=12)
			plt.gcf().text(0.32, 0.43, "Distance IQR = "+str(int(pmdistiqr))+" pc", fontsize=12)
			
#A subplot illustrating the cluster RA and DEC.  To attempt to place the polar coordinates onto a meaningful 2D plane, the RA are 
#transformed relative to the central RA and scaled to the cosine of declination.  Lastly, the median cluster RA and DEC and their 
#IQR in pc space are plotted. 
			plt.subplot(222)
			plt.ylabel("DEC",fontsize=15)
			plt.xlabel("RA Scaled",fontsize=15)
			ramin=min(ClusteredData["ratransform"])
			ramax=max(ClusteredData["ratransform"])
			decmin=min(ClusteredData["dec"])
			decmax=max(ClusteredData["dec"])
			raiqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["ratransform"]])*3.14159/180*distcen
			ramedian=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["ra"]])
			plt.text(ramin, decmax+(decmax-decmin)*0.16, "RA = "+str(round(ramedian,3)), fontsize=12)
			plt.text(ramin+(ramax-ramin)*0.6, decmax+(decmax-decmin)*0.16, "RA IQR = "+str(round(raiqr,2))+" pc", fontsize=12)
			decmedian=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["dec"]])
			deciqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["dec"]])*3.14159/180*distcen
			plt.text(ramin, decmax+(decmax-decmin)*0.08, "DEC = "+str(round(decmedian,3)), fontsize=12)
			plt.text(ramin+(ramax-ramin)*0.6, decmax+(decmax-decmin)*0.08, "DEC IQR = "+str(round(deciqr,2))+" pc", fontsize=12)
			plt.scatter(ClusteredData["ratransform"],ClusteredData["dec"], s=sizes, alpha=1,color=colors)
		
#A subplot illustrating the cluster proper motions.  Additionally, the median cluster RA PM and DEC PM and their IQR in km/s space are plotted. 
			plt.subplot(221)
			plt.ylabel("DEC PM (mas/yr)",fontsize=15)
			plt.xlabel("RA PM (mas/yr)",fontsize=15)
			pmramedian=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmra"]])
			pmdecmedian=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmdec"]])
			plt.gcf().text(0.08, 0.978, "RA PM = "+str(round(pmramedian,2))+" mas/yr", fontsize=12)
			plt.gcf().text(0.08, 0.948, "DEC PM = "+str(round(pmdecmedian,2))+" mas/yr", fontsize=12)
			plt.gcf().text(0.28, 0.978, "RA PM IQR = "+str(round(pmraiqr,2))+" km/s", fontsize=12)
			plt.gcf().text(0.28, 0.948, "DEC PM IQR = "+str(round(pmdeciqr,2))+" km/s", fontsize=12)
			plt.ylim(pmdecmedian-10,pmdecmedian+10)
			plt.xlim(pmramedian-10,pmramedian+10)
			plt.scatter(ClusteredData["pmra"],ClusteredData["pmdec"], s=sizes, alpha=1,color=colors)
			
#A subplot illustrating the cluster CMD.  Additionally, the cluster number and the total count of cluster members are plotted.
			plt.subplot(224)
			plt.xlabel("BP-RP",fontsize=15)
			plt.ylabel("M$_{\mathrm{G}}$",fontsize=15)
			ymin=min(ClusteredData["M_G"])
			ymax=max(ClusteredData["M_G"])
			plt.ylim(ymax,ymin-1)
			plt.xlim(-0.6,3.5)
			number=len(ClusteredData.loc[ClusteredData["clusternum"]==step])
			plt.text(1.95, ymin+0.5, "Cluster Number "+str(step), fontsize=12)
			plt.text(1.95, ymin+1.68, "Cluster Count = "+str(number), fontsize=12)
			plt.scatter(ClusteredData["bp_rp"],ClusteredData["M_G"], s=sizes, alpha=1,color=colors)
			plt.tight_layout()
		
	plt.show()
	