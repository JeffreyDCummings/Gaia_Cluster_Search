def clusterplot(ClusteredData,MembershipPlot):
	import seaborn as sns
	import matplotlib.pyplot as plt
	import scipy.stats as spy
	import numpy as np
	import pandas as pd
	import matplotlib.colors as pltc

#We count the number of found clusters, set up the color palette, and set tick label font size to 12 for clarity.

	uniquegroup=ClusteredData.clusternum.unique()
	palette = sns.color_palette('bright', max(uniquegroup)+1)
	plt.rcParams['xtick.labelsize'] = 12
	plt.rcParams['ytick.labelsize'] = 12

#If plotting a membership probability color map is wanted, the input dataframe must be separated into a different dataframe for
#each cluster.  Doing this with a dictionary is fast and efficient.  A second dictionary is made to set the color map for each cluster
#number.  If plotting membership probability is not needed, we set a unique color to each cluster and set nonclustered field stars to black. 
	if MembershipPlot==1:
		ClusteredDict={elem: pd.DataFrame for elem in uniquegroup}
		for key in ClusteredDict.keys():
			ClusteredDict[key] = ClusteredData[:][ClusteredData.clusternum == key]
		cmapdict={-1:"Greys",0:"Blues",1:"Oranges",2:"Greens",3:"Reds",4:"Purples",5:"PuBu",6:"YlGn",7:"Greys",8:"Purples",9:"Blues",10:"Greens",11:"Reds",12:"Oranges",13:"PuBu",14:"YlGn",15:"Greys"}
		normalize = pltc.Normalize(vmin=0, vmax=1)
	else:
		colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in ClusteredData["clusternum"]]

#This loop iterates through all clustered groups found with hdbscan and plots the requested clusters that pass the checks.
	for step in range(max(uniquegroup)+1):

#The median stellar distance is found (in parsecs), which is used to transform proper motions and RA and DEC into spatial scales.
#Because these distributions are not always Gaussian, the interquartile ranges are used to analyze parameter variations throughout a cluster.
		distcen=np.median(ClusteredData.loc[ClusteredData["clusternum"]==step,["distance"]])
		pmraiqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmra"]])*3.14159/180*distcen*0.271795
		pmdeciqr=spy.iqr(ClusteredData.loc[ClusteredData["clusternum"]==step,["pmdec"]])*3.14159/180*distcen*0.271795

#Hdbscan will typically create a cluster from the overdensity in the field disk stars, producing apparent clusters but with large PM variation.
#Therefore, we require that the transformed PM variation be less than 3 km/s in each axis to be displayed as a valid cluster.  If a legitimate 
#cluster is identified but retains a lot of noise/non-members, a 3 km/s cut may remove it from display.  In these cases, you may increase the 
#cut to 4 or 5 km/s, etc., to first analyze its output, but in the end I recommend adjusting the hdbscan min_samples or the parameter scalings.
		if (pmraiqr < 3 and pmdeciqr < 3):

#To set up the plots, we set a separate Figure for each cluster.  Forr when membership probability is not shown, we set all clustered data to 
#a large size and nonclustered data to small data size.  Lastly, set tick label font size to 12 for clarity.
			plt.figure(step,figsize=(14, 10))
			if MembershipPlot != 1:
				sizes = [4 if x == step else 0.003 for x in ClusteredData["clusternum"]]

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
			
#A subplot illustrating the cluster scaled RA and DEC.  The median cluster RA and DEC and their IQR in pc space are plotted. 
#Whether or not the membership probability colormap was desired is checked, and a separate color map for each cluster is set 
#or the simple setup is adopted.
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
			if MembershipPlot==1:
				for g in uniquegroup:
					if g==-1:
						plt.scatter(ClusteredDict[g]["ratransform"],ClusteredDict[g]["dec"],s=0.005, alpha=1,color="black")
					elif g==step:
						plt.scatter(ClusteredDict[g]["ratransform"],ClusteredDict[g]["dec"],s=6, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)
						cbar=plt.colorbar()
						cbar.set_label('Membership Probability')
					else:
						plt.scatter(ClusteredDict[g]["ratransform"],ClusteredDict[g]["dec"], s=0.02, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)				
			else:
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
			if MembershipPlot==1:
				for g in uniquegroup:
					if g==-1:
						plt.scatter(ClusteredDict[g]["pmra"],ClusteredDict[g]["pmdec"],s=0.005, alpha=1,color="black")
					elif g==step:
						plt.scatter(ClusteredDict[g]["pmra"],ClusteredDict[g]["pmdec"],s=5, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)
						cbar=plt.colorbar()
						cbar.set_label('Membership Probability')
					else:
						plt.scatter(ClusteredDict[g]["pmra"],ClusteredDict[g]["pmdec"], s=0.02, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)				
			else:
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
			if MembershipPlot==1:
				xlabelpos=1.57
				for g in uniquegroup:
					if g==-1:
						plt.scatter(ClusteredDict[g]["bp_rp"],ClusteredDict[g]["M_G"],s=0.005, alpha=1,color="black")
					elif g==step:
						plt.scatter(ClusteredDict[g]["bp_rp"],ClusteredDict[g]["M_G"],s=5, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)
						cbar=plt.colorbar()
						cbar.set_label('Membership Probability')
					else:
						plt.scatter(ClusteredDict[g]["bp_rp"],ClusteredDict[g]["M_G"], s=0.02, alpha=1,c=ClusteredDict[g]["clusterprob"],cmap=cmapdict[g],norm=normalize)				
			else:
				xlabelpos=1.95
				plt.scatter(ClusteredData["bp_rp"],ClusteredData["M_G"], s=sizes, alpha=1,color=colors)
			plt.text(xlabelpos, ymin+0.5, "Cluster Number "+str(step), fontsize=12)
			plt.text(xlabelpos, ymin+1.68, "Cluster Count = "+str(number), fontsize=12)
			plt.tight_layout()
		
	plt.show()
	