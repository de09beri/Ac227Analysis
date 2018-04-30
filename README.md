# Ac227Analysis

**(1)** Set all environment variables in your personal .bashrc. *See Examples/ac227examplebashrc.txt*  
	The ones that will be important are   
	$P2X_PHYSDATA, the calibrated PhysPulse data that you want to analyze  
	$P2X_ANALYZED, where the Ac227TreePlugin results live

	$AD_AC227ANALYSIS_SCRIPTS, the /AnalysisScripts directory in this repository   

	$AD_AC227ANALYSIS_DATA_PLOTS, where you want the analysis plots to live

**(2)** To run the Ac227 Plugin:

	Make sure you have the most updated version of P2x

	Edit the config file to run the **Ac227TreePlugin** with the desired PSD, energy, and position cuts. *See Examples/ADAnalyzer.cfg*

	Move to $P2X_ANALYSIS_CODE/ControlScripts

	Run `./LaunchBatchScript.py --anacal <Data type> --mode <Config file> --jmax <# of cores>`

		where <Data type> = WetCommissioning, 180316_Rampdown, or 180316_Background	

	The results will be placed in $P2X_ANALYZED 

**(3)** To combine the root trees created by the plugin:

	Move to Ac227Analysis

	Run `./CombineTrees.sh`

	This will combine all trees in a series into one tree per series 

	Resulting tree will be placed where the series directory lives

**(4)** To create plots:

	Move to Ac227Analysis

	Edit RunReadAc227Trees.sh 

		numTrees = 590/n, where n is the number of time periods that you would like to look at. n = 1, will look at all data 

		PLOTFLAG = <1,2>, will determine which plots are created

			     = 1 will create plots of PSD vs. Energy, <Dt, Energy, PSD, Dz> distributions, <Rate, Efficiency, Energy Mean, PSD Mean, Dz Mean,...> per cell. 

				   This is suggested for using when n = 1

				 = 2 will create plots of <Rate, Efficiency, Energy Mean, PSD Mean, Dz Mean, Position Sigma> over time 

				   This is suggested for using when n > 1

	Run `./RunReadAc227Trees.sh`

	Resulting plots will be placed in $AD_AC227ANALYSIS_DATA_PLOTS	

