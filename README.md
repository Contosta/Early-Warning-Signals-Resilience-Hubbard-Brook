# Early-Warning-Signals-Resilience-Hubbard-Brook
Data and code for analyzing changing resilience at the Hubbard Brook Experimental Forest using early warning signals. These materials form the basis for reproducing the results presented in a manuscript currently under review.

There are five files containing code for performing the analysis. While source files exist to run each script ad hoc, they are shown below in the order in which the analysis was conducted. 

**HB-Data_Preparation.R** does initial data processing. **Data inputs** are from the _Original Data_ folder. **Data outputs** are preprocessed data that have been reformatted so that all data sources are aggregated at annual timesteps. They are located in the _Preprocessed Data_ folder. 

**HB-Sensitivity.R** performs sensitivity analysis on EWS calculated from Hubbard Brook data, examining Kendall's tau and p-value statistics for different window lengths and smoothing bandwiths using the _earlywarnings_ R package. **Data inputs** are from the _Preprocessed Data_ folder. **Outputs** include diagnostic plots and related tables of Kendall's Tau and p-values that result from using different smoothing bandwiths for detrending data and different window lengths for calculating early warning signals. **Outputs** also include _Figure_S3_, which shows the distribution of p-values for change over time for each of four early warning signals (standard deviation, autocorrelation, skewness, and kurtosis) over five possible bandwiths for smoothing and detrending the data and over five possible sliding windows for calculating early warning signals. **Data outputs** are in the _Sensitivity Analysis Figures_ and _Sensitivity Analysis Tables_ folders for the results of the sensivitiy analyses and in the _EWS Summary Figures_ for _Figure_S3_.  

**HB-EWS.R** calculates generic early warning signals using the earlywarnings R package. **Data inputs** are preprocessed data located in the _Preprocessed Data_ folder. **Data outputs** are figures and associated tables of generic early warning signals. **Outputs** are located in the _Early Warning Signals_ and _Early Warning Tables_
folders. 

**HB-Trends.R** calculates trends over time for early warning signals (standard deviation, autocorrelation, skewness, and kurtosis) calculated from detrended data. **Input data** are calculated early warning signals located in the _Early Warning Signals_ Tables folder. **Outputs** are results from Mann-Kendall analysis of change over time and are written as two tables: _alltrends.csv_ is in long format and is an input for further analysis in the the **HB-Synthesis script**; _tabs2.csv_ is in wide format and is the basis for _Table S2_ in the manuscript. Both output tables are written to the _EWS Summary Tables folder._

**HB-Synthesis.R** contains code for producing most of the figures and tables presented in the manuscript. **Inputs** include hypothesized drivers of climate, precipitation chemistry, and three climate oscillation indices located in the _Original Data_ folder. **Inputs** also consist of calculated early warning signals for each response variable from the _Early Warning Signals Tables_ folder, as well as trends over time in early warning signals (_alltrends.csv_) in the _EWS Summary Tables_ folder. **Outputs** include _Figure_2, Figure_3, Figure_4, and _Figure_S4_ in the manuscript, which are written to the _EWS Summary Figures_ folder. Outputs also consist of Kendall's Tau and associated p-values that compare early warning signals to hypothesized drivers of climate and precipitation chemistry (_EWS_regtab.csv_) that form the basis of _Tables S3-S6,_ as well as p- and r2 values derived by fitting hypothesized drivers onto nonmetric multidimensional scaling (NMDS) ordinations of early warning signals (_allenvfits.csv_). Both output tables are located in the _EWS Summary Tables_ folder.


