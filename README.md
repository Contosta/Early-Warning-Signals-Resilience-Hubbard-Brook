# Early-Warning-Signals-Resilience-Hubbard-Brook
Data and code for analyzing changing resilience at the Hubbard Brook Experimental Forest using early warning signals. These materials form the basis for reproducing the results presented in a manuscript currently in press in Environmental Research Letters.

There are five files containing code for performing the analysis. While source files exist to run each script ad hoc, they are shown below in the order in which the analysis was conducted. 

**HB-Data_Preparation.R** does initial data processing. **Data inputs** are from the _Original Data_ folder. **Data outputs** are preprocessed data that have been reformatted so that all data sources are aggregated at annual timesteps. They are located in the _Preprocessed Data_ folder. Original source files can also be obtained from the Hubbard Brook Data Catalog at https://hubbardbrook.org/data-catalog/. Dataset citations with associated DOIs are:

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021a). Continuous precipitation and stream chemistry data, Hubbard Brook Ecosystem Study, 1963 – present. ver 6. Environmental Data Initiative. https://doi.org/10.6073/pasta/ee9815b41b79c134fd714736ce98676a (Accessed 2021-08-04).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021b). Hubbard Brook Experimental Forest: Chemistry of Streamwater – Monthly Fluxes, Watershed 1, 1963 - present ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/a912e690ebe20d6a10352143c0c8a24d (Accessed 2021-08-12).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021c). Hubbard Brook Experimental Forest: Chemistry of Streamwater – Monthly Fluxes, Watershed 2, 1963 - present ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/1a22e1a37f6802c632426f137d7f9049 (Accessed 2021-08-12).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021d). Hubbard Brook Experimental Forest: Chemistry of Streamwater – Monthly Fluxes, Watershed 4, 1963 - present ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/59ec30d409d3ef05e5f76f48eb18d325 (Accessed 2021-08-05).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021e). Hubbard Brook Experimental Forest: Chemistry of Streamwater – Monthly Fluxes, Watershed 5, 1963 - present ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/a471b540dd141e361b137bad8fc92389 (Accessed 2021-08-05).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021f). Hubbard Brook Experimental Forest: Chemistry of Streamwater – Monthly Fluxes, Watershed 6, 1963 - present ver 17. Environmental Data Initiative. 

Groffman, P.M. and Martel, L. (2020). Long-term measurements of microbial biomass and activity at the Hubbard Brook Experimental Forest, 1994 - present ver 21. Environmental Data Initiative. https://doi.org/10.6073/pasta/3f608226a1ed499e8fa3cd188e70757c (Accessed 2021-08-12).

Halman, J.M., Schaberg, P.G., Hawley, G.J., Hansen, C.F. and Fahey, T.J. (2008). Sugar maple and American beech tree cores from NuPert plots at the Hubbard Brook Experimental Forest, NH. https://www.uvm.edu/femc/dendro/project/17 (Accessed  2021-08-12).

Ayres, M.P. and Holmes, R.T. (2020). Long-term trends in abundance of Lepidoptera larvae at Hubbard Brook Experimental Forest and three additional northern hardwood forest sites, 1986-2018 ver 8. Environmental Data Initiative.  https://doi.org/10.6073/pasta/5d2a8c67c5a3278032b2b14d66c09a7f (Accessed 2020-10-09).

Holmes, R.T. (2019). Bird Abundances at the Hubbard Brook Experimental Forest (1969-present) and on three replicate plots (1986-2000) in the White Mountain National Forest ver 7. Environmental Data Initiative. https://doi.org/10.6073/pasta/fa7380ace83abcb62c35872522d118d4 (Accessed 2020-10-09).

USDA Forest Service, Northern Research Station. (2020). Hubbard Brook Experimental Forest: Daily Temperature Record, 1955 – present ver 9. Environmental Data Initiative. https://doi.org/10.6073/pasta/e7c793b98b895de2bb5e505f9ff5e0cb (Accessed 2021-08-12).

USDA Forest Service, Northern Research Station. (2021). Hubbard Brook Experimental Forest: Daily Precipitation Rain Gage Measurements, 1956 - present ver 16. Environmental Data Initiative. https://doi.org/10.6073/pasta/1391b5fa8c342a30df1957a78b38096e (Accessed 2021-08-12).

Hubbard Brook Watershed Ecosystem Record (HBWatER). (2021f). Hubbard Brook Experimental Forest: Chemistry of Precipitation – Monthly Fluxes, Watershed 6, 1964 - present ver 11. Environmental Data Initiative. https://doi.org/10.6073/pasta/39887003e9c00b21953f9f5f03b558e7 (Accessed 2021-08-12).

**HB-Sensitivity.R** performs sensitivity analysis on EWS calculated from Hubbard Brook data, examining Kendall's tau and p-value statistics for different window lengths and smoothing bandwiths using the _earlywarnings_ R package. **Data inputs** are from the _Preprocessed Data_ folder. **Outputs** include diagnostic plots and related tables of Kendall's Tau and p-values that result from using different smoothing bandwiths for detrending data and different window lengths for calculating early warning signals. **Outputs** also include _Figure_S3_, which shows the distribution of p-values for change over time for each of four early warning signals (standard deviation, autocorrelation, skewness, and kurtosis) over five possible bandwiths for smoothing and detrending the data and over five possible sliding windows for calculating early warning signals. **Data outputs** are in the _Sensitivity Analysis Figures_ and _Sensitivity Analysis Tables_ folders for the results of the sensivitiy analyses and in the _EWS Summary Figures_ for _Figure_S3_.  

**HB-EWS.R** calculates generic early warning signals using the earlywarnings R package. **Data inputs** are preprocessed data located in the _Preprocessed Data_ folder. **Data outputs** are figures and associated tables of generic early warning signals. **Outputs** are located in the _Early Warning Signals Figures_ and _Early Warning Tables_ folders. 

**HB-Trends.R** calculates trends over time for early warning signals (standard deviation, autocorrelation, skewness, and kurtosis) calculated from detrended data. **Input data** are calculated early warning signals located in the _Early Warning Signals_ Tables folder. **Outputs** are results from Mann-Kendall analysis of change over time and are written as two tables: _alltrends.csv_ is in long format and is an input for further analysis in the **HB-Synthesis script**; _tabs2.csv_ is in wide format and is the basis for _Table S2_ in the manuscript. Both output tables are written to the __EWS Summary Tables_ folder.

**HB-Synthesis.R** contains code for producing most of the figures and tables presented in the manuscript. **Inputs** include hypothesized drivers of climate, precipitation chemistry, and three climate oscillation indices located in the _Original Data_ folder. **Inputs** also consist of calculated early warning signals for each response variable from the _Early Warning Signals Tables_ folder, as well as trends over time in early warning signals (_alltrends.csv_) in the _EWS Summary Tables_ folder. **Outputs** include __Figure_2, Figure_3, Figure_4_, and _Figure_S4_ in the manuscript, which are written to the _EWS Summary Figures_ folder. Outputs also consist of Kendall's Tau and associated p-values that compare early warning signals to hypothesized drivers of climate and precipitation chemistry (_EWS_regtab.csv_) that form the basis of _Tables S3-S6,_ as well as p- and r2 values derived by fitting hypothesized drivers onto nonmetric multidimensional scaling (NMDS) ordinations of early warning signals (_allenvfits.csv_). Both output tables are located in the _EWS Summary Tables_ folder.

**EWS Example.R** produces a figure (_Figure_S5_) showing original NO3 flux data, detrended NO3 flux data, and  distributions of detrended NO3 flux data for the biogeochemical reference watershed (Watershed 6) and the whole-tree harvested watershed (Watershed 5) at Hubbard Brook as an example of how early warning signals such as standard deviation, skewness, and kurtosis can change over time in different directions. **Inputs** include original and detrended NO3 flux data in the _Detrended Data_ folder.  **Outputs** include _Figure_S5_ in the manuscript, which is written to the _EWS Summary Figures_ folder. 
