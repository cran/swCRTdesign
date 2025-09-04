\name{Release Notes}
\alias{Release Notes}

\title{Release notes for v4.1}

\description{
Bug fixes: \cr
1) Documentation for swPwr has been corrected. icc, cac, iac are not part of the returned object for swPwr \cr
2) swGlmPwr now includes zeta in the returned object. \cr
3) A bug was fixed in swPlot - formerly, if by.wave=FALSE and combined.plot=FALSE and clusters were not ordered by wave, the clusters were plotted in the wrong wave. This has been fixed. \cr
4) swPwr and swGlmPwr could fail for incomplete designs in which some parameters are not estimable. These functions now drop parameters that are not estimable (with a message) and provide results for the remaining parameters. \cr
5) In the results of swSimPwr for closed cohort designs the label for the individual variance component (zeta) in the object sdcor was missing; it has now been added. \cr
6) if swSummary is called with data with missing cluster or time, it no longer fails. Instead, it deletes the records with missing cluster or time and runs on the remaining data and a warning is issued. This fix is inherited by swPlot since it calls swSummary.

New features: \cr
1) A new function, swSim2, is available. swSim2 extends the functionality of swSim by allowing one to simulate data with an exponentially decaying autocorrelation structure as described in Kasza et al. (2019). In addition, swSim2 uses the glmmTMB function simulate_new to simulate datasets. swSim has been retained for backwards compatability. \cr
2) Argument ar has been added to swPwr, swGlmPwr and swSimPwr, allowing one to specify an exponentially decaying autocorrelation structure as described in Kasza et al. (2019). Note that swSimPwr calls swSim2 (see above) to generate data with an exponentially decaying autocorrelation structure and swSim in all other cases. Thus, for models without an exponentially decaying autocorrelation structure, the behavior of swSimPwr has not changed. \cr
3) Two wrapper functions, swSiz and swGlmSiz, have been added. These functions compute the number of individuals per cluster-period required to achieve a given power for a given design. The functions call swPwr and swGlmPwr, respectively. \cr 
4) In swPwr, if multiple treatment groups and random treatment effects are included in the model, the user can specify different variances for each random treatment effect. \cr
5) swSummary and swPlot now allow arbitrary (ie character) labels for the cluster id. \cr
6) swSummary now returns swDsn.waves - the wave for each cluster - as part of its return object. \cr
7) A new argument - imputeNA - is available in swSummary and swPlot. If imputeNA is to set to TRUE, then swSummary and swPlot attempt to impute the intervention status for cluster-periods with 0 non-missing observations. If FALSE (default, to be consistent with previous versions of swSummary), cluster-periods with 0 non-missing observations are assigned NA in the returned design matricies. See details in the swSummary and swPlot documentation.\cr
8) Cluster names are now included on plots in swPlot when appropriate. \cr
9) The legend giving waves and clusters printed on the plot when by.wave=FALSE and combined.plot=TRUE has been removed.
}

