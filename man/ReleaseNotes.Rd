\name{Release Notes}
\alias{Release Notes}

\title{Release notes for v4.0}

\description{
Changes from prior versions:

1) A bug in the way that swGlmPwr handled fractional treatment effects has been fixed.

2) In SwDsn, the extra.time argument has been deprecated and replaced by extra.trt.time. If you enter extra.time it will be used as extra.trt.time.

3) In swPwr, the warning "The order of variance component arguments has changed for swPwr (in version 2.2.0, it was tau, eta, rho, sigma); please modify existing code if necessary." has been removed.

4)	In swPwr, the closed form power calculation (only applicable to a special subset of situations) is no longer provided.

5)  In swPwr and swSim variance components now default to 0 instead of NULL so you only have to enter non-zero variance components. However, you cannot mix the two types of specifications (ie either use the tau, eta, gamma, zeta form or the icc, cac, iac form)

6)	In swSim, the retTimeOnTx argument has been deprecated and timeOnTx.var is always returned

New features:
1)	In swDsn, you can now directly specify an arbitrary stepped wedge design using the swBlk argument. This also allows one to specify multiple treatment levels

2)	In swPwr and swGlmPwr, closed cohort designs may be specified by entering a non-zero value for either zeta or iac. If zeta/iac is 0 (the default) then a cross-sectional sampling design is assumed.

3)	In swPwr and swGlmPwr, ETI-type estimators (Kenny et al) may be specified using the H argument. If the elements of H do not sum to 1 they are automatically renormalized.

4)	swSim can simulate cross-sectional or closed cohort data, data with multiple treatment levels and data with an exposure-time-varying treatment effect.

5)	swPlot allows specification of x and y axis limits.

6)  A new function, swSimPwr, is included to compute power by simulation. Analyses are based on linear mixed models (lmer) or generalized linear mixed models (glmer).
}

