\name{Release Notes}
\alias{Release Notes}

\title{Release notes for v4.2}

\description{
This new version has been released because I received an email indicating 
that glmmTMB might be removed from CRAN.

This affects the functions swSimPwr (which calls glmmTMB) and swSim2 (which calls 
simulate_new in the glmmTMB package) in the swCRTdesign package. glmmTMB was used
because of bugs in the earlier lme4 implementation of the ar() function. These
bugs have been fixed and I have now updated the swSimPwr function to call lmer/glmer 
instead of glMMTMB and swSim2 to call simulate instead of simulate_new. 

If you are using autoregressive models (for power or simulation) it is important 
that you have the latest version of the lme4 library, which can be downloaded 
from GitHub:

remotes::install_github("lme4/lme4")
}

