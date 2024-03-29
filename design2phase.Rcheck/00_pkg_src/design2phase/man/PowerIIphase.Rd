% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PowerIIphase.R
\name{PowerIIphase}
\alias{PowerIIphase}
\title{Power estimation via simulation for two-phase sampling designs}
\usage{
PowerIIphase(
  pBM,
  betaBM,
  pstrata = NULL,
  betastrata = NULL,
  acc.aux = NULL,
  design2p = NULL,
  N,
  n,
  cens = 0.1,
  tau = 2,
  lambda = 0.1,
  k = 0.9,
  B = 1000,
  seed = NULL
)
}
\arguments{
\item{pBM}{expected prevalence of the biomarker of interest in the cohort (phase I).}

\item{betaBM}{expected beta coefficient (ln(HR))  of the new (bio)marker of interest on the time-to-event endpoint.}

\item{pstrata}{prevalence of the stratum variable to be used in the design. Default is NULL for not stratified sampling.}

\item{betastrata}{expected beta coefficient (ln(HR))   of the stratum to be used in the design on the time-to-event endpoint. Default is NULL for not stratified sampling.}

\item{acc.aux}{vector of the expected sensibility and specificity (accuracy) of the auxiliary with respect to the biomarker (default NULL for no auxiliary variable used in the design). See the details for the definition.}

\item{design2p}{a list of extra sampling designs performed during the sampling process, with SRS and Case-Control (CC) fitted by default. Optional choice are "PPS", "NCC" and "CM". "CM" is possible only if auxiliary variable is available. If stratum or auxiliary variable are available, CC and PPS stratified by these variables are automatically performed. See the details for the definition.}

\item{N}{a number providing the sample size of the full cohort (phase I).}

\item{n}{a vector providing the sample size of the subsample (phase II).}

\item{cens}{rate of censoring from exponential distribution. Parameter useful for simulating censoring in the cohort (phase I). Default is a rate=0.1.}

\item{tau}{maximum follow-up time in the cohort. Parameter useful for simulating time-to-event in the cohort (phase I). Default is 2.}

\item{lambda}{scale parameter of Weibull baseline hazard of event. Parameter used for simulating time-to-event in the cohort (phase I). Default is 0.1.}

\item{k}{shape parameter of Weibull baseline hazard of event. Parameter used for simulating time-to-event in the cohort (phase I). Default is 0.9.}

\item{B}{number of datasets generated by the simulation. Default is 1000.}

\item{seed}{of the random number generation. Default is NULL.}
}
\value{
PowerIIphase() returns an object of class PowerIIphase containing:
\itemize{
 \item{\code{phaseI_events} expected number of events for phase I;}
 \item{\code{designs} sampling designs;}
 \item{\code{n} planned sample sizes of phase II (n);}
 \item{\code{PhaseII_performance} performance over B simulations. The data frame includes: \itemize{
  \item number of simulations (\code{Blength})
  \item actual sample size of phase II (\code{sample})
  \item number of events included in phase II (\code{nevent})
  \item performances of the sampling designs (\code{Power} and \code{deff})
}}}
}
\description{
Power and efficiency estimation (via simulation) in assessing the influence of a novel biomarker on time-to-event outcome using two-phase sampling approach.
The code is useful to plan a subsampling (called phase II) of a cohort/register (called phase I), on which to measure a new marker. Sampling design could be stratified for a variable available in the cohort (phase I e.g. by event or by possible counfounder/risk factor). If a stratification by an auxiliary variable of the new (bio)marker is considered, expected sensibility and specificity (as measures of accuracy of the auxiliary variable with respect to the new (bio)marker) in \code{acc.aux} argument need to be included.
Marker and stratum are assumed to be dichotomous variables.
The mandatory arguments are \code{betaBM}, \code{pBM}, \code{N} and \code{n}.
}
\details{
The time-to-event is generated from a Weibull hazard model with the scale parameter lambda and shape parameter k. Censoring is generated from an exponential function. Expected number of events in phase I is stored in the \code{PowerIIphase} object as \code{phaseI_events}.

Power is estimated as the proportion of simulated datasets for which the null hypothesis of no effect of the biomarker on the time-to-event  is rejected. The test is the Wald test from a two-phase Cox model weighted by the inverse of the empirical inclusion probability to account for the specific sampling design.

Design effect (\code{deff}) is defined as the ratio between the estimated variance of biomarker in each sampling design and the one in Simple Random Sampling (SRS).

All possible designs are the following: Simple Random Sampling (SRS), Case-Control (CC), Probability Proportional to Size (PPS), Nested Case Control (NCC),Countermatching (CM).

   According with the parameters provided for the arguments (\code{pstrata}, \code{betastrata}, \code{acc.aux}), the following designs are fitted: \itemize{
 \item{\code{pstrata=NULL}, \code{betastrata=NULL}, \code{acc.aux=NULL} : SRS and Case-Control (CC) stratified by event. (PPS and NCC in \code{design2p} argument are optional)}
 \item{if stratum is available: SRS, CC stratified by event, and CC stratified by event and stratum. (PPS, NCC in \code{design2p} argument are optional)}
 \item{if auxiliary is available: SRS, CC stratified by event and CC stratified by auxiliary (CM in \code{design2p} argument is optional)}
 }


The accuracy of the auxiliary with respect to the biomarker is defined by expected sensibility and specificity. Sensibility is the expected probability to have a 'positive' auxiliary variable in subjects with a 'positive' biomarker. Specificity is the expected probability to have a 'negative' auxiliary variable in subjects with a 'negative' biomarker.

An equal number of cases and controls are sampled (balanced 1:1). Sample size of phase II in NCC and CM designs might results slightly lower than the fixed one due to the lack of subjects event-free at the time of event of the cases.
Graziano, F., Valsecchi, M.G. & Rebora, P. Sampling strategies to evaluate the prognostic value of a new biomarker on a time-to-event end-point. BMC Med Res Methodol 21, 93 (2021).\url{https://doi.org/10.1186/s12874-021-01283-0}
}
\examples{
\dontrun{
library(design2phase)
# no stratum and no auxiliary
perfBM <- PowerIIphase(betaBM=0.91,  pBM=0.25,  N=400, n=c(80,100,120) , seed=467)
PlotPower(perfBM)

}
}
\seealso{
use \url{PlotPower} function to provide plot power curve estimates
}
