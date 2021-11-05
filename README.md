# design2phase package

## Power estimation via simulation for two-phase sampling designs
 
This package provides power and efficiency estimation (via simulation) in assessing the influence of a novel biomarker on time-to-event outcome using two-phase sampling approach.

## Installation

This is the development version of the package. It can be installed from github with

```r
devtools::install_github("Fgraziano/design2phase")
library("design2phase")
```

## Basic usage

The core function is called ```PowerIIphase()``` and includes the following arguments:
```r
PowerIIphase(pBM,betaBM,pstrata=NULL,betastrata=NULL,acc.aux=NULL,design2p=NULL,N,n,cens=0.1, tau=2,lambda=0.1,k=0.9,B=1000,seed=NULL)
```

All possible sampling designs are the following: 

- Simple Random Sampling (SRS), 
- Case-Control (CC), 
- Probability Proportional to Size (PPS), 
- Nested CaseControl (NCC),
- Countermatching (CM).

Default sampling designs includes SRS and CC designs. 
The main function has few mandatory parameters (```pBM, betaBM, N, n```); other parameters are set with default, but can be modified according to the setting of interest and the information available from the phase I. 

To perform the simplest scenario:

```r
perfBM <- PowerIIphase(betaBM=0.91,  pBM=0.25,  N=400, n=c(80,100,120) , seed=467)
```

The output in ```PowerIIphase object``` is:

```r
$PhaseI_events
[1] 82

$designs
[1] "SRS"            "CC strat event"

$n
[1]  80 100 120

$PhaseII_Performance
$PhaseII_Performance$`n=80`
               Blength sample nevent Power   deff
SRS               1000     80     16 0.473 1.0000
CC strat event    1000     80     40 0.608 1.5287

$PhaseII_Performance$`n=100`
               Blength sample nevent Power   deff
SRS               1000    100     21 0.566 1.0000
CC strat event    1000    100     50 0.702 1.4349

$PhaseII_Performance$`n=120`
               Blength sample nevent Power   deff
SRS               1000    120     25 0.599 1.0000
CC strat event    1000    120     60 0.761 1.4283


attr(,"class")
[1] "PowerIIphase"
```
The output includes:

- Expected number of events for phase I ```$PhaseI_events```
- Explored sampling designs ```$designs```
- Planned sample sizes of phase II (n) ```$n```
- Performance over B simulations ```$PhaseII_Performance```

It is possible to plot the power curves if the hypothetical sample sizes of the subsample n are at least 2 

```r
PlotPower(perfBM)
```
## Stratification Process

### stratify for the risk factor/confonder
To consider the stratification for the risk factor, ```betastrata``` (expected beta coefficient, ln(HR)) and ```pstrata``` (the prevalence of the stratum variable) arguments need to be specify in ```PowerIIphase``` function. 

### stratify for the auxiliary variable
To consider the stratification for auxiliary variaible, ```acc.aux``` ( expected sensibility and specificity (accuracy)) argument need to be specify in ```PowerIIphase``` function. 

## For more details

```r
?PowerIIphase
?PlotPower
```

## References

Graziano, F., Valsecchi, M.G. & Rebora, P. Sampling strategies to evaluate the prognostic value of a new biomarker on a time-to-event end-point.BMC Med Res Methodol 21, 93 (2021). https://doi.org/10.1186/s12874-021-01283-0

## MIT License

Copyright (c) 2020 design2phase

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
