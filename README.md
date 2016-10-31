# dst-rho-sens 

[![Build Status](https://travis-ci.org/dcervenkov/dst-rho-sens.svg?branch=master)](https://travis-ci.org/dcervenkov/dst-rho-sens)

This is a collection of programs used for a generator-level time-dependent angular CP violation sensitivity study of the B → D* + ρ decay. See my [diploma thesis](http://www-ucjf.troja.mff.cuni.cz/~cervenkov/diploma_thesis/dip_thesis.pdf) for more information.

### DSRhoDataMining
This program reads a ROOT file created by translating a gen/mdst file (panther table) created by EvtGen to a ROOT file (using panther2root) and calculates the transversity basis angular variables θ<sub>t</sub>, θ<sub>b</sub>, φ<sub>t</sub>, Δt and decay type. It then saves these variables in a plain text file (I know).

### DSRhoFit
Does the actual time-dependent angular fit. The following parameters are fitted:
- transversity amplitudes (ap, apa, a0, ata) 
- weak phase (phiw = φ<sub>w</sub> = 2φ<sub>1</sub> + φ<sub>3</sub>) 
- suppressed/favored amplitude ratios (rp, r0, rt)
- strong phases (sp, s0, st)

### DSRhoFitCartesian
Also performs a time-dependent angular fit, but a different set of parameters is used:
- transversity amplitudes (ap, apa, a0, ata) 
- cartesian coordinates (xp, x0, xt, yp, y0, yt, xpb, x0b, xtb, ypb, y0b, ytb)

The problem of the variables r and s is that when relative errors on rs are are large we cannot establish errors on ss. When r = 0, we have no s sensitivity. This is a tricky situation and must be handled in a much more sophisticated manner. Therefore we extract the simple cartesian coordinates that have Gaussian errors and which are to be processed later.

### DSRhoGraphs
Creates various plots (pull, residual, etc.) from the results of DSRhoFitCartesian.

### DSRhoGraphsPolar
Creates various plots (pull, residual, etc.) from the results of DSRhoFit.

### DSRhoRecover
An abandoned naive attempt to extract the physical observables (r,s) from cartesian ones.
