<a name="top"></a>

[![Join the chat at ncimm.slack.com](https://img.shields.io/badge/slack-Join%20Chat-blue.svg)](https://ncimm.slack.com)

[![License](https://img.shields.io/badge/license-Public Domain-blue.svg)]()

### SWMMengine, pure Fortran engine for the U.S. EPA Stormwater Management Model


This project is building a new parallel time-marching solver for unsteady flow solutions in the U.S. EPA Stormwater Management Model (SWMM). The prior link-node hydraulic solver has been demonstrated to have problems in maintaining mass conservation and parallelization typically is saturated with 10 or 15 processors. The new computational engine is being designed to work with the existing user interface and input files of SWMM (version 5.1). Key advances in the new engine are: 

+ Finite-volume sub-discretization of links;
+ Representation of nodes (junctions) as physical elements with multiple faces and different possible heads on each face; 
+ Use of artificial compressibility for surcharged pipe flow; 
+ Explicit time-marching of Saint-Venant equation for simplified parallel solution,
+ Parallelization through coarray Fortran (i.e., Fortran 2008 or later) that provides run-time parallelization on any number of processors;
+ SWMMengine is Fortran 2008+ standard compliant;
+ SWMMengine supports _ascii_, _binary_ (Base64 encoding) and _raw_ file formats;
+ SWMMengine is a _Free_, _Open Source_ Project under _Public Domain_ license.

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-pass%20(v6.0.1+)-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/Intel-pass%20(v16.x+)-brightgreen.svg)]()

---

[Main features](#main-features) | [Copyrights](#copyrights) | [Documentation](#documentation) | [Acknowledment](#acknowledment) |

---

## Main features

### Parallel Support

SWMMengine can be safely used in parallel *environments*, handling multiple procesors. We are currently working on making it parallel using the Fortran Coarrays and make it available as a next version to the community.

## Copyrights

SWMMengine is an open source project, it is distributed under a public domain:

+ [Public Domain](https://opensource.org/node/878).

Anyone is interest to use, to develop or to contribute to SWMMengine is welcome, feel free to select the license that best matches your soul!

Go to [Top](#top)

## Documentation

Detailed documentation of SWMMengine is contained into this README file. Detailed documentation of SWMM can be found at [U.S. EPA](https://www.epa.gov/water-research/storm-water-management-model-swmm) website.

## Acknowledment

This project is being  developed by the National Center for Infrastructure Modeling and Management (NCIMM), under funding from the U.S. Environmental Protection Agency with Cooperative Agreement No. 83595001 awarded  to the University of Texas at Austin.

Go to [Top](#top)
