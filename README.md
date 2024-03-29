![License Badge](https://img.shields.io/github/license/saumitrank/disordered-wires?style=plastic)
![Docs](https://img.shields.io/website?down_message=offline&label=docs&up_message=online&url=https%3A%2F%2Fsaumitrank.github.io%2Fdisordered-wires%2F)

# Disordered Wires

[Github pages site with documentation](https://saumitrank.github.io/disordered-wires/)

## Overview

This repo contains code to create various graphene ribbons using [KWANT](https://kwant-project.org/) and study their transport properties. 
There is also an implementation of the Recursive Green's Function (RGF) technique ([see this article](https://doi.org/10.1007/BF01328846)) to calculate the density of states and other transport quantities for a generic quasi-1D system with on-site or hopping-disorder. 


## Usage

The main requirement is the [KWANT](https://kwant-project.org/) package, which can be installed by following the instructions on its site. Apart from this the usual scientific libraries such as numpy, scipy, and matplotlib are required. 

The documentation can be found in the [Github pages site](https://saumitrank.github.io/disordered-wires/) for the repo.


## Results

The density of states and transport data generated using this code has resulted in two publications so far:

### Two-parameter scaling

In the [first paper](https://doi.org/10.48550/arXiv.2112.09748), it was shown that the transport distributions near the topological phase transition of the SSH chain obey a two-parameter scaling paradigm. 
This is in contrast to the usual one-parameter scaling found in the traditional Fokker-Planck approach. 
The dataset for this work can be [found here](https://hdl.handle.net/11299/229873). 
This is the "phase diagram" summarizing the transport regimes:
<p align="center">
  <img alt="Phase Diagram" src="/images/phase-diagram.jpg" width="400" height="300">
</p>

### Multi-criticality in zigzag graphene

In the [second paper](https://doi.org/10.48550/arXiv.2208.05529), the transport distributions near the critical point found previously were used to show the existence of a topological multi-critical point in zigzag graphene. 
Here is the dispersion of zigzag graphene, highlighting the edge states and flat dispersion of the lowest sub-band:
<p align="center">
  <img alt="Phase Diagram" src="/images/zigzag-dispersion.jpg" width="400" height="300">
</p>
And here is a phase diagram illustrating the existence of multiple topological phases as the hopping parameters are tuned, with graphene tuned to the multi-critical point in the middle:
<p align="center">
  <img alt="Phase Diagram" src="/images/phase-diagram-zigzag.jpg" width="400" height="300">
</p>
