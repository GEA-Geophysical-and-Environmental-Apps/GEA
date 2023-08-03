<p align="center">
    <img alt="GEA" src="docs/gea_logo.png" width="200" />
</p>

## GEA - Geophysical and Environmental Applications ##

<p align="center">
    <a href="https://www.gnu.org/licenses/lgpl-3.0" target="_blank">
        <img alt="Software License" src="https://img.shields.io/badge/License-LGPL%20v3-blue.svg">
    </a>
</p>

### 0. Introduction
**GEA** is an implementation in [**OpenFOAM**](https://www.openfoam.com) of several atmosphere and ocean models.

### 1. Prerequisites
**GEA** requires
* [**OpenFOAM**](https://www.openfoam.com)
* [**swak4Foam**](https://openfoamwiki.net/index.php/Contrib/swak4Foam)

### 2. Installation and usage
First of all you need to source the bashrc file of your installation of **OpenFOAM**. This is of course depending on the location of your OpenFOAM installation and of your particular version of OpenFOAM. Then navigate to the folder where you want to install GEA. Now you can clone the **GEA** repository inside the selected folder
```
git clone https://github.com/GEA-Geophysical-and-Environmental-Apps/GEA
```
and you can compile the solvers of your interest by navigating inside the src folder and using wmake.

**GEA** has been tested on Centos Stream 8 and OpenFOAM v2106 but it can be compiled on any Linux distribution with a compiled version of OpenFOAM. 

### 3. Tutorials
Tutorials are provided in the [**tutorials** subfolder](atmosphere/tutorials).
* [**Tutorial 1**](atmosphere/tutorials/01-densityCurrent): In this tutorial it is implemented the numerical solution of a non-linear density current. 

### 4. Authors and contributors
**GEA** is currently developed and mantained at [SISSA mathLab](http://mathlab.sissa.it/) by [Dr. Michele Girfoglio](mailto:mgirfogl@sissa.it) in collaboration with [Prof. Annalisa Quaini](mailto:quaini@math.uh.edu) under the supervision of [Prof. Gianluigi Rozza](mailto:gianluigi.rozza@sissa.it)

Contact us by email for further information or questions about **GEA**, or open an ''Issue'' on this website. **GEA** is at an early development stage, so contributions improving either the code or the documentation are welcome, both as patches or merge requests on this website. More to come!

### 5. How to cite
Most of the theoretical aspects behind GEA are deeply explained in the following works. For this reason, if you use this software, please consider to cite them.

* Girfoglio, Quaini, Rozza. *A POD-Galerkin reduced order model for the Navier–Stokes equations in stream function-vorticity formulation*. Computers & Fluids, vol. 244, p. 105536, 2022. [[DOI](https://doi.org/10.1016/j.compfluid.2022.105536)] [[arXiv](https://arxiv.org/abs/2201.00756)].

* Girfoglio, Quaini, Rozza. *A novel Large Eddy Simulation model for the Quasi-Geostrophic Equations in a Finite Volume setting*. Journal of Computational and Applied Mathematics, vol. 418, p. 114656, 2023. [[DOI](https://doi.org/10.1016/j.cam.2022.114656)] [[arXiv](https://arxiv.org/abs/2202.00295)].

* Girfoglio, Quaini, Rozza. *A linear filter regularization for POD-based reduced-order models of the quasi-geostrophic equations*. Comptes Rendus Mècanique, p. 1-21, 2023. [[DOI](https://doi.org/10.5802/crmeca.183)] [[arXiv](https://arxiv.org/abs/2211.16851)].

* Girfoglio, Quaini, Rozza. *Validation of an OpenFOAM-based solver for the Euler equations with benchmarks for mesoscale atmospheric modeling*. AIP Advances, vol. 13, p. 055024, 2023. [[DOI](https://doi.org/10.1063/5.0147457)] [[arXiv](https://arxiv.org/abs/2302.04836)].

### 6. License
**GEA** is freely available under the GNU LGPL, version 3.
