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

### 3. [Tutorials]
Tutorials are provided the [**tutorials** subfolder](atmosphere/tutorials/01-densityCurrent).
* [**Tutorial 1**] In this tutorial it is implemented the numerical solution of a non-linear density current. 

### 4. Authors and contributors
**GEA** is currently developed and mantained at [SISSA mathLab](http://mathlab.sissa.it/) by [Dr. Michele Girfoglio](mailto:mgirfogl@sissa.it) in collaboration with [Prof. Annalisa Quaini](mailto:quaini@math.uh.edu) under the supervision of [Prof. Gianluigi Rozza](mailto:gianluigi.rozza@sissa.it)

Contact us by email for further information or questions about **GEA**, or open an ''Issue'' on this website. **GEA** is at an early development stage, so contributions improving either the code or the documentation are welcome, both as patches or merge requests on this website. More to come!

### 5. How to cite
Most of the theoretical aspects behind GEA are deeply explained in the following works. For this reason, if you use this software, please consider to cite them. More to come!

### 6. License
**GEA** is freely available under the GNU LGPL, version 3.
