# k-Wave: A MATLAB toolbox for the time-domain simulation of acoustic wave fields

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](LICENSE.txt)

## Overview

k-Wave is an open source MATLAB toolbox designed for the time-domain simulation of propagating acoustic waves in 1D, 2D, or 3D [1]. The toolbox has a wide range of functionality, but at its heart is an advanced numerical model that can account for both linear and nonlinear wave propagation, an arbitrary distribution of heterogeneous material parameters, and power law acoustic absorption.

The numerical model is based on the solution of three coupled first-order partial differential equations which are equivalent to a generalised form of the Westervelt equation [2]. The equations are solved using a k-space pseudospectral method, where spatial gradients are calculated using a Fourier collocation scheme, and temporal gradients are calculated using a k-space corrected finite-difference scheme. The temporal scheme is exact in the limit of linear wave propagation in a homogeneous and lossless medium, and significantly reduces numerical dispersion in the more general case.

Power law acoustic absorption is accounted for using a linear integro-differential operator based on the fractional Laplacian [3]. A split-field perfectly matched layer (PML) is used to absorb the waves at the edges of the computational domain. The main advantage of the numerical model used in k-Wave compared to models based on finite-difference time domain (FDTD) schemes is that fewer spatial and temporal grid points are needed for accurate simulations. This means the models run faster and use less memory. A detailed description of the model is given in the k-Wave User Manual and the references below.

1. B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields," _J. Biomed. Opt._, vol. 15, no. 2, p. 021314, 2010. [https://doi.org/10.1117/1.3360308](https://doi.org/10.1117/1.3360308)
2.  B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear ultrasound propagation in heterogeneous media with power law absorption using a k-space pseudospectral method," _J. Acoust. Soc. Am._, vol. 131, no. 6, pp. 4324-4336, 2012. [https://doi.org/10.1121/1.4712021](https://doi.org/10.1121/1.4712021)
3.  B. E. Treeby and B. T. Cox, "Modeling power law absorption and dispersion for acoustic propagation using the fractional Laplacian," _J. Acoust. Soc. Am._, vol. 127, no. 5, pp. 2741-2748, 2010. [https://doi.org/10.1121/1.3377056](https://doi.org/10.1121/1.3377056)

## Installation Instructions

The k-Wave toolbox is installed by adding the `k-Wave` folder to the MATLAB path. 

To update the path interactively, use the `Set Path` dialog box which is accessed by typing `pathtool` at the MATLAB command line. This dialog box can also be accessed using the `Set Path` button on the ribbon bar. Once the dialog box is open, the toolbox is installed by clicking `Add Folder`, selecting the k-Wave toolbox folder, and clicking `Save`. The toolbox can be uninstalled in the same fashion. 

For Linux users, using the `Set Path` dialog box requires write access to `pathdef.m`. This file can be found under `<matlabroot>/toolbox/local`. To find where MATLAB is installed, type `matlabroot` at the MATLAB command line.

To update the path programmatically, add the line 

```matlab
addpath('<pathname>/k-Wave');
```
    
to the `startup.m` file, where `<pathname>` is replaced with the location of the toolbox, and the slashes should be in the direction native to your operating system. If no `startup.m` file exists, create one, and save it in the MATLAB startup directory.

## Using The C++ Codes

Accelerated versions of the simulation functions written in C++/CUDA are also available. These are not included in the MATLAB toolbox. To use the C++/CUDA codes, the appropriate binaries (and library files if using Windows) should be downloaded from [k-wave.org/download.php](http://www.k-wave.org/download.php) and placed in the `k-Wave/binaries` folder of the toolbox.

## Documentation And Examples

After installation, you should see the k-Wave help files in the MATLAB help browser by selecting **k-Wave Toolbox** from the list of **Supplemental Software** on the contents page. If you can't see **k-Wave Toolbox** in the contents list of the MATLAB help browser, try typing `help k-Wave` at the command prompt to see if the toolbox has been installed correctly.

After installation, to make the k-Wave documentation searchable from within the MATLAB help browser, run

```matlab
builddocsearchdb(getkWavePath('helpfiles'));
```
    
Note, the created database file will only work with the version of MATLAB used to create it.

## Getting Started

Regardless of your intended application for the k-Wave Toolbox, the easiest way to get started is to work through the **Initial Value Problems** examples in the MATLAB help browser, in particular the **Homogeneous Propagation Medium** example. This gives a step-by-step introduction to the way the simulation functions within k-Wave work. Each of the examples comes with an accompanying m-file which can be opened or run from within the help file.

There is additional information on the functions and algorithms used in k-Wave in the k-Wave Manual (this can be downloaded from [k-wave.org/documentation.php](http://www.k-wave.org/documentation.php)). The manual includes a general introduction to the governing equations and numerical methods used in the main simulation functions in k-Wave. It also provides a basic overview of the software architecture and a number of canonical examples. The manual has a different emphasis to the MATLAB documentation, thus it can be beneficial when starting with k-Wave to read both in parallel.

## Getting Help

If you have any feedback or comments, or wish to report bugs, obtain help, or request new features, please visit the k-Wave forum at [k-wave.org/forum](http://www.k-wave.org/forum).

## License

The k-Wave toolbox is distributed by the copyright owners under the terms of the GNU Lesser General Public License (LGPL) which is a set of additional permissions added to the GNU General Public License (GPL). The full text of both licenses is included with the toolbox in the folder 'license'.

The licence places copyleft restrictions on the k-Wave toolbox. Essentially, anyone can use the software for any purpose (commercial or non-commercial), the source code for the toolbox is freely available, and anyone can redistribute the software (in its original form or modified) as long as the distributed product comes with the full source code and is also licensed under the LGPL. You can make private modified versions of the toolbox without any obligation to divulge the modifications so long as the modified software is not distributed to anyone else. The copyleft restrictions only apply directly to the toolbox, but not to other (non-derivative) software that simply links to or uses the toolbox. 

k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details (http://www.gnu.org/licenses/lgpl.html). 

If you find the toolbox useful for your academic work, please consider citing one or more of the following papers:

1. Overview of the toolbox with applications in photoacoustics: 

   [B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields," _J. Biomed. Opt._, 15(2), 021314, 2010.](https://doi.org/10.1117/1.3360308)

2. Nonlinear ultrasound model and the C++ code:

   [B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear ultrasound propagation in heterogeneous media with power law absorption using a k-space pseudospectral method," _J. Acoust. Soc. Am._, 131(6), 4324-4336, 2012.](https://doi.org/10.1121/1.4712021)

3. Modelling sources using the `kWaveArray` class:

   [E. S. Wise, B. T. Cox, J. Jaros, and B. E. Treeby, "Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods," _J. Acoust. Soc. Am._, 146(1), 278-288, 2019.](https://doi.org/10.1121/1.5116132)

4. Axisymmetric model:

   [B. E. Treeby, E. S. Wise, F. Kuklis, J. Jaros, and B. T. Cox, "Nonlinear ultrasound simulation in an axisymmetric coordinate system using a k-space pseudospectral method," _J. Acoust. Soc. Am._, 148(4), 2288-2300, 2020.](https://doi.org/10.1121/10.0002177)

5. Elastic wave model:

   [B. E. Treeby, J. Jaros, D. Rohrbach, B. T. Cox , "Modelling elastic wave propagation using the k-Wave MATLAB toolbox," _IEEE International Ultrasonics Symposium_, pp. 146-149, 2014.](https://doi.org/10.1109/ULTSYM.2014.0037)

6. Accelerated C++ code:

   [J. Jaros, A. P. Rendell, and B. E. Treeby, "Full-wave nonlinear ultrasound simulation on distributed clusters with applications in high-intensity focused ultrasound," _The International Journal of High Performance Computing Applications_, 30(2), 137-155, 2016.](https://doi.org/10.1177/1094342015581024)

## Contributing To k-Wave

k-Wave version 1.4.x is the last feature release of k-Wave version 1.x. Development work is now focused on k-Wave 2.x. Pull requests for k-Wave 1.x won't be accepted. Please direct any bug reports to the k-Wave forum.
