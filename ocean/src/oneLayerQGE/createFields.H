/*---------------------------------------------------------------------------*\
                           ██████╗ ███████╗ █████╗
                          ██╔════╝ ██╔════╝██╔══██╗
                          ██║  ███╗█████╗  ███████║
                          ██║   ██║██╔══╝  ██╔══██║
                          ╚██████╔╝███████╗██║  ██║
                           ╚═════╝ ╚══════╝╚═╝  ╚═╝

                 GEA: Geophysical and Environmental Applications

    Author:        Michele Girfoglio
    Institution:   University of Palermo

------------------------------------------------------------------------------
License
    This file is part of GEA.

    GEA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

Info<< "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Re
(
    "Re",
    physicalProperties.lookup("Re")
);

dimensionedScalar Ro
(
    "Ro",
    physicalProperties.lookup("Ro")
);

Info<< "Reading field q\n" << endl;
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField zero
(
    IOobject
    (
        "zero",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0.0)
);

const vector i(1, 0, 0);
const vector j(0, 1, 0);
const vector k(0, 0, 1);

volVectorField psiVector
(
    IOobject
    (
        "psiVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimless, vector(0, 0, 0))
);

dimensionedScalar one("one", dimensionSet(0,1,0,0,0,0,0), 1.0);
volScalarField forcingTerm("forcingTerm", sin(constant::mathematical::pi*mesh.C().component(1)/one));

#include "createPhiCurl.H"
#include "createFvOptions.H"
