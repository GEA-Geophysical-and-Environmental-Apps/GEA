/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    createPhiCurl

Description
    Creates and initialises the interpolation of the curl of the stream
    function on the faces of the control volumes and enforcing null normal
    velocity on the boundary patches

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Creating a vector of the form (0, 0, psi) for the computation of the curl

Info<< "Calculating face flux field phi\n" << endl;

volVectorField curlPsi(fvc::curl(psiVector));

// finding patch ID
label vertWallsIndex = mesh.boundaryMesh().findPatchID("verticalWalls");
label horWallsIndex = mesh.boundaryMesh().findPatchID("horizontalWalls");

// forcing null velocity on the boundary patch
forAll(mesh.boundaryMesh()[vertWallsIndex], faceI)
{
    curlPsi.boundaryFieldRef()[vertWallsIndex][faceI].component(0) = 0;
}

forAll(mesh.boundaryMesh()[horWallsIndex], faceI)
{
    curlPsi.boundaryFieldRef()[horWallsIndex][faceI].component(1) = 0;
}

surfaceScalarField phi("phi", fvc::flux(curlPsi));


// ************************************************************************* //
