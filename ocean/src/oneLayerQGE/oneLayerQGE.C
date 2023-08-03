/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    myQGESolver

Description
    This is a transient solver for the one-layer quasigeostrophic equations:
        \[
            ddt(q) 
            + div(curl(psiVector)q)
            - 1/Re \Delta q 
            = F
        \] 
    where F = sin(pi*y). The forcing term F may be changed in CreateFields.H. 
    The corticity and stream function are related throught the kinetic relationship:
        \[
            -Ro*laplacian(psi) 
            + y
            = q
        \]
    q : non-dimensional potential vorticity
    psi : non-dimensional stream function
    Re : Reynolds number
    Ro : Rossby number
    F : forcing term

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vector.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const volVectorField& C = mesh.C();

    Info<< "\nSolving the one-layer quasi-geostrophic equations\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Computing potential vorticity: ddt(q) + \nabla \cdot ((\nabla\times\psi)q) - 1/Re \Delta q = F where F = sin(pi*y)

        psiVector = i*zero + j*zero + k*psi; // creating a vector (0, 0, psi)

        #include "createPhiCurl.H"

        fvScalarMatrix qEqn 
        (
            fvm::ddt(q)
          + fvm::div(phi, q)
          - fvm::laplacian(1./Re, q)
         ==
            forcingTerm
          + fvOptions(q)
        );

        qEqn.relax();

        fvOptions.constrain(qEqn);

        qEqn.solve();

        fvOptions.correct(q);

        // Computing stream function: q = -Ro\Delta\psi + y

        fvScalarMatrix psiEqn
        (
          - fvm::laplacian(Ro,psi)
          + C().component(1)
        );

        solve(psiEqn == q);
        
        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
