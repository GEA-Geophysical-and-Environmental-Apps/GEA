/*---------------------------------------------------------------------------*\
GEOPHYSICAL AND ENVIRONMENTAL APPLICATIONS (GEA)
License
    This file is part of GEA (https://github.com/GEA-Geophysical-and-Environmental-Apps)
    GEA is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GEA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

Application
    twoLayerQGE

Description
    This is a transient solver for the two-layer quasigeostrophic equations:
        \[
            ddt(q1)
            + div((curl(psi1))q1)
            + (Fr/Re*delta)*laplacian(psi2 - psi1)
            - (1/Re)*laplacian(q1)
            = F1
        \]
    where F1 = sin(pi*y). The forcing term F1 may be changed in CreateFields.H
    or by using fvOptions,
        \[
            ddt(q2)
            + div((curl(psi2))q2)
            + (Fr/Re*(1-delta))*laplacian(psi1 - psi2)
            - (1/Re)*laplacian(q2)
            + sigma*laplacian(psi2)
            = F2
        \]
    with the kinematic relationships:
        \[
            Ro*laplacian(psi1)
            + y
            + Fr/delta*(psi2 - psi1)
            = q1
        \]
    and
        \[
            Ro*laplacian(psi2)
            + y
            + Fr/(1-delta)*(psi1 - psi2)
            = q2
        \]
    where
    q_i : non-dimensional potential vorticity of the ith layer 
    psi_i : stream function on the ith layer 
    Re : Reynolds number
    Ro : Rossby number
    Fi : forcing term on the ith layer
    Fr : Froude number
    delta : aspect ratio of two layer's thicknesses H1/(H1+H2)

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
    #include "createPhiCurl.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const volVectorField& C = mesh.C();         // extracting mesh information

    Info<< "\nSolving the two-layer quasi-geostrophic equations\n" << endl;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "2QGEsolver.C"

        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
