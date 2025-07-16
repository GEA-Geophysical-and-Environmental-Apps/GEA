/*---------------------------------------------------------------------------*\
License
    This file is part of GEA

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
    dryAir

Description
    Transient solver for compressible Navier-Stokes equations with gravity

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pointMesh.H"
#include "interpolation.H"
#include "volPointInterpolation.H"
#include "primitivePatchInterpolation.H"

#include "dynamicFvMesh.H"              // aggiunta per AMR classico 
#include "dynamicRefineFvMesh.H"        // aggiunta per AMR classico
#include "CorrectPhi.H"                 // aggiunta per correggere il flusso


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"                 // creato: runTime
    #include "addCheckCaseOptions.H"
    //#include "createMesh.H"               // modifica
    #include "createDynamicFvMesh.H"        // aggiunta. creato: mesh, reference a oggetto mesh dinamica
    #include "createDyMControls.H"          // aggiunta. crea tre boolean "correctPhi" - "checkMeshCourantNo" - "moveMeshOuterCorrectors"
    #include "initContinuityErrs.H"
    //#include "createControl.H"            // modifica
    #include "createFields.H"
    #include "createFieldRefs.H"            // creato: psi (come const object!)
    //#include "createTimeControls.H"       // modifica
    #include "createRhoUfIfPresent.H"       // aggiunta. crea surfaceVectorField rhoUf = fvc::interpolate(rho*U) 
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"
        #include "readDyMControls.H"         // aggiunta
         
        // aggiunto blocco che calcola divrhoU. originariamente era un autoPtr<volScalarField> divrhoU;
        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence

        volScalarField divrhoU
        (
            IOobject
            (
                "divrhoU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::div(fvc::absolute(phi, rho, U)) 
        );

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //#include "rhoEqn.H" spostata nella sola prima iterazione del pimple

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {   

            if(pimple.firstIter())
            {
                #include "rhoEqn.H"
            }
            
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"

                /*
                Info << "Magnitude of rAU: " << max(mag(rAU)) << endl;
                Info << "Magnitude of psip0: " << max(mag(psip0)) << endl;
                Info << "Magnitude of rhorAUf: " << max(mag(rhorAUf)) << endl; 
                Info << "Magnitude of HbyA: " << max(mag(HbyA)) << endl; 
                Info << "Magnitude of phiHbyA: " << max(mag(phiHbyA)) << endl;
                Info << "Magnitude of phig: " << max(mag(phig)) << endl;
                */  
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

        }

        rho = thermo.rho(); 

        // Update Exner function
        Exner = Foam::pow(p/pRef, (thermo.Cp() - thermo.Cv())/thermo.Cp());

        // --- Computation of the potential temperature field
        theta = thermo.T()/Exner;

        // --- Computation of the gradient of theta
        gradTheta = fvc::grad(theta);
        
        // --- Computation of the magnitude of gradTheta
        magGradTheta = mag(gradTheta);
        
        // --- Computation of the normalized gradTheta
        alpha = magGradTheta / (max(magGradTheta)+epsilon);	

        // --- Computatino of the magnitude of velocity
        magU = mag(U);

        // --- Computation of beta
        beta = magU / (max(magU) + epsilon2);
        
        /*
        Info << "Magnitude of pressure: " << max(mag(p)) << endl;
        Info << "Magnitude of velocity: " << max(mag(U)) << endl;
        Info << "Magnitude of density: " << max(mag(rho)) << endl;
        Info << "Magnitude of psi: " << max(mag(psi))<< endl;
        Info << "Magnitude of p_rgh: " << max(mag(p_rgh)) << endl;
        Info << "Magnitude of phi: " << max(mag(phi)) << endl;
        Info << "Magnitude of Exner: " << max(mag(Exner)) << endl;
        Info << "Magnitude of Theta: " << max(mag(theta)) << endl;
        Info << "Magnitude of gradTheta: " << max(mag(gradTheta)) << endl;
        */
        
        

        // ------------------ mesh update block -------------
        volVectorField rhoU
        (
            IOobject
            (
                "rhoU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rho*U
        );
            
        //fvModels.preUpdateMesh(); da errore, in effetti e' models

        mesh.update();

        if (mesh.changing())
        {

            gh = (g & mesh.C()) - ghRef;
            ghf = (g & mesh.Cf()) - ghRef;  

            MRF.update();
            
            if (correctPhi)
            {

                #include "correctPhi.H"                    

            }
        }
        // ------------------------------------------------------
        

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
