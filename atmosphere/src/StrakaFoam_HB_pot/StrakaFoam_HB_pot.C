/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    buoyantPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of compressible fluids for
    ventilation and heat-transfer.

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pointMesh.H"                    // aggiunta perche' presente in CE1
#include "interpolation.H"                // aggiunta perche' presente in CE1
#include "volPointInterpolation.H"        // aggiunta perche' presente in CE1
#include "primitivePatchInterpolation.H"  // aggiunta perche' presente in CE1

#include "dynamicFvMesh.H"                // aggiunta per AMR classico 
#include "dynamicRefineFvMesh.H"          // aggiunta per AMR classico
#include "CorrectPhi.H"                   // aggiunta per correggere il flusso


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    //#include "createMesh.H"
    #include "createDynamicFvMesh.H"        // aggiunta. creato: mesh, reference a oggetto mesh dinamica
    #include "createDyMControls.H"          // aggiunta. crea tre boolean "correctPhi" - "checkMeshCourantNo" - "moveMeshOuterCorrectors"
    //#include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"       // aggiunta. crea surfaceVectorField rhoUf = fvc::interpolate(rho*U) 
    #include "initContinuityErrs.H"
    //#include "createTimeControls.H"
    #include "compressibleCourantNo.H" 
    #include "setInitialDeltaT.H"


    turbulence->validate();

    dimensionedScalar mass0 = fvc::domainIntegrate(rho_hyd);

    Info<< "\nStarting time loop\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    while (runTime.run())
    {   

        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"
        #include "readDyMControls.H"         // aggiunta

        // genero divrhoU per aggiustare i flussi
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

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "thetaEqn.H"


            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

        }


        rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*Foam::pow(p/pRef, thermo.Cv()/thermo.Cp());

        // --- Computation of the gradient of theta, look at fvSchemes, should be Gauss linear
        gradTheta = fvc::grad(theta);

        // --- Computation of the magnitude of gradTheta
        magGradTheta = mag(gradTheta);
        
        // --- Computation of the normalized gradTheta
        alpha = magGradTheta / (max(magGradTheta)+epsilon);	


        // ------------------ Computation enriched gradient basato sui volumi e numero di vicini--------------------------
        
        forAll(mesh.cells(), cellI)   // mesh.cells() restituisce una labelList con gli indici di ogni cella, cellI e' la label della cella corrente
        {
            const labelList& neighborCells = mesh.cellCells()[cellI];  // lista di label. indici delle celle vicine alla cella corrente

            const double numberOfNeighbor = neighborCells.size() + 1;

            const double w_j = (1/numberOfNeighbor);                   // reciproco numero di celle vicine

            scalar totalVolume = mesh.V()[cellI];                      // comincio settando il volume totale pari a quello della mesh di partenza

            auto sumR = w_j*totalVolume*gradTheta[cellI];              // setto la sommatoria di R a: w_j*|T|*gradTheta per la cella considerata


            forAll(neighborCells, indexI) // looping over the indeces of the neighbor
            {   
                const scalar currentVolume = mesh.V()[indexI];  // volume cella vicina corrente

                totalVolume = totalVolume + currentVolume;      // aggiorno il volume totale

                sumR = sumR + w_j*currentVolume*gradTheta[indexI];  // aggiorno la sommatoria di R
            }

            enrichedGrad[cellI] = (1/totalVolume)*sumR;

        }


        // ------------------ Tentativo ricostruzione della soluzione --------------------------
        
        forAll(mesh.cells(), cellI)   // mesh.cells() restituisce una labelList con gli indici di ogni cella, cellI e' la label della cella corrente
        {
            const labelList& neighborCells = mesh.cellCells()[cellI];  // lista di label. indici delle celle vicine alla cella corrente

            const double numberOfNeighbor = neighborCells.size() + 1;

            const double w_j = (1/numberOfNeighbor);                   // reciproco numero di celle vicine

            scalar totalVolume = mesh.V()[cellI];                      // comincio settando il volume totale pari a quello della mesh di partenza

            auto sumR = w_j*totalVolume*theta[cellI];              // setto la sommatoria di R a: w_j*|T|*gradTheta per la cella considerata


            forAll(neighborCells, indexI) // looping over the indeces of the neighbor
            {   
                const scalar currentVolume = mesh.V()[indexI];  // volume cella vicina corrente

                totalVolume = totalVolume + currentVolume;      // aggiorno il volume totale

                sumR = sumR + w_j*currentVolume*theta[indexI];  // aggiorno la sommatoria di R
            }

            enrichedTheta[cellI] = (1/totalVolume)*sumR;

        }
        
        
        magEnrichedGrad = mag(enrichedGrad);

        gamma = magEnrichedGrad / (max(magEnrichedGrad) + epsilon_3);

        diffGrad = gradTheta - enrichedGrad;

        diffSol = theta - enrichedTheta;

        eta = mag(diffGrad);

        eta_2 = mag(diffSol);

        delta = (mag(eta) / (max(mag(eta) )+ epsilon_4));
        


        // ------------------ mesh update block ------------------
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
        
        
        // FUNZIONE CHE GENERA .TXT PER THETA PRIMO E Uy
        
        scalar theta_primo_max = max(theta).value()-300;
        scalar Uymax = max(U.component(1)).value();

        std::ofstream file;
        file.open ("theta_primo_Uy.txt", std::ofstream::out | std::ofstream::app);
        if (Pstream::master())
        {
            file << runTime.timeName() << "\t" << theta_primo_max << "\t" << Uymax <<std::endl << "\n";
        }

        // ----------------------------------------------------------



        /* FUNZIONI PER GENERARE ERRORI OUTPUT
        dimensionedScalar mass1 = fvc::domainIntegrate(rho - rho_hyd);

        std::ofstream file;
        file.open ("err_mass.txt", std::ofstream::out | std::ofstream::app);
        if (Pstream::master())
        {
            file << mass1.value()/mass0.value() << std::endl << "\n"; //"\t" << vavg << std::endl << "\n";
        }


        scalar thetamax1 = max(theta).value();
        scalar thetamax2 = max(theta.boundaryField());
        scalar thetamax = thetamax2;
        if (thetamax1 > thetamax2)
        thetamax = thetamax1;
        scalar Uymax = max(U.component(1)).value();
        std::ofstream file2;
        file2.open ("results_BDF2.txt", std::ofstream::out | std::ofstream::app);
        if (Pstream::master())
        {
            file2 << runTime.timeName() << "\t" << thetamax << "\t" << Uymax << "\t" <<  std::endl << "\n";
        }*/


        runTime.write();

        runTime.printExecutionTime(Info);

    }
    
    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //
