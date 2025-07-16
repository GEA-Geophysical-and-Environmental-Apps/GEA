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
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H" 
    #include "setInitialDeltaT.H"

    turbulence->validate();


    int ind = 0;
    double err_mass = 0;
    //ass0 = 0;
    //double energy0 = 0;
    //double err_energy = 0;

//E0 = rho*thermo.he() - p + rho*gh + rho*K;


 /*forAll (rho_hyd.internalField(), cellI)
        {
             mass0 += (rho_hyd[cellI])*mesh.V()[cellI];
             
             //Volume += mesh.V()[cellI];
        }
*/

dimensionedScalar mass0 = fvc::domainIntegrate(rho_hyd);

    Info<< "\nStarting time loop\n" << endl;

    Info<< "QUESTO E' IL SOVER CE3 ORIGINALE SENZA NESSUNA MODIFICA"<<endl;
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            //#include "EEqn.H"
            #include "thetaEqn.H"
       //theta = theta - (Da + Dg)*(theta - theta0)/oneLength - Dx1*(theta - theta0)/oneLength - Dx2*(theta - theta0)/oneLength;

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

         //rho = psi*p; // rho = thermo.rho();
       rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*Foam::pow(p/pRef, thermo.Cv()/thermo.Cp());

/*forAll (rho.internalField(), cellI)
        {
             err_mass += (rho[cellI] -rho_hyd[cellI])*mesh.V()[cellI];
             //err_energy += (rho[cellI]*he[cellI] - p[cellI] + rho[cellI]*gh[cellI] + rho[cellI]*K[cellI] -E0[cellI])*mesh.V()[cellI];
        }

       reduce(err_mass, sumOp<scalar> ());
*/

       dimensionedScalar mass1 = fvc::domainIntegrate(rho - rho_hyd);
       //err_mass = mass1/mass0;
       //reduce(err_energy, sumOp<scalar> ());
       //err_energy = err_energy/energy0;

       //err_emergy


std::ofstream file;
        file.open ("err_mass.txt", std::ofstream::out | std::ofstream::app);
if (Pstream::master())
        {
            file << mass1.value()/mass0.value() << std::endl << "\n"; //"\t" << vavg << std::endl << "\n";
        }
       //U = U - (Da + Dg)*(U - vector(20, 0, 0)*oneLength/oneTime)/oneLength - Dx1*(U - vector(20, 0, 0)*oneLength/oneTime)/oneLength - Dx2*(U - vector(20, 0, 0)*oneLength/oneTime)/oneLength;
        
        //theta = thermo.T() - gh/thermo.Cp() - theta0;


        //scalar vmax = max(mag(U.component(1))).value();

        //rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*Foam::pow(p/pRef, thermo.Cv()/thermo.Cp());



/*ind = ind + 1;
if (ind == 1 || ind == 43200 || ind == 86400 || ind == 129600 || ind ==  172800 || ind == 216000 || ind ==  259200 || ind == 302400 || ind == 345600 || ind == 388800 || ind ==  432000 || ind == 475200 || ind ==  518400|| ind == 561600 || ind ==  604800 || ind == 648000 || ind == 691200 || ind == 734400 || ind ==  777600 || ind == 820800 || ind ==  864000 || ind == 907200 || ind ==  950400 || ind == 993600 || ind ==  1036800 || ind == 1080000 || ind ==  1123200 || ind == 1166400 || ind ==  1209600 || ind == 1252800 || ind == 1296000 || ind == 1339200 || ind ==  1382400 || ind == 1425600 || ind ==  1468800 || ind == 1512000 || ind ==  1555200 || ind == 1598400 || ind ==  1641600 || ind == 1684800 || ind ==  1728000 || ind == 1771200 || ind == 1814400 || ind == 1857600 || ind ==  1900800 || ind == 1944000 || ind ==  1987200 || ind == 2030400 || ind ==  2073600 || ind == 2160000)
{
        
std::ofstream file;
        file.open ("results.txt", std::ofstream::out | std::ofstream::app);
if (Pstream::master())
        {
            file << vmax << std::endl << "\n"; //"\t" << vavg << std::endl << "\n";
        }

}
*/

//scalar vmax = max(U.component(1)).value();
//scalar thetamax = max(theta).value();
        /*scalar vavg = 0;
        scalar Volume = 0;

        forAll (U.internalField(), cellI)
        {
             vavg += (U[cellI] & vector (0, 1, 0))*mesh.V()[cellI];
             Volume += mesh.V()[cellI];
        }

        reduce(vavg, sumOp<scalar> ());
        reduce(Volume, sumOp<scalar> ());

        vavg = vavg/Volume;*/
        /*scalar nuttt = 0;
        scalar alphattt = 0;
        scalar Volume = 0;

const volScalarField& nutt  = turbulence -> nut();
//const volScalarField& alphatt  = turbulence -> alphat();

        forAll (nutt.internalField(), cellI)
        {
             nuttt += rho[cellI]*nutt[cellI]*mesh.V()[cellI];
             //alphattt += alphatt[cellI]*mesh.V()[cellI];
             Volume += mesh.V()[cellI];
        }

        reduce(nuttt, sumOp<scalar> ());
        //reduce(alphattt, sumOp<scalar> ());
        reduce(Volume, sumOp<scalar> ());

        nuttt = nuttt/Volume;*/
        //alphattt = alphattt/Volume;


        
//std::ofstream file;
//        file.open ("results.txt", std::ofstream::out | std::ofstream::app);
//if (Pstream::master())
//       {
//            file << runTime.timeName() << "\t" << thetamax << "\t" << vmax << std::endl << "\n";
//        }

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
        }

        scalar theta_primo_max = max(theta).value()-300;
        scalar Uymax_new = max(U.component(1)).value();

        std::ofstream file_new;
        file_new.open ("theta_primo_Uy.txt", std::ofstream::out | std::ofstream::app);
        if (Pstream::master())
        {
            file_new << runTime.timeName() << "\t" << theta_primo_max << "\t" << Uymax_new <<std::endl << "\n";
        }
//IDRO
        /*scalar vmax = max(U.component(1)).value();
        scalar vavg = 0;
        scalar Volume = 0;

        forAll (U.internalField(), cellI)
        {
             vavg += (U[cellI] & vector (0, 1, 0))*mesh.V()[cellI];
             Volume += mesh.V()[cellI];
        }

        reduce(vavg, sumOp<scalar> ());
        reduce(Volume, sumOp<scalar> ());

        vavg = vavg/Volume;*/
        /*scalar nuttt = 0;
        scalar alphattt = 0;
        scalar Volume = 0;

const volScalarField& nutt  = turbulence -> nut();
//const volScalarField& alphatt  = turbulence -> alphat();

        forAll (nutt.internalField(), cellI)
        {
             nuttt += rho[cellI]*nutt[cellI]*mesh.V()[cellI];
             //alphattt += alphatt[cellI]*mesh.V()[cellI];
             Volume += mesh.V()[cellI];
        }

        reduce(nuttt, sumOp<scalar> ());
        //reduce(alphattt, sumOp<scalar> ());
        reduce(Volume, sumOp<scalar> ());

        nuttt = nuttt/Volume;*/
        //alphattt = alphattt/Volume;


        
/*std::ofstream file;
        file.open ("results.txt", std::ofstream::out | std::ofstream::app);
if (Pstream::master())
        {
            file << runTime.timeName() << "\t" << vmax << "\t" << vavg << std::endl << "\n";
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }*/


/*scalar thetamax1 = max(theta).value();
scalar thetamax2 = max(theta.boundaryField());
scalar thetamax = thetamax2;
if (thetamax1 > thetamax2)
thetamax = thetamax1;
scalar Uymax = max(U.component(1)).value();
std::ofstream file;
        file.open ("results.txt", std::ofstream::out | std::ofstream::app);
if (Pstream::master())
        {
            file << runTime.timeName() << "\t" << thetamax << "\t" << Uymax << "\t" <<  std::endl << "\n";
        }*/

        runTime.write();

        runTime.printExecutionTime(Info);

}
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
