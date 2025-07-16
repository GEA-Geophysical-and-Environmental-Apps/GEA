//volScalarField p0("p0", 0*pRef*pow(1 + gh/(thermo.Cp()*theta0), thermo.Cp()/(thermo.Cp()-thermo.Cv())));

if (pimple.dict().lookupOrDefault("hydrostaticInitialization", false))
{
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    volScalarField& ph_rgh_hyd = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh_hyd",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (equal(runTime.value(), 0))
    {
        p_hyd = pRef + ph_rgh_hyd;
        
        //thermo.correct();
        //scalar R = Foam::constant::thermodynamic::RR/28.9; ///molWeight_;
        //T00 = 1.0/(R*T0);
        rho_hyd = (pRef/(thermo.Cp()-thermo.Cv())/theta0)*Foam::pow(p_hyd/pRef, thermo.Cv()/thermo.Cp());

        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );

        for (label i=0; i<nCorr; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho_hyd));

            surfaceScalarField phig
            (
                "phig",
               rhof*rhof*(g & mesh.Sf())
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh_hyd, rho_hyd, U, phig, rhof);

            fvScalarMatrix ph_rgh_hydEqn
            (
                fvm::laplacian(rhof, ph_rgh_hyd) == fvc::div(phig)
            );

            ph_rgh_hydEqn.solve();

            p_hyd = pRef + ph_rgh_hyd;
            rho_hyd = (pRef/(thermo.Cp()-thermo.Cv())/theta0)*Foam::pow(p_hyd/pRef, thermo.Cv()/thermo.Cp());
            //thermo.correct();
            //rho = thermo.rho();

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh_hyd) - min(ph_rgh_hyd)).value() << endl;
        }

        //ph_rgh_hyd.write();

        //p_rgh = 0*ph_rgh;
        p0 = ph_rgh_hyd;

    }


    if (equal(runTime.value(), 0))
    {
        p = ph_rgh + pRef;
        //thermo.correct();
        //rho = thermo.rho();

        //rho = p*(1.0/((thermo.Cp() - thermo.Cv())*(theta + gh/thermo.Cp())));
        rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp());
        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );

        for (label i=0; i<nCorr; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));

            surfaceScalarField phig
            (
                "phig",
               rhof*rhof*(g & mesh.Sf())
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
            );

            ph_rghEqn.solve();

            p = ph_rgh + pRef;
            rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp());
            //thermo.correct();
            //rho = thermo.rho();

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;
        }

        //ph_rgh.write();

        p_rgh = p - p0 - pRef; //0*ph_rgh;
        //p_rgh.write();//p0 = ph_rgh;

    }

}


//volScalarField ptemp("ptemp", pRef*pow(1 + gh/(thermo.Cp()*theta0), thermo.Cp()/(thermo.Cp()-thermo.Cv())));

//p_rgh = 

//volScalarField p0("p0", pRef*pow(1 + gh/(thermo.Cp()*theta0), thermo.Cp()/(thermo.Cp()-thermo.Cv())));
//p = p0;
//thermo.correct();
//volScalarField rho0("rho0", thermo.rho());


//theta = thermo.T() - gh/thermo.Cp();
//volScalarField ptemp("ptemp", pRef*pow(1 + gh/(thermo.Cp()*theta), thermo.Cp()/(thermo.Cp()-thermo.Cv())));

//p = ptemp;
//thermo.correct();
//rho = thermo.rho();

//p_rgh = ptemp - p0;

//volScalarField T0("T0", theta0 + gh/thermo.Cp());
//volScalarField rho0("rho0", p0/(thermo.Cp()-thermo.Cv())/T0);

//p0.write();
//T0.write();
//rho0.write();

//p = p0 + p_rgh;
//thermo.correct();
//rho = thermo.rho();

/*theta = thermo.T() - gh/thermo.Cp();
volScalarField ptemp("ptemp", pRef*pow(1 + gh/(thermo.Cp()*theta), thermo.Cp()/(thermo.Cp()-thermo.Cv())));

p = ptemp;
thermo.correct();
rho = thermo.rho();

p_rgh = ptemp - p0;*/


