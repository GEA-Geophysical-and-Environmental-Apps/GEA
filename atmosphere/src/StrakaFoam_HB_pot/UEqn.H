    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)

     ==
       g*(rho - rho_hyd)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         == 
            fvc::reconstruct
            (
                (
                  //- ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh) 
                )*mesh.magSf()

		//+fvc::interpolate(rho - rho_hyd)*(g & mesh.Sf())
            )
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
