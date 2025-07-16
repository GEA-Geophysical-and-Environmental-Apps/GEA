{
    //volScalarField& he = thermo.he();

    fvScalarMatrix ThetaEqn
    (
        fvm::ddt(rho, theta) + fvm::div(phi, theta)
      //+ fvc::ddt(rho, K) + fvc::div(phi, K)
      //+ (
      //      he.name() == "e"
      //    ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
      //     : -dpdt
     //   )
      //+ thermophysicalTransport->divq(he)
        - fvm::laplacian(turbulence->alphaEff(), theta)
       //- fvm::laplacian(rho*visc*oneLength*oneLength/oneTime, theta)
     //==
      //- rho*(Da + Dg)*(fvm::Sp(1.0/oneTime/oneLength, theta) - theta0/oneTime/oneLength) - rho*Dx1*(fvm::Sp(1.0/oneTime/oneLength, theta) - theta0/oneTime/oneLength) - rho*Dx2*(fvm::Sp(1.0/oneTime/oneLength,theta) - theta0/oneTime/oneLength)
    //    rho*(U&g)
    //  + fvOptions(rho, he)
    );

    ThetaEqn.relax();

    //fvOptions.constrain(EEqn);

    ThetaEqn.solve();

    //fvOptions.correct(he);

    //he.write();
    //he.oldTime().write();

    //thermo.correct();
    //T = theta + gh/thermo.Cp() + theta0;
    //psi = 1.0/((thermo.Cp() - thermo.Cv())*(theta + gh/thermo.Cp()));
    psi = (1.0/(thermo.Cp() - thermo.Cv())/theta)*Foam::pow(pRef/p,(1 - thermo.Cv()/thermo.Cp()));
            //rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp());
    

    /*
    // --- Computation of the gradient of theta
    gradTheta = fvc::grad(theta);

    // --- Computation of the magnitude of gradTheta
    magGradTheta = mag(gradTheta);
    
    // --- Computation of the normalized gradTheta
    alpha = magGradTheta / (max(magGradTheta)+epsilon);

    // --- mesh size parameter
    const dimensionedScalar delta_quadro
    (
      "delta_quadro",
      dimensionSet(0, 2, 0, 0, 0, 0, 0),  
      10  // value               
    );

    
    // --- computation for thetaMean
    fvScalarMatrix ThetaMeanEqn
    (
      -fvm::laplacian(delta_quadro*alpha, thetaMean) + fvm::Sp(1.0,thetaMean) - theta
    );
    
    //ThetaMeanEqn.relax();
    ThetaMeanEqn.solve();
    */

}
