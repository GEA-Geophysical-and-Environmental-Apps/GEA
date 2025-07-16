    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);


    //SPONGE LAYER
   forAll(ycord,faceI)
   {           
       ycord[faceI] = mesh.C()[faceI].component(1);
   } 

   forAll(xcord,faceI)
   {           
       xcord[faceI] = mesh.C()[faceI].component(0);
   }   


   forAll(ycord, faceII)
   {

       if (ycord[faceII] > 25000)
       {           
             Cy[faceII] = 1.0; 
       }  

       else
       {           
             Cy[faceII] = 0.0; 
       }  
   }

forAll(xcord, faceII)
   {

       if (xcord[faceII] > 110000)
       {           
             Cx1[faceII] = 1.0; 
       }  

       else
       {           
             Cx1[faceII] = 0.0; 
       }  
   }


forAll(xcord, faceII)
   {

       if (xcord[faceII] < -110000)
       {           
             Cx2[faceII] = 1.0; 
       }  

       else
       {           
             Cx2[faceII] = 0.0; 
       }  
   }
      


    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
       g*(rho - rho_hyd) 
- Cy*rho*(U - Uref)*(alpha*Foam::sin(0.5*(ycord - 25000)/(30000 - 25000))*Foam::sin(0.5*(ycord - 25000)/(30000 - 25000))) - Cx1*rho*(U - Uref)*(alpha*Foam::sin(0.5*(xcord - 110000)/(120000 - 110000))*Foam::sin(0.5*(xcord -  110000)/(120000 - 110000))) -  Cx2*rho*(U - Uref)*(alpha*Foam::sin(0.5*(xcord + 110000)/(-120000 + 110000))*Foam::sin(0.5*(xcord + 110000)/(-120000 + 110000)))
//fvOptions(rho, U)
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
