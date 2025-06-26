    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);


    /*scalar yu1 = 15000;
    scalar yu2 = 19500;
    scalar xl1 = -25000;
    scalar xl2 = -50000;
    scalar xd1 = 25000;
    scalar xd2 = 50000;
    scalar yd = 15000;*/

    //SPONGE LAYER
   /*forAll(ycord,faceI)
   {           
       ycord[faceI] = mesh.C()[faceI].component(1);
   } 

   forAll(xcord,faceI)
   {           
       xcord[faceI] = mesh.C()[faceI].component(0);
   }   


   forAll(ycord, faceII)
   {

       if (ycord[faceII] > yu1)
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

       if (xcord[faceII] > xd1 && ycord[faceII] < 15000)
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

       if (xcord[faceII] < xl1 && ycord[faceII] < 15000)
       {           
             Cx2[faceII] = 1.0; 
       }  

       else
       {           
             Cx2[faceII] = 0.0; 
       }  
   }
      
*/

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
       //- fvm::laplacian(rho*visc*oneLength*oneLength/oneTime, U)
     ==
       g*(rho - rho_hyd) //- rho*(Da + Dg)*(fvm::Sp(1.0/oneTime/oneLength,U) - vector(Uref.value(), 0, 0)*oneLength/oneTime/oneTime/oneLength) - rho*Dx1*(fvm::Sp(1.0/oneTime/oneLength,U) - vector(Uref.value(), 0, 0)*oneLength/oneTime/oneTime/oneLength) - rho*Dx2*(fvm::Sp(1.0/oneTime/oneLength,U) - vector(Uref.value(), 0, 0)*oneLength/oneTime/oneTime/oneLength)
//+ fvm::Sp(-Cy*rho*(alpha*Foam::sin(1.57*(ycord - yu1)/(yu2 - yu1))*Foam::sin(1.57*(ycord - yu1)/(yu2 - yu1))), U) //
//+ Cy*rho*Uref*(alpha*Foam::sin(1.57*(ycord - yu1)/(yu2 - yu1))*Foam::sin(1.57*(ycord - yu1)/(yu2 - yu1))) //
//-  Cx1*rho*(U - Uref)*(alpha*Foam::sin(1.57*(xcord - xd1)/(xd2 - xd1))*Foam::sin(1.57*(xcord -  xd1)/(xd2 - xd1))) -  ////Cx2*rho*(U - Uref)*(alpha*Foam::sin(1.57*(xcord - xl1)/(xl2-xl1))*Foam::sin(1.57*(xcord - xl1)/(xl2 - xl1)))
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
         == //-fvc::grad(p)
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
