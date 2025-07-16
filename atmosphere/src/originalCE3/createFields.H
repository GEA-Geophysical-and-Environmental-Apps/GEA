//#include "createRDeltaT.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();
volScalarField& p_hyd = thermo.p();
volScalarField& T = thermo.T();

/*dimensionedScalar theta0
(
    "theta0",
    dimTemperature,
    transportProperties.lookup("theta0")
);*/

/*dimensionedScalar alpha
(
    "alpha",
    //dimTemperature,
    transportProperties.lookup("alpha")
);
*/

/*dimensionedScalar Uref
(
    "Uref",
    //dimTemperature,
    transportProperties.lookup("Uref")
);
*/

/*volScalarField ycord
(
    IOobject
    (
        "ycord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField xcord
(
    IOobject
    (
        "xcord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/

/*volVectorField Uref
(
    IOobject
    (
        "Uref",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Cy
(
    IOobject
    (
        "Cy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Cx1
(
    IOobject
    (
        "Cx1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Cx2
(
    IOobject
    (
        "Cx2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField theta0
(
    IOobject
    (
        "theta0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField rho_hyd
(
    IOobject
    (
        "rho_hyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

/*Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);*/


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p0
(
    IOobject
    (
        "p0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}
mesh.setFluxRequired(p_rgh.name());


/*hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);*/

#include "phrghEqn.H"


Info<< "Creating field dpdt\n" << endl;
/*volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);*/

#include "createDpdt.H"

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvOptions.H"
//#include "createFvConstraints.H"
#include "createRadiationModel.H"

dimensionedScalar oneT
(
    "oneT",
    dimensionSet(0, 0, 0, 1, 0, 0 ,0),
    1.0
);

dimensionedScalar oneP
(
    "oneP",
    dimensionSet(1, -1, -2, 0, 0, 0 ,0),
    1.0
);

dimensionedScalar oneRho
(
    "oneRho",
    dimensionSet(1, -3, 0, 0, 0, 0 ,0),
    1.0
);

dimensionedScalar oneCp
(
    "oneCp",
    dimensionSet(0, 2, -2, -1, 0, 0 ,0),
    1.0
);

dimensionedScalar oneTime
(
    "oneTime",
    dimensionSet(0, 0, 1, 0, 0, 0 ,0),
    1.0
);

dimensionedScalar oneLength
(
    "oneLength",
    dimensionSet(0, 1, 0, 0, 0, 0, 0),
    1.0
);

dimensionedScalar sourceTheta
(
    "sourceTheta",
    dimensionSet(0, 0, -1, 1, 0, 0 ,0),
    1.0
);



volScalarField xcord ("xcord", mesh.C().component(0));
volScalarField ycord ("ycord", mesh.C().component(1));
volScalarField Dx1 ("Dx1", xcord*0.0);
volScalarField Dx2 ("Dx2", xcord*0.0);
volScalarField Da ("Da", ycord*0.0);
volScalarField Dg ("Dg", ycord*0.0);

volScalarField visc ("visc", ycord*0.0/oneLength);
volScalarField visc2 ("visc2", ycord*0.0/oneLength);


dimensionedScalar yu1
(
    "yu1",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("yu1")
    //1.0
);

dimensionedScalar yu2
(
    "yu2",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("yu2")
    //1.0
);

dimensionedScalar xd1
(
    "xd1",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("xd1")
    //1.0
);

dimensionedScalar xd2
(
    "xd2",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("xd2")
    //1.0
);

dimensionedScalar xs1
(
    "xs1",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("xs1")
    //1.0
);

dimensionedScalar xs2
(
    "xs2",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("xs2")
    //1.0
);

dimensionedScalar Uref
(
    "Uref",
    dimensionSet(0, 0, 0, 0, 0, 0 ,0),
    transportProperties.lookup("Uref")
);



//scalar yu1 = 15000; //15000;
//scalar yu2 = 30000; //19500;
//scalar xd1 = 40000; //32000; //40000;
//scalar xs1 = -40000; //32000; //40000;
//scalar xd2 = 120000; //72000; //120000;
//scalar xs2 = -120000; //-72000; //-120000;

scalar dx = 500;
//scalar dz = 0.2*(yu2 - yu1);
scalar alpha = 0.1;
scalar alpha2 = 0.1;
//SPONGE LAYER

scalar Kvmax = 3;
scalar Kvmin = 1;

scalar za = 0;
scalar zb = 20000;



//forAll(visc2, face_visc)
   //{

       /*if ((ycord[face_visc]) > 0  && (ycord[face_visc]) <= za)
       {           
             visc[face_visc] = Kvmax;
       } */ 

       //if ((ycord[face_visc]) > za  && (ycord[face_visc]) <= zb) // && (xcord[face_visc]) > xs1 && (xcord[face_visc]) < xd1)
       //{           
      //       visc[face_visc] = Kvmax*((zb - ycord[face_visc])/(zb - za)) + Kvmin*((ycord[face_visc] - za)/(zb - za));
       //}  
  // }

/*forAll(visc, face_visc)
{

       if ((ycord[face_visc]) > 0  && (ycord[face_visc]) <= za)
       {           
             visc[face_visc] = Kvmax;
       }  

       else if ((ycord[face_visc]) > za  && (ycord[face_visc]) <= zb && (xcord[face_visc]) > xs1 && (xcord[face_visc]) < xd1)
       {           
            visc[face_visc] = Kvmax*((zb - ycord[face_visc])/(zb - za)) + Kvmin*((ycord[face_visc] - za)/(zb - za));
       }  
 }*/


forAll(Dg, faceII)
   {

         //if ((ycord[faceII] - yu1)/(yu2 - yu1) >= 0 && (ycord[faceII] - yu1)/(yu2 - yu1) <= 0.5)
       if ((ycord[faceII]) >= yu1.value())
       {           
             //Dg[faceII] = (alpha/2)*(1 - Foam::cos((ycord[faceII] - yu1)*3.14/(yu2 - yu1)));
          Dg[faceII] = alpha*Foam::sin(1.57*(ycord[faceII] - yu1.value())/(yu2.value() - yu1.value()))*Foam::sin(1.57*(ycord[faceII] - yu1.value())/(yu2.value() - yu1.value()));
       }  

       else //if ((ycord[faceII] - yu1)/(yu2 - yu1) > 0.5 && (ycord[faceII] - yu1)/(yu2 - yu1) <= 1)
       {           
             Dg[faceII] = 0.0; //Dg[faceII] = (alpha/2)*(1 + Foam::cos(((ycord[faceII] - yu1)/(yu2 - yu1) - 0.5)*3.14));
       }  
   }


/*forAll(Da, faceIII)
   {
       if (ycord[faceIII] > yu2 - 3*dz)
                      
             Da[faceIII] = (1 - Foam::tanh((yu2 - ycord[faceIII])/dz))/Foam::tanh((yu2 - ycord[faceIII])/dz);
 
   }

   */
forAll(Dx1, faceIx)
   {

       //if ((xcord[faceIx] - xd1)/(xd2 - xd1) >= 0 && (xcord[faceIx] - xd1)/(xd2 - xd1) <= 0.5)
       if ((xcord[faceIx]) >= xd1.value())
       {           
             Dx1[faceIx] =  alpha*Foam::sin(1.57*(xcord[faceIx] - xd1.value())/(xd2.value() - xd1.value()))*Foam::sin(1.57*(xcord[faceIx] - xd1.value())/(xd2.value() - xd1.value()));
       }  

       else //if ((xcord[faceIx] - xd1)/(xd2 - xd1) > 0.5 && (xcord[faceIx] - xd1)/(xd2 - xd1) <= 1)
       {           
             Dx1[faceIx] = 0.0;
       }  

       /*if (xcord[faceIx] <= xs1)
                      
             Dx1[faceIx] = (1 - Foam::tanh((xs2 - xcord[faceIx])/dx))/Foam::tanh((xs2 - xcord[faceIx])/dx);
       */
  }


forAll(Dx2, faceIIx)
   {


        //if ((xcord[faceIIx] - xs1)/(xs2 - xs1) >= 0 && (xcord[faceIIx] - xs1)/(xs2 - xs1) <= 0.5)
      if ((xcord[faceIIx]) <= xs1.value())
       {           
             Dx2[faceIIx] =  alpha*Foam::sin(1.57*(xcord[faceIIx] - xs1.value())/(xs2.value() - xs1.value()))*Foam::sin(1.57*(xcord[faceIIx] - xs1.value())/(xs2.value() - xs1.value()));
       }   

       else //if ((xcord[faceIIx] - xs1)/(xs2 - xs1) > 0.5 && (xcord[faceIIx] - xs1)/(xs2 - xs1) <= 1)
       {           
             Dx2[faceIIx] = 0.0; //(alpha2/2)*(1 + Foam::cos(((xcord[faceIIx] - xs1)/(xs2 - xs1) - 0.5)*3.14));
       }  
      /*if (xcord[faceIIx] >= xd1)
                      
             Dx2[faceIIx] = (1 - Foam::tanh((xd2 - xcord[faceIIx])/dx))/Foam::tanh((xd2 - xcord[faceIIx])/dx);
       */
   }


