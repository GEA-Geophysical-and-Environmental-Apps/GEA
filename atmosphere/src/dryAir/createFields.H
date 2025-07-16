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

/*dimensionedScalar theta0
(
    "theta0",
    dimTemperature,
    transportProperties.lookup("theta0")
);*/

volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// aggiunta. campo gradiente di theta
volVectorField gradTheta
(
    IOobject
    (
        "gradTheta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(theta)
);

// aggiutno compo magnitude gradiente di theta
volScalarField magGradTheta
(
    IOobject
    (
        "magGradTheta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(gradTheta)
);

// aggiuntno sclare SMALL con le dimensioni fisiche giuste
const dimensionedScalar epsilon
(
    "epsilon",
    magGradTheta.dimensions(), 
    SMALL                       
);

// aggiunto il gradiente normalizzato
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    magGradTheta / (max(magGradTheta) + epsilon)
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



// aggiunto campo scalare magnetudo velocita
volScalarField magU
(
    IOobject
    (
        "magU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U)
);

// aggiuntno sclare SMALL con le dimensioni fisiche giuste
const dimensionedScalar epsilon2
(
    "epsilon2",
    magU.dimensions(), 
    SMALL                       
);

// aggiunto campo scalare magnetudo velocita normalizzata
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    magU / (max(magU) + epsilon2)
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

volScalarField Exner("Exner", Foam::pow(p/pRef, (thermo.Cp()-thermo.Cv())/thermo.Cp()));

#include "phrghEqn.H" // risolve il sistema idrostatico

//Info<< "Creating field dpdt\n" << endl; modifica, prima venica creata dpdt dentro a questo include, ora in modo esplicito
//#include "createDpdt.H"

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvOptions.H"
#include "createRadiationModel.H"
