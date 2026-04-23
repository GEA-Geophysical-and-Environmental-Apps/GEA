#include "fvCFD.H"

Info<< "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Re
(
    "Re",
    physicalProperties.lookup("Re")
);

dimensionedScalar Ro
(
    "Ro",
    physicalProperties.lookup("Ro")
);

dimensionedScalar Fr
(
    "Fr",
    physicalProperties.lookup("Fr")
);

dimensionedScalar sigma
(
    "sigma",
    physicalProperties.lookup("sigma")
);

dimensionedScalar delta
(
    "delta",
    physicalProperties.lookup("delta")
);

Info<< "Reading field q1\n" << endl;
volScalarField q1
(
    IOobject
    (
        "q1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field q2\n" << endl;
volScalarField q2
(
    IOobject
    (
        "q2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field psi1\n" << endl;
volScalarField psi1
(
    IOobject
    (
        "psi1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field psi2\n" << endl;
volScalarField psi2
(
    IOobject
    (
        "psi2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField zero
(
    IOobject
    (
        "zero",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0.0)
);

const vector k(0, 0, 1);

// Forcing Term
dimensionedScalar one("one", dimensionSet(0,1,0,0,0,0,0), 1.0);
volScalarField forcingTerm1("forcingTerm1", sin(constant::mathematical::pi*mesh.C().component(1)));

#include "createFvOptions.H"
