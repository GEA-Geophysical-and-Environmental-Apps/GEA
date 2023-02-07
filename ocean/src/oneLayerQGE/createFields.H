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

Info<< "Reading field q\n" << endl;
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
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

const vector i(1, 0, 0);
const vector j(0, 1, 0);
const vector k(0, 0, 1);

volVectorField psiVector
(
    IOobject
    (
        "psiVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimless, vector(0, 0, 0))
);

#include "createPhiCurl.H"
#include "createFvOptions.H"