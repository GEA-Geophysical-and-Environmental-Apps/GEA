
IOdictionary 
turbProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dictionary
FilteringProperties (turbProperties.subDict("FilteringProperties"));

dimensionedScalar alphaR("fRADIUS",dimLength, readScalar(FilteringProperties.lookup("FilteringRadius")));

 

volScalarField aU
(
  IOobject
  (
    "aU",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
 mesh	
);

volScalarField aE
(
 IOobject
 (
  "aE",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh
);


volVectorField Ufi
(
 IOobject
 (
   "U_filtered",
   runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::AUTO_WRITE
 ),
 U
);


volScalarField hefi
(
 IOobject
 (
  "he_filtered",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 thermo.he()
);

autoPtr<Ifunction> 
pIndicatorFunction 
(
 Ifunction::New(
  word(FilteringProperties.lookup("FilterType")),
  alphaR,
  aU,
  aE,
  U,
  thermo.he()
  )
);


