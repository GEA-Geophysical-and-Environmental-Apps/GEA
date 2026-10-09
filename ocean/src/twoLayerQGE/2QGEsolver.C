
// Step 1: Compute for the potential vorticity of the top layer.

fvScalarMatrix q1Eqn
(
    fvm::ddt(q1)
  + fvm::div(phi, q1)
  + fvc::laplacian(Fr/(Re*delta), psi2)
  - fvc::laplacian(Fr/(Re*delta), psi1)
  - fvm::laplacian(1.0/Re, q1) 
 ==
    forcingTerm1
);

q1Eqn.solve();

// Step 2: Compute for the stream function of the top layer.

fvScalarMatrix psi1Eqn
(
    fvm::laplacian(dimensionedScalar("1", dimless, 1), psi1)
   + C().component(1)/Ro
  + (Fr/delta/Ro)*psi2
  - fvm::Sp(Fr/delta/Ro, psi1)
 == 
    q1/Ro
);

psi1Eqn.solve();

// Step 2.5: Update the flux phi associated to the convective term of the PDE for q1.

forAll(psiVector1, faceI)
{   
    psiVector1[faceI].component(0) = 0;
    psiVector1[faceI].component(1) = 0;
    psiVector1[faceI].component(2) = psi1[faceI];
}

curlPsi1 = fvc::curl(psiVector1);

// Forcing a null normal velocity on the boundaries.

forAll(mesh.boundaryMesh()[vIndex], faceI)
{
    curlPsi1.boundaryFieldRef()[vIndex][faceI].component(0) = 0;
}

forAll(mesh.boundaryMesh()[hIndex], faceI)
{
    curlPsi1.boundaryFieldRef()[hIndex][faceI].component(1) = 0;
}

phi = fvc::flux(curlPsi1);

// Step 3: Compute for the potential vorticity of the bottom layer.

fvScalarMatrix q2Eqn 
(
    fvm::ddt(q2)
  + fvm::div(phi2, q2)
  + fvc::laplacian(Fr/(Re*(1.0-delta)), psi1)
  - fvc::laplacian(Fr/(Re*(1.0-delta)), psi2)
  - fvm::laplacian(1.0/Re, q2)
 == 
    fvc::laplacian(-sigma, psi2)
);

q2Eqn.solve();

// Step 4: Compute for the stream function of the bottom layer.

fvScalarMatrix psi2Eqn
(
    fvm::laplacian(dimensionedScalar("1", dimless, 1), psi2)
  + C().component(1)/Ro
  + (Fr/Ro/(1.0-delta))*psi1
  - fvm::Sp(Fr/Ro/(1.0-delta), psi2)
 ==
    q2/Ro
);

psi2Eqn.solve();

// Step 4.5: Update the flux phi2 associated to the convective term of the PDE for q2.

forAll(psiVector2, faceI)
{           
    psiVector2[faceI].component(0) = 0;
    psiVector2[faceI].component(1) = 0;
    psiVector2[faceI].component(2) = psi2[faceI];
}

curlPsi2 = fvc::curl(psiVector2);

// Forcing a null normal velocity on the boundaries.

forAll(mesh.boundaryMesh()[vIndex], faceI)
{
    curlPsi2.boundaryFieldRef()[vIndex][faceI].component(0) = 0;
}

forAll(mesh.boundaryMesh()[hIndex], faceI)
{
    curlPsi2.boundaryFieldRef()[hIndex][faceI].component(1) = 0;
}

phi2 = fvc::flux(curlPsi2);


