Info<< "Calculating face flux field for the curls of psi1 and psi2\n" << endl;

// Creating a vectors of the form (0, 0, psi) for the computation of the curl
volVectorField psiVector1(k*psi1);
volVectorField psiVector2(k*psi2);

volVectorField curlPsi1(fvc::curl(psiVector1));
volVectorField curlPsi2(fvc::curl(psiVector2));

// finding patch IDs
label vIndex = mesh.boundaryMesh().findPatchID("verticalWalls");
label hIndex = mesh.boundaryMesh().findPatchID("horizontalWalls");

// enforcing null normal velocity on the boundary patches
forAll(mesh.boundaryMesh()[vIndex], faceI)
{
    curlPsi1.boundaryFieldRef()[vIndex][faceI].component(0) = 0;
    curlPsi2.boundaryFieldRef()[vIndex][faceI].component(0) = 0;
}

forAll(mesh.boundaryMesh()[hIndex], faceI)
{
    curlPsi1.boundaryFieldRef()[hIndex][faceI].component(1) = 0;
    curlPsi2.boundaryFieldRef()[hIndex][faceI].component(1) = 0;
}

surfaceScalarField phi(fvc::flux(curlPsi1));
surfaceScalarField phi2(fvc::flux(curlPsi2));


// ************************************************************************* //
