//if (!mesh.steady() && !pimple.simpleRho())
{
   // rho = thermo.rho();
      //rho = p*psi;
     rho = (pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp());
       //rho = p*psi;
}
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
//const volScalarField psip0(psi*p);
const volScalarField psip0((pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp()));


const volScalarField rAU("rAU", 1.0/UEqn.A());
const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
 // + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

bool adjustMass = mesh.steady() && adjustPhi(phiHbyA, U, p_rgh);

//surfaceScalarField phig(rhorAUf*fvc::interpolate(rho - rho_hyd)*(g & mesh.Sf()));

//phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvc::makeRelative(phiHbyA, rho, U);

fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p_rgh)*phiHbyA/fvc::interpolate(rho);

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA) + fvm::div(phid, p_rgh)
     //==
    //fvOptions(psi, p_rgh, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        p_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAUf, p_rgh);

        // Relax the pressure equation to ensure diagonal-dominance
        p_rghEqn.relax();

        p_rghEqn.setReference
    (
        pRefCell,
        compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
    );
    
    

        p_rghEqn.solve();
    }
}
else
{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
     ==
        fvOptions(psi, p_rgh, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        p_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAUf, p_rgh);

        p_rghEqn.setReference
    (
        pRefCell,
        compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
    );

        p_rghEqn.solve();
    }
}

phi = phiHbyA + p_rghEqn.flux();

//p = p_rgh + rho*gh + pRef;
p = p_rgh + p0 + pRef;
if (mesh.steady())
{

    //#include "incompressible/continuityErrs.H"
}
else
{
    //const bool constrained = fvOptions.correct(p); //fvConstraints.constrain(p);

    // Thermodynamic density update
    thermo.correctRho((pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp()) - psip0);
    //thermo.correctRho((pRef/(thermo.Cp()-thermo.Cv())/theta)*pow(p/pRef, thermo.Cv()/thermo.Cp()) - psip0);
    Info << "provaaaa";
    //if (constrained)
    //{
    //    rho = thermo.rho();
    //}

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}

// Explicitly relax pressure for momentum corrector
p_rgh.relax();

// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
//U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux())/rhorAUf);
U.correctBoundaryConditions();
        fvOptions.correct(U);
K = 0.5*magSqr(U);

/*if (mesh.steady())
{
    fvOptions.correct(p);
}*/

// For steady closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass && !thermo.incompressible())
{
    /*p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    //p_rgh = p - rho*gh - pRef;
    p_rgh = p - p0 - pRef;
    p.correctBoundaryConditions();*/
    cout << "diagnostica";
}

/*if (mesh.steady() || pimple.simpleRho() || adjustMass)
{
    rho = thermo.rho();
}*/

/*if ((mesh.steady() || pimple.simpleRho()) && !pimple.transonic())
{
    rho.relax();
}*/

// Correct rhoUf if the mesh is moving
//fvc::correctRhoUf(rhoUf, rho, U, phi);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p_rgh);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
