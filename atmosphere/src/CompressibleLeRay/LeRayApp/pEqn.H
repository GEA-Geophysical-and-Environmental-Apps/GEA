
rho = thermo.rho();

dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

const volScalarField psip0(psi*p);

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
);

// MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// bool adjustMass = mesh.steady() && adjustPhi(phiHbyA, U, p_rgh);

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvc::makeRelative(phiHbyA, rho, U);

fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);



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


phi = phiHbyA + p_rghEqn.flux();

p = p_rgh + rho*gh + pRef;


thermo.correctRho(psi*p - psip0);
#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"


// Explicitly relax pressure for momentum corrector
p_rgh.relax();

// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
