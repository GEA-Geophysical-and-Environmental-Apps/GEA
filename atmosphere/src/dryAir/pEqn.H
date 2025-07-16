rho = thermo.rho();

dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
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
    //+ MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

bool adjustMass = mesh.steady() && adjustPhi(phiHbyA, U, p_rgh); // aggiunto

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvc::makeRelative(phiHbyA, rho, U);

fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);

// caso non transonico
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




// aggiunta blocco if-else
//if (mesh.steady())
//{
    //#include "incompressible/continuityErrs.H"
//}
//else
//{
    p = p_rgh + rho*gh + pRef;

    //const bool constrained = fvConstraints.constrain(p); cambiato con quello sotto
    const bool constrained = true;      //questa condizione da errori fvOptions.correct(p);

    // Thermodynamic density update
    thermo.correctRho(psi*p - psip0);

    if (constrained)
    {
        rho = thermo.rho();
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
//}

// Explicitly relax pressure for momentum corrector
p_rgh.relax();

//p = p_rgh + rho*gh + pRef; //SE AMR SPENTO, LUI FA SMINCHIARE IL RISULTATO

//thermo.correctRho(psi*p - psip0); spostata nel blocco if-else

//#include "compressibleContinuityErrs.H" spostato nel blocco if-else


// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

// aggiunto
//if (mesh.steady())
//{
    fvOptions.correct(p); // anche qua non uso model ma option
//}

// aggiunto i 3 blocchi if seguenti

// For steady closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass && !thermo.incompressible()) //HO COME L'IMPRESSIONE CHE ANCHE LUI ROVINI LA SOLUZIONE
{
    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    p_rgh = p - rho*gh - pRef;
    p_rgh.correctBoundaryConditions();
}

//if (mesh.steady() || adjustMass) // modificato, prima c'era anche pimple.simpleRho()
//{
    rho = thermo.rho();
//}

//if (mesh.steady() && !pimple.transonic()) // stessa cosa di sopra   SOLUZIONE MIGLIORE SE LASCIO COMMENTATI GLI ULTIMI DUE IF
//{
    rho.relax();
//}

// aggiunto
// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);

/* cambiato in quello sotto
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
*/

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
