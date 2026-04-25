#include "GradIndFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(GradIndFunction, false);
addToRunTimeSelectionTable(Ifunction, GradIndFunction, dict);


tmp<volScalarField> GradIndFunction::EvalViscosity(){
   
  // Evaluating the viscosity ( multiplyed by deltaT )
  const volVectorField& Uorg = U_;

  volScalarField::Internal SqrAbsGradU =((fvc::grad(Uorg)) && (fvc::grad(Uorg)));
  volScalarField::Internal AbsGradU = sqrt(SqrAbsGradU);
   
  dimensionedScalar AbsGradULinf = max(AbsGradU);
  dimensionedScalar one("one",AbsGradULinf.dimensions(),1.0);

  dimensionedScalar NormalizingFactor = max(one,AbsGradULinf);

  a_U.ref()  = AbsGradU/NormalizingFactor;

 // Set the boundary field to zero in order to avoid spurious oscillations // 

  forAll(a_U.boundaryField(),patchi)
  {
   a_U.boundaryFieldRef()[patchi] = scalar(0.0);
  }

 return FilteringViscosityU();

 }

void GradIndFunction::FilteringStep(volVectorField& Ufi,volScalarField& Efi)
{
  
  tmp<volScalarField> FilteringViscosityStep = EvalViscosity(); 

  // Copy for security
  const volVectorField& Uorg = U_;
  const volScalarField& Eorg = E_;

  fvVectorMatrix UfilteredEqn
  (
   -fvm::laplacian(FilteringViscosityStep(),Ufi)+fvm::Sp(dimensionedScalar("one",dimless,1),Ufi)
  );

  solve( UfilteredEqn == Uorg );
  
  fvScalarMatrix hefilteredEqn
  (
   fvm::Sp(dimensionedScalar("one",dimless,1),Efi)-fvm::laplacian(FilteringViscosityStep(),Efi)
  );

  solve ( hefilteredEqn == Eorg );
	

} //end method


}// end namespace
