#include "LinIndFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(LinIndFunction, false);
addToRunTimeSelectionTable(Ifunction, LinIndFunction, dict);

tmp<volScalarField> LinIndFunction::EvalViscosity(){
   
  // Evaluating the viscosity for the linear-step
  
  return FilteringViscosityU();
  
 }

void LinIndFunction::FilteringStep(volVectorField& Ufi,volScalarField& Efi)
{
 tmp<volScalarField> FilterViscosityStep = FilteringViscosityU();

 fvVectorMatrix UfiEqn
 (
   -fvm::laplacian(FilterViscosityStep(),Ufi) + fvm::Sp(dimensionedScalar("one",dimless,1),Ufi)
 );

 solve( UfiEqn == U_ );
 

 fvScalarMatrix EfiEqn
 (
   fvm::Sp(dimensionedScalar("one",dimless,1),Efi)-fvm::laplacian(FilterViscosityStep(),Efi)
 );
 solve( EfiEqn == E_ );


}


}// end namespace
