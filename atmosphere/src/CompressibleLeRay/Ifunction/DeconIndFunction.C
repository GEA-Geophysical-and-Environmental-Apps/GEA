#include "DeconIndFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(DeconIndFunction, false);
addToRunTimeSelectionTable(Ifunction, DeconIndFunction, dict);


tmp<volScalarField> DeconIndFunction::EvalViscosity(){
   
  // Evaluating the viscosity ( multiplyed by deltaT ) based on the velocity ONLY 
  
  const volVectorField& Uorg(U_);
  
  volVectorField Udummy("Udummy", Uorg ); //Needed for the linear-step. 

  Udummy.ref()*=scalar(0.0);
  Udummy.boundaryFieldRef()*=scalar(0.0);

  dimensionedScalar SqrRadius  = pow(alphaRadius(),2) ;

  fvVectorMatrix UfiEqnLinear
  (
    fvm::Sp(dimensionedScalar("one",dimless,1),Udummy)-fvm::laplacian(SqrRadius,Udummy)
  );

  solve( UfiEqnLinear == Uorg );

  /* ------- Evaluation of the indicator function ---- */ 
 
  dimensionedScalar one("one",Uorg.dimensions(),1.0);
  volVectorField::Internal deltaU =( Udummy.internalField() - Uorg.internalField() );
  volScalarField::Internal Udiff = mag(deltaU);
  dimensionedScalar maxDiff = max( Udiff );
  dimensionedScalar NormFact = max( one, maxDiff);
  
  a_U.ref() = Udiff/NormFact;

  forAll(a_U.boundaryField(),patchi)
  {
    a_U.boundaryFieldRef()[patchi] = scalar(0.0);
  }

 return FilteringViscosityU();

 }// End EvalViscosity

void DeconIndFunction::FilteringStep(volVectorField& Ufi,volScalarField& Efi)
{
  // Evaluating the viscosity thanks to a linear step in the EvalViscosity() method.
  
  const volVectorField Uorg(U_);
  const volScalarField Eorg(E_);

  tmp<volScalarField> FilteringViscosityStep = EvalViscosity(); 
  
  fvVectorMatrix UfiEqn
  (
   fvm::Sp(dimensionedScalar("one",dimless,1),Ufi) - fvm::laplacian(FilteringViscosityStep(),Ufi)
  );

  solve
  (
    UfiEqn == Uorg
  );
  
  fvScalarMatrix EfiEqn
  (
    fvm::Sp(dimensionedScalar("one",dimless,1),Efi) - fvm::laplacian(FilteringViscosityStep(),Efi)
  );
  solve
  (
   EfiEqn == Eorg
  );

} //end FilteringStep


}// end namespace
