// Linear Indicator Function

#ifndef DeconIndFunction_H
#define DeconIndFunction_H

#include "Ifunction.H"

namespace Foam
{

class DeconIndFunction:
 public Ifunction
 {
  public:
   TypeName("DeconIndFunction");
   
   // Default Constructor
   DeconIndFunction
   (
    dimensionedScalar& _alphaRadius_,
    volScalarField& _aU,
    volScalarField& _aE,
    volVectorField& _U,
    volScalarField& _E
   ):
   Ifunction
   (
    _alphaRadius_,
    _aU,
    _aE,
    _U,
    _E
   )
  {};
  
  virtual tmp<volScalarField> EvalViscosity(); 
  
  virtual void FilteringStep( volVectorField& Ufi,volScalarField& Efi );

  virtual ~DeconIndFunction(){};
   
 };// end class definition

}// end namespace



#endif 
