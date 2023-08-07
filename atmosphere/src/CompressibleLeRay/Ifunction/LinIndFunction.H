// Linear Indicator Function

#ifndef LinIndFunction_H
#define LinIndFunction_H

#include "Ifunction.H"

namespace Foam
{

class LinIndFunction:
 public Ifunction
 {
  public:
   TypeName("LinIndFunction");
   
   // Default Constructor
   LinIndFunction
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


  virtual ~LinIndFunction(){};
 
  // Adding a method to filter the quantity into the class itself   

 };// end class definition

}



#endif 
