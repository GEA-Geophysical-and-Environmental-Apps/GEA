/*---------------------------------------------------------------------------*\
License
    This file is part of GEA
    GEA is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GEA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
    

Description
    
    General class useful to filter the compressible Navier Stokes equations:
     
     alpha_Radius_ : Filtering Radius 
     a_U           : Indicator Function for the velocity
     a_E           : Indicator Function for the Energy

    
 
\*---------------------------------------------------------------------------*/

#ifndef Ifunction_H
#define Ifunction_H

#include "fvCFD.H"


namespace Foam
{

class Ifunction
{
 protected:

  dimensionedScalar& alphaRadius_; //FilteringRadius
  volScalarField& a_U; //Indicator function for the velocity
  volScalarField& a_E; //Indicator Function for the energy
  volVectorField& U_; // Velocity
  volScalarField& E_; // Energy 

  public:
  
  TypeName("Ifunction");
  
  declareRunTimeSelectionTable
  (
   autoPtr,
   Ifunction,
   dict,
   (
    dimensionedScalar& _alphaRadius_,
    volScalarField& _a_U,
    volScalarField& _a_E,
    volVectorField& _U_,
    volScalarField& _E_
   ),
   (_alphaRadius_,_a_U,_a_E,_U_,_E_)
  );

  // Constructor:

  Ifunction
   (
    dimensionedScalar& _alphaRadius_,
    volScalarField& _aU,
    volScalarField& _aE,
    volVectorField& _U,
    volScalarField& _E
    ):
    alphaRadius_(_alphaRadius_),
    a_U(_aU),
    a_E(_aE),
    U_(_U),
    E_(_E)
  {};
      
    //Selectors
    static autoPtr<Ifunction> New
    (
      const word& name,
      dimensionedScalar& _alphaRadius_,
      volScalarField& _a_U,
      volScalarField& _a_E,
      volVectorField& _U_,
      volScalarField& _E_
    );

  Ifunction(const Ifunction&) = delete;

  dimensionedScalar&  alphaRadiusRef(); // Filtering radius 

  const dimensionedScalar& alphaRadius() const;
  
   // Filtering-viscosities:

   
   tmp<volScalarField> FilteringViscosityU() const;
  
   tmp<volScalarField> FilteringViscosityE() const;

  const volScalarField& aU() const;

  const volScalarField& aE() const;

  virtual tmp<volScalarField>  EvalViscosity() = 0; 
   
  virtual void FilteringStep( volVectorField& Ufi,volScalarField& Efi ) = 0;
   
    
  virtual ~Ifunction(){};
  

 };
}

#endif

// ************************************************************************* //
