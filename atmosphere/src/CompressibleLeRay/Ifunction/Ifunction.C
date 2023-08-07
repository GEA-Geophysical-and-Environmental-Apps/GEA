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
*/

// Class implementation

#include "Ifunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam{

defineTypeNameAndDebug(Ifunction,0);
defineRunTimeSelectionTable(Ifunction, dict);


autoPtr<Ifunction> Ifunction::New
(
   const word& name,
   dimensionedScalar& _alphaRadius_,
   volScalarField& _a_U,
   volScalarField& _a_E,
   volVectorField& _U_,
   volScalarField& _E_

)
{
  Info<< "Selecting indicator-function type:  " << name << endl;

  dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(name);

  if( cstrIter == dictConstructorTablePtr_->end())
  {
    FatalErrorInFunction <<
	 "Unknown filter-indicator-function"<<
	 name << nl <<nl 
	 << "Avalaiable indicator functions are" << 
	 dictConstructorTablePtr_->sortedToc() << exit(FatalError);
  }

  return autoPtr<Ifunction>
  (
   cstrIter()( _alphaRadius_ , _a_U, _a_E, _U_,_E_)
  );

}//end Selector definition

dimensionedScalar& Ifunction::alphaRadiusRef(){
	return alphaRadius_;
}
const dimensionedScalar& Ifunction::alphaRadius() const {
	return alphaRadius_;
}

tmp<volScalarField> Ifunction::FilteringViscosityU() const {
	dimensionedScalar SqRadius = pow(alphaRadius_,2);	
	tmp<volScalarField> AlphaaU(SqRadius*a_U);
	return AlphaaU;
}

tmp<volScalarField> Ifunction::FilteringViscosityE() const 
{
	dimensionedScalar SqRadius = pow(alphaRadius_,2);
	tmp<volScalarField> AlphaaT(SqRadius*a_E);
	return AlphaaT;
}

const volScalarField& Ifunction::aU() const{
	return a_U;
}

const volScalarField& Ifunction::aE() const{
	return a_E;
}

// Ifunction::~Ifunction(){};

} // end Namespace foam
// ************************************************************************* //
