/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/
#include "MooseMesh.h"
#include "MaterialConstants.h"


template<>
InputParameters validParams<MaterialConstants>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription( "Object prototype to calculate material constants");
  return params;
}

MaterialConstants::MaterialConstants(const InputParameters & parameters)
: GeneralUserObject(parameters)
{
}

void MaterialConstants::initialize()
{
}

void MaterialConstants::execute()
{
}

void MaterialConstants::finalize()
{
}

Real MaterialConstants::absorb(int a,int b,std::string str1,std::string str2,double cc,int dd,int e) const{

return 0;//need overwrite

}

Real MaterialConstants::emit(int a,int b,double c,std::string str1,std::string str2,int cc,int d) const{

return 0;//need overwrite

}

Real MaterialConstants::disl_ksq(int a,std::string str1,double cc,int d) const{

return 0;//need overwrite

}

Real MaterialConstants::diff(int a,std::string str1,double cc) const{

return 0;//need overwrite

}
