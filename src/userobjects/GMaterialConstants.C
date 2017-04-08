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
#include "GMaterialConstants.h"


template<>
InputParameters validParams<GMaterialConstants>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<Real>("dislocation",0.0,"dislocation density in #/um^2");
  params.addParam<Real>("i_disl_bias",1.1,"dislocation bias for intersitials");
  params.addParam<Real>("v_disl_bias",1.0,"dislocation bias for vacancies");
  params.addParam<Real>("atomic_vol",0.0,"atomic volume");
  params.addClassDescription( "Object prototype to calculate material constants");
  return params;
}

GMaterialConstants::GMaterialConstants(const InputParameters & parameters)
: GeneralUserObject(parameters),
  atomic_vol(getParam<Real>("atomic_vol")),
  _rho_d(getParam<Real>("dislocation")),
  _i_bias(getParam<Real>("i_disl_bias")),
  _v_bias(getParam<Real>("v_disl_bias"))
{
}

void GMaterialConstants::initialize()
{
}

void GMaterialConstants::execute()
{
}

void GMaterialConstants::finalize()
{
}

Real GMaterialConstants::absorb(int a,int b,std::string str1,std::string str2,double cc,int dd,int e) const{

return 0;//need overwrite

}

Real GMaterialConstants::absorbVV(int a,int b,int dd,double e) const{

return 0;//need overwrite

}

Real GMaterialConstants::absorbVI(int a,int b,int dd,double e) const{

return 0;//need overwrite

}

Real GMaterialConstants::absorbII(int a,int b,int dd,double e) const{

return 0;//need overwrite

}

Real GMaterialConstants::emit(int a,int b,double c,std::string str1,std::string str2,int cc,int d) const{

return 0;//need overwrite

}

Real GMaterialConstants::disl_ksq(int a,std::string str1,double cc,int d) const{

return 0;//need overwrite

}

Real GMaterialConstants::diff(int a,std::string str1,double cc) const{

return 0;//need overwrite

}
