/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "DislocationSinkRate.h"

template<>
InputParameters validParams<DislocationSinkRate>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<std::string>("Diffusivity","The diffusivity to use in this calculation");
//  params.addRequiredParam<Real>("Diffusivity","The diffusivity to use in this calculation");
  params.addParam<std::string>("DislocationDensity","","The location-dependent density of dislocations from material properties");
  params.addCoupledVar("VariedDislocation",1.0e-4,"The location-dependent density of dislocations");
  params.addParam<Real>("DislocationCoreSize", 10, "The core size of a dislocation, relative to this defect [nm]");
  return params;
}


DislocationSinkRate::DislocationSinkRate(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
   _prop_name_D(getParam<std::string>("Diffusivity")),
   _D_species(getMaterialProperty<Real>(_prop_name_D)),
//   _D_species(getParam<Real>("Diffusivity")),
   _prop_name_DD(getParam<std::string>("DislocationDensity")),
   _dislocation_density(getMaterialProperty<Real>(_prop_name_DD)),
   _varied_dislocation(coupledValue("VariedDislocation")),
   _dislocation_core_size(getParam<Real>("DislocationCoreSize"))

{}

Real
DislocationSinkRate::computeValue()
{ 
  if (!_prop_name_DD.empty() && isParamValid("VariedDislocation"))
     mooseError("Error: too many dislocation density provided" );
  if (_prop_name_DD.empty() && !(isParamValid("VariedDislocation")))
     mooseError("Error: dislocation density needed" );
  double rho_disl = (isParamValid("VariedDislocation")? _varied_dislocation[_qp] :  _dislocation_density[_qp]);
  Real Half_Dislocation_Distance = std::pow((1 / (3.14159 * rho_disl)), 0.5);

  Real K_dislocations = (2 * 3.14159 * rho_disl * _D_species[_qp])//_D_species[_qp]
                         / std::log(Half_Dislocation_Distance / _dislocation_core_size);

  return K_dislocations;
}
