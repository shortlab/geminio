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

#include "DislocationSink.h"

template<>
InputParameters validParams<DislocationSink>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("Diffusivity","The diffusivity to use in this calculation");
  params.addParam<std::string>("DislocationDensity","","The location-dependent density of dislocations from material properties");
  params.addParam<Real>("DislocationDensityValue",-1,"Direct input of dislocation density");
  params.addCoupledVar("VariedDislocation","The location-dependent density of dislocations");
  params.addParam<Real>("DislocationCoreSize", 1e-9, "The core size of a dislocation, relative to this defect");
  params.addParam<Real>("Coef", 1.0, "The efficiency factor");
  params.addCoupledVar("ConcentrationCorrection", 0, "The correction to subtract from this defect concentration, such as a thermal vacancy concentration");
  return params;
}


DislocationSink::DislocationSink(const
                                   InputParameters & parameters)
  :Kernel(parameters),
   _D_species(getParam<Real>("Diffusivity")),
   _prop_name_DD(getParam<std::string>("DislocationDensity")),
   _dislocation_density(_prop_name_DD.empty() ? NULL : &getMaterialProperty<Real>(_prop_name_DD)),
   _varied_dislocation(isCoupled("VariedDislocation")? coupledValue("VariedDislocation"):_zero),
   _rho_disl(getParam<Real>("DislocationDensityValue")),
   _dislocation_core_size(getParam<Real>("DislocationCoreSize")),
   _coef(getParam<Real>("Coef"))
{
  if (_prop_name_DD.empty() && !(isCoupled("VariedDislocation")) && _rho_disl<0)
     mooseError("Error: valid dislocation density needed" );
}

Real
DislocationSink::computeQpResidual()
{
//  return _test[_i][_qp] * _sink_rate[_qp] * _sink_concentration[_qp] * (_u[_qp] - _sink_concentration_correction[_qp]); // Positive sign because negative source from weak form PDE
  Real rho_disl = 0.0;
  if(_rho_disl>0)
      rho_disl = _rho_disl;
  else rho_disl = (isCoupled("VariedDislocation")? _varied_dislocation[_qp] :  (*_dislocation_density)[_qp]);

/* one form of dislocation sink
  Real Half_Dislocation_Distance = std::pow((1 / (3.14159 * rho_disl)), 0.5);

  Real K_dislocations = (2 * 3.14159 * rho_disl * _D_species)//_D_species[_qp]
                         / std::log(Half_Dislocation_Distance / _dislocation_core_size);
  //printf("dislocation sink rate: %f\n",K_dislocations);
*/


// one form of dislocation sink
  Real K_dislocations = rho_disl * _coef *_D_species;
  //printf("dislocation sink rate: %f\n",K_dislocations);
  return _test[_i][_qp] * K_dislocations * (_u[_qp]); // Positive sign because negative source from weak form PDE
  
}

Real
DislocationSink::computeQpJacobian()
{
//  return _test[_i][_qp] * _sink_rate[_qp] * _sink_concentration[_qp] * _phi[_j][_qp]; // Positive sign because negative source from weak form PDE
  Real rho_disl = 0.0;
  if(_rho_disl>0)
      rho_disl = _rho_disl;
  else rho_disl = (isCoupled("VariedDislocation")? _varied_dislocation[_qp] :  (*_dislocation_density)[_qp]);

/* one form of dislocation sink
  Real Half_Dislocation_Distance = std::pow((1 / (3.14159 * rho_disl)), 0.5);

  Real K_dislocations = (2 * 3.14159 * rho_disl * _D_species)//_D_species[_qp]
                         / std::log(Half_Dislocation_Distance / _dislocation_core_size);
*/

// one form of dislocation sink
  Real K_dislocations = rho_disl * _coef * _D_species;
  //printf("dislocation sink rate: %f\n",K_dislocations);
  return _test[_i][_qp] * K_dislocations * (_u[_qp]); // Positive sign because negative source from weak form PDE
  
}

Real
DislocationSink::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0;
}
