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

#ifndef GMOBILE_H
#define GMOBILE_H

#include "Kernel.h"
#include "GGroup.h"

//Forward Declarations
class GMobile;


template<>
InputParameters validParams<GMobile>();

class GMobile : public Kernel
{
public:
  
  GMobile(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  int getGroupNumber(std::string);
  double getConcBySize(int i);

private:
  int _number_v;
  int _number_i;
  int _max_mobile_v; 
  int _max_mobile_i; 
  const GGroup & _gc;
  std::vector<unsigned int> _no_v_vars;
  std::vector<const VariableValue *> _val_v_vars;
  std::vector<unsigned int> _no_i_vars;
  std::vector<const VariableValue *> _val_i_vars;
  int _cur_size;
  int max_v,max_i;
};
#endif 
