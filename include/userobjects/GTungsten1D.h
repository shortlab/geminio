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

#ifndef GTungsten1D_H
#define GTungsten1D_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

class GTungsten1D : public GMaterialConstants
{
public:
  GTungsten1D(const InputParameters & parameters);

  ~GTungsten1D(){}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  Real absorbVV(int,int,int,double) const;
  Real absorbVI(int,int,int,double) const;
  Real absorbII(int,int,int,double) const;
  Real emit(int,int,double,std::string,std::string,int,int) const;
  Real disl_ksq(int,std::string,double,int=1) const;
  double energy(int,std::string,std::string) const;
  double D_prefactor(int,std::string="") const;
  double diff(int, std::string,double) const;

private:
  double Ev_formation,Ei_formation,Evb2,Eib2;//description in cpp file
  double Ei_binding_factor;
  double Ev_binding_factor;
};

template<>
InputParameters validParams<GTungsten1D>();

#endif
