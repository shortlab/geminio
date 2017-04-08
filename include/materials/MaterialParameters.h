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

#ifndef MATERIALPARAMETERS_H
#define MATERIALPARAMETERS_H
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include "Moose.h"

double absorb(int,int,std::string,std::string,double,int,int);
double emit(int,int,double,std::string,std::string,int,int);
double energy(int,std::string,std::string);
double D_prefactor(int,std::string);
double diff(int, std::string,double);

#endif
