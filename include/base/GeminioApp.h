#ifndef GEMINIOAPP_H
#define GEMINIOAPP_H

#include "MooseApp.h"

class GeminioApp;

template<>
InputParameters validParams<GeminioApp>();

class GeminioApp : public MooseApp
{
public:
  GeminioApp(const  InputParameters & parameters);
  virtual ~GeminioApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* GEMINIOAPP_H */
