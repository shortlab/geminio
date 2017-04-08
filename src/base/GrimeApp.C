#include "GrimeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// Include kernel header files
#include "DefectSink.h"
#include "DefectSource.h"
#include "DefectRecombination.h"
#include "MatPropDiffusion.h"
#include "FuncCoefVariable.h"
//for cluster
#include "VariableProduct.h"
#include "CoeffDiffusion.h"
#include "SingleVariable.h"
#include "DislocationSink.h"
#include "UserObjectVariableProduct.h"
#include "UserObjectSingleVariable.h"
#include "UserObjectDiffusion.h"
#include "MobileDefects.h"
#include "ImmobileDefects.h"

// Include auxkernel header files
#include "DislocationSinkRate.h"
#include "DefectRecombinationRateConstant.h"
#include "VoidSinkRate.h"
#include "GVoidSwelling.h"
#include "GSumSIAClusterDensity.h"
#include "ClusterDensity.h"


// Next come the materials
#include "RadiationMaterial.h"
#include "ArcMaterial.h"

// ********************actions********************* //
#include "AddLotsOfVariableProduct.h"
#include "AddLotsOfTimeDerivative.h"
#include "AddLotsOfCoeffDiffusion.h"
#include "AddLotsOfSingleVariable.h"
#include "AddLotsOfSource.h"
#include "AddClusterICAction.h"
#include "AddLotsOfSink_disl.h"
#include "AddLotsOfVariableAction.h"
#include "AddUserObjectVariableProduct.h"
#include "AddUserObjectSingleVariable.h"
#include "AddUserObjectDiffusion.h"
#include "AddUserObjectDislocationSink.h"
#include "AddMobileDefects.h"
#include "AddImmobileDefects.h"
#include "AddClusterDensity.h"
#include "AddGVoidSwelling.h"
#include "AddGSumSIAClusterDensity.h"

//*************Postprocessors**********************//
#include "NodalConservationCheck.h"
#include "TotalDefectLoss.h"

//*************UserObjects**************************//
#include "GroupConstant.h"
#include "MaterialConstants.h"
#include "TestProperty.h"
#include "GroupingTest.h"
#include "BCCIronProperty.h"




/***************grouping method start*********************/
//###################kernels################//
#include "GImmobileL0.h"
#include "GImmobileL1.h"
#include "GMobile.h"
#include "ConstantKernel.h"
//#####################Actions##############//
#include "AddGVariable.h"
#include "AddGImmobile.h"
#include "AddGMobile.h"
#include "AddGTimeDerivative.h"
#include "AddGConstantKernels.h"
//#####################user objects##########//
#include "GGroup.h"
#include "GMaterialConstants.h"
#include "GGroupingTest.h"
#include "GIron.h"
#include "GTungsten.h"
/***************grouping method end*********************/


template<>
InputParameters validParams<GrimeApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

GrimeApp::GrimeApp(const InputParameters & parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  GrimeApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  GrimeApp::associateSyntax(_syntax, _action_factory);
}

GrimeApp::~GrimeApp()
{
}

void
GrimeApp::registerApps()
{
  registerApp(GrimeApp);
}

void
GrimeApp::registerObjects(Factory & factory)
{
  // Register kernels
  registerKernel(DefectSink);
  registerKernel(DefectSource);
  registerKernel(DefectRecombination);
  registerKernel(MatPropDiffusion);
  registerKernel(FuncCoefVariable);
  registerKernel(VariableProduct);
  registerKernel(CoeffDiffusion);
  registerKernel(SingleVariable);
  registerKernel(UserObjectSingleVariable);
  registerKernel(UserObjectDiffusion);
  registerKernel(DislocationSink);
  registerKernel(UserObjectVariableProduct);
  registerKernel(MobileDefects);
  registerKernel(ImmobileDefects);

  // Register auxkernels
  registerAux(DislocationSinkRate);
  registerAux(DefectRecombinationRateConstant);
  registerAux(VoidSinkRate);
  registerAux(GVoidSwelling);
  registerAux(ClusterDensity);
  registerAux(GSumSIAClusterDensity);


  // Register materials classes
  registerMaterial(RadiationMaterial);
  registerMaterial(ArcMaterial);

  //register postprocessors
  registerPostprocessor(NodalConservationCheck);
  registerPostprocessor(TotalDefectLoss);

  // Register UserObjects
  registerUserObject(GroupConstant);
  registerUserObject(MaterialConstants);

  registerUserObject(TestProperty);
  registerUserObject(GroupingTest);
  registerUserObject(BCCIronProperty);


//grouping method
  //register kernels
  registerKernel(GMobile);
  registerKernel(GImmobileL0);
  registerKernel(GImmobileL1);
  registerKernel(ConstantKernel);
  //register userobjects
  registerUserObject(GGroup);
  registerUserObject(GMaterialConstants);
  registerUserObject(GGroupingTest);
  registerUserObject(GIron);
  registerUserObject(GTungsten);


}

void
GrimeApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
//actions
  registerAction(AddLotsOfVariableProduct, "add_kernel");
  registerAction(AddLotsOfTimeDerivative,"add_kernel");
  registerAction(AddLotsOfCoeffDiffusion,"add_kernel");
  registerAction(AddLotsOfSingleVariable,"add_kernel");
  registerAction(AddUserObjectSingleVariable,"add_kernel");
  registerAction(AddUserObjectDiffusion,"add_kernel");
  registerAction(AddUserObjectDislocationSink,"add_kernel");
  registerAction(AddLotsOfSource,"add_kernel");
  registerAction(AddUserObjectVariableProduct, "add_kernel");
  registerAction(AddClusterICAction, "add_ic");
  registerAction(AddLotsOfSink_disl,"add_kernel");
  registerAction(AddLotsOfVariableAction,"add_variable");
  registerAction(AddLotsOfVariableAction,"add_ic");
  registerAction(AddLotsOfVariableAction,"add_bc");
  registerAction(AddMobileDefects,"add_kernel");
  registerAction(AddImmobileDefects,"add_kernel");
  registerAction(AddClusterDensity,"add_aux_kernel");

//syntax
  syntax.registerActionSyntax("AddLotsOfVariableProduct", "LotsOfVariableProduct/*");
  syntax.registerActionSyntax("AddUserObjectVariableProduct", "LotsOfUserObjectVariableProduct/*");
  syntax.registerActionSyntax("AddLotsOfTimeDerivative", "LotsOfTimeDerivative/*");
  syntax.registerActionSyntax("AddLotsOfCoeffDiffusion", "LotsOfCoeffDiffusion/*");
  syntax.registerActionSyntax("AddLotsOfSingleVariable", "LotsOfSingleVariable/*");
  syntax.registerActionSyntax("AddUserObjectSingleVariable", "LotsOfUserObjectSingleVariable/*");
  syntax.registerActionSyntax("AddUserObjectDiffusion", "LotsOfUserObjectDiffusion/*");
  syntax.registerActionSyntax("AddUserObjectDislocationSink", "LotsOfUserObjectDislocationSink/*");
  syntax.registerActionSyntax("AddLotsOfSource", "LotsOfSource/*");
  syntax.registerActionSyntax("AddClusterICAction", "ClusterIC/*");
  syntax.registerActionSyntax("AddLotsOfSink_disl", "LotsOfSink_disl/*");
  syntax.registerActionSyntax("AddLotsOfVariableAction", "LotsOfVariables/*");
  syntax.registerActionSyntax("AddMobileDefects", "MobileDefects/*");
  syntax.registerActionSyntax("AddImmobileDefects", "ImmobileDefects/*");
  syntax.registerActionSyntax("AddClusterDensity", "ClusterDensity/*");

//***grouping method
//actions
  registerAction(AddGVariable,"add_variable");
  registerAction(AddGVariable,"add_ic");
  registerAction(AddGVariable,"add_bc");
  registerAction(AddGMobile,"add_kernel");
  registerAction(AddGImmobile,"add_kernel");
  registerAction(AddGTimeDerivative,"add_kernel");
  registerAction(AddGConstantKernels,"add_kernel");
  registerAction(AddGVoidSwelling,"add_aux_kernel");
  registerAction(AddGSumSIAClusterDensity,"add_aux_kernel");
//syntax
  syntax.registerActionSyntax("AddGVariable","GVariable/*");
  syntax.registerActionSyntax("AddGMobile","GMobile/*");
  syntax.registerActionSyntax("AddGImmobile","GImmobile/*");
  syntax.registerActionSyntax("AddGTimeDerivative", "GTimeDerivative/*");
  syntax.registerActionSyntax("AddGConstantKernels", "Sources/*");
  syntax.registerActionSyntax("AddGVoidSwelling", "GVoidSwelling/*");
  syntax.registerActionSyntax("AddGSumSIAClusterDensity", "GSumSIAClusterDensity/*");

}
