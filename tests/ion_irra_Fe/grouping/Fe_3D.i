# UNITS: um,s,#/um^3
# METHODS: grouping

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.
  number_v = 120    #number of vacancy variables i.e. total groups
  number_single_v = 20  #max size with group size 1
  max_mobile_v = 1

  number_i = 120      #number of interstitial variables i.e. total groups
  number_single_i = 20  #max size with group size 1
  max_mobile_i = 3

  temperature = 723  #temperature [K]
  source_v_size = '1 2 3 4 5 6 7 8 9 10'
  source_i_size = '1 2 3 4 5 6 7 8 9 10'
[]

[Mesh]
  # 1D ion source
  type = GeneratedMesh
  xmin = 0
  xmax = 1.6 #um with spatical dependence
  dim = 1
  nx = 50
[]

# define defect variables, set variables and boundadry condition as 0 where appropriate
[GVariable]
  [./groups]
    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = dirichlet
    #IC_v_size = '1 2 3 4'
    #IC_v = '2000.0 4000.0 1000.0 250.0' #'3.9            2.323' #thermal equil
    IC_v_size = ''
    IC_v = '' #'3.9            2.323' #thermal equil
    #initial concentration for species with value NON-ZERO
    IC_i_size = ''
    IC_i = '' #thermal equil
  [../]
[]

[GTimeDerivative]
  [./groups]
  [../]
[]

[GMobile]
  [./groups]
    group_constant = group_constant
  [../]
[]

[GImmobile]
  [./groups]
    group_constant = group_constant
  [../]
[]

[LotsOfFunction]
  [./func_defect_] 
#for now only 1D, and data file should be in column, postion+v_size+i_size; defectsi1 should include injected interstitials. Make sure unit is in um^3 or um
#data should be the same with [sources] block and the sub_block name [func_defect_]
    type = PiecewiseLinearTimeLimit
    tlimit = 2.0e8 #limit the time with sources [s]
    data_file = spatial_defect_cluster_production_powerlaw_v10i10.txt
    axis = 0 #x axis
    format = columns
    scale_factor = 1.0e-2 #scaling factor in 
  [../]
[]
[LotsOfSource]
  #Attention: add 0 (e.g. groups0) after defects for L0 term
  [./groups0] 
    func_pre_name = func_defect_
    # either provide data file to construct functions or direct input of constant source for each size.
  [../]
[]

[LotsOfUserObjectDiffusion]
  #Attention: add 0 (e.g. groups0) after defects for L0 term
  [./groups0]
    group_constant = group_constant
  [../]
[]

[AuxVariables]
  active = 'void_swelling'
  [./void_swelling]
  [../]
[]

[GVoidSwelling]
  [./groups]
    aux_var = void_swelling
    group_constant = group_constant
    lower_bound = 360 #lower bound cluster size, 360->2.0nm in diameter
  [../]
[]

[UserObjects]
  [./material]
    type = GIron   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 100 #dislocation density 1.0 /um^2
  [../]

  [./group_constant]
    type = GGroup
    material = 'material'
    #GroupScheme = Uniform
    GroupScheme = RSpace
    dr_coef = 0.5
    update = false #updating grouping scheme during solving process, set to false only
    execute_on = initial
  [../]
[]

[Postprocessors]
  [./FluxChecker-V]
    type = NodalVariableValue
    nodeid = 1
    variable = groups0v1
  [../]
  [./FluxChecker-I]
    type = NodalVariableValue
    nodeid = 1
    variable = groups0i1
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value =  'hypre    boomeramg  81'
  l_max_its =  30
  nl_max_its =  40
  nl_abs_tol=  1e-10  #Question: why change to 1e-12 not work!!!
  nl_rel_tol =  1e-7
  l_tol =  1e-5
  num_steps = 1000
  start_time = 0
  end_time = 1.0e2  #certain dpa with dose rate defined by source term: 4.6e-3dpa/s in this case
  dtmin = 1.0e-10 
  dtmax = 50
  active = 'TimeStepper'
  [./TimeStepper]
      cutback_factor = 0.4
      dt = 1e-9
      growth_factor = 2
      type = IterationAdaptiveDT
  [../]
[]

[Outputs]
  #file_base = out
  exodus = true
  csv = true
  console = true
[]
