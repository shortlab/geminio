#UNITS: um,s,/um^3
#consider only vacancy cluster for tungsten
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 300    #number of vacancy variables i.e. total_groups
  number_single_v = 52  #max size with group size 1
  max_mobile_v = 8
  mobile_v_size = '1 2 3 4 5 6 7 8'

  number_i = 300      #number of interstitial variables, set to 0
  number_single_i = 52  #max size with group size 1
  max_mobile_i = 4
  mobile_i_size = '1 2 3 4'

  temperature = 723  #temperature [K]
  source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
  source_i_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'
[]

[Mesh]
  type = GeneratedMesh
# 1D ion source
  xmin = 0
  xmax = 2 #2000 change to 1, uniform source for simplicity, no spatical dependence
  dim = 1
  nx = 100
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
    data_file = spatial_defect_cluster_production.txt
    axis = 0 #x axis
    format = columns
    scale_factor = 0.25
  [../]
[]
[LotsOfSource]
  [./groups0] 
#attention to this: add 0 after defects for L0 term
    func_pre_name = func_defect_
    # either provide data file to construct functions or direct input of constant source for each size.
  [../]
[]

[AuxVariables]
  active = 'void_swelling'
  [./void_swelling]
  [../]
  [./SIA_density]
  [../]
[]

[GVoidSwelling]
  [./groups]
    aux_var = void_swelling
    group_constant = group_constant
  [../]
[]

#[GSumSIAClusterDensity]
##sum up of SIA cluster density in range [lower_bound,upper_bound]
#  [./groups]
#    aux_var = SIA_density 
#    group_constant = group_constant
#    lower_bound = 60
#  [../]
#[]

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
    dr_coef = 1.0
    update = false
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

#[Preconditioning]
#  active = smp
#  [./smp]
#    type = SMP
#    full = true
#  [../]
#[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  solve_type = 'PJFNK'
#  petsc_options =  '-snes_mf_operator'
  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value =  'hypre    boomeramg  81'
#  petsc_options_iname =  '-pc_type -sub_pc_type -ksp_gmres_restart'
#  petsc_options_value =  'bjacobi ilu  81'
  #trans_ss_check = true
  #ss_check_tol = 1.0e-14
  l_max_its =  30
  nl_max_its =  40
  nl_abs_tol=  1e-10  #Question: why change to 1e-12 not work!!!
  nl_rel_tol =  1e-7
  l_tol =  1e-5
  num_steps = 1000
  start_time = 0
  end_time = 7563.7 #15127.4  #35 dpa with 4.6e-3dpa/s
  #dt = 1.0e-2
  dtmin = 1.0e-10 
  dtmax = 50
  active = 'TimeStepper'
  [./TimeStepper]
      cutback_factor = 0.8
      dt = 1e-9
      growth_factor = 2
      type = IterationAdaptiveDT
  [../]
[]

[Debug]
#    show_top_residuals=1
#    show_var_residual_norms=1
[]


[Outputs]
  #output_linear = true
  #file_base = out
  exodus = true
  csv = true
  console = false
[]
