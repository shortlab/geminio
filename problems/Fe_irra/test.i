#UNITS: um,s,/um^3
#solve discrete cluster dynamics equation and compare with results in T.Jourdan Efficient simulation of kinetics of radiation induced defects (figure 5)
#no grain boundary sink is considered
[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.
  number_v = 1000
  number_i = 5000
  mobile_v_size = '1' #'1 2 3 4'
  mobile_i_size = '1' #'1 2 3'
  temperature = 600  #temperature [K]
  source_v_size = '1 2 3 4'
  source_i_size = '1 2 3'
[]

[Mesh]
  type = GeneratedMesh
# 1D ion source
  xmin = 0
  xmax = 2 #2000 change to 1, uniform source for simplicity, no spatical dependence
  nx = 50
  dim = 1
[]


[LotsOfVariables]
# define defect variables, set variables and boundadry condition as 0
  [./defects]
    scaling = 1.0  #important factor, crucial to converge
    boundary_value = 0.0
    bc_type = dirichlet
  [../]
[]

[LotsOfFunction]
#checked: no problem in adding many functions here
  [./func_defect_] 
#for now only 1D, and data file should be in column, postion+v_size+i_size; defectsi1 should include injected interstitials. Make sure unit is in nm or nm^3
#data should be the same with [sources] block and the sub_block name [func_defect_]
    type = PiecewiseLinearTimeLimit
    tlimit = 2.0e8 #limit the time with sources [s]
    data_file = new_data.txt
    axis = 0 #x axis
    format = columns
    scale_factor = 1.0
  [../]
[]

[ClusterIC]
  [./defects]
    IC_v_size = '1'
    IC_i_size = ''
    #initial concentration for species with value NON-ZERO
    IC_v = '0.0'        #'3.9            2.323'
    IC_i = ''    #'2.79            2.34        0.57' #corresponds to IC_i_size
  [../]
[]

[LotsOfVariableProduct]
  [./defects]
    custom_input = false
    #set this flag to be true to accept mannully input values for coefficients
  [../]
[]

[LotsOfSingleVariable]
  [./defects]
    custom_input = false 
    #set true if mannully input coefficients
  [../]
[]

[LotsOfTimeDerivative]
  [./defects]
  [../]
[]

[LotsOfSource]
  [./defects]
    func_pre_name = func_defect_
    # either provide data file to construct functions or direct input of constant source for each size.
  [../]
[]

[LotsOfSink_disl]
#essentially variable product. separate this for clearness
  [./defects]
#selectively add kernels by providing mobile species size, consistent with diffusion block, but not a good practice if there are too many mobile species
    custom_input = false
    #dislocation = dislocation
    dislocation_density = 1.0 
    i_disl_bias = 1.1
    v_disl_bias = 1.0
    #Rvd = 1.2e-3 #um capture radius of vacancies by dislocations
    #Rid = 3.6e-3 #um capture radius of interstitials by dislocations
  [../]
[]

[Postprocessors]
  
  [./FluxChecker-V]
    type = NodalVariableValue
    nodeid = 1
    variable = defectsv1
  [../]
#  [./Total_Number_v]
#    type = NodalConservationCheck
#    nodeid = 1
#    var_prefix = 'defectsv'
#    size_range = '1 200'
#  [../]
#  [./Total_Number_i]
#    type = NodalConservationCheck
#    nodeid = 1
#    var_prefix = 'defectsi'
#    size_range = '1 1000'
#  [../]
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
  scheme = bdf2
#  solve_type = NEWTON
  solve_type = 'PJFNK'
#  petsc_options =  '-snes_mf_operator'
  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value =  'hypre    boomeramg  31'
#  petsc_options_iname =  '-pc_type -sub_pc_type -ksp_gmres_restart'
#  petsc_options_value =  'bjacobi ilu  81'
  l_max_its =  30
  nl_max_its =  40
  nl_abs_tol=  1.0e-10
  nl_rel_tol =  1.0e-6
  l_tol =  1.0e-6
  num_steps = 400
  start_time = 0
  end_time = 1.0e5
 # dt = 1e0
  dtmin = 1.0e-10 
  #dtmax = 1000
  #dt = 1.0e-2
  active = 'TimeStepper'

  [./TimeStepper]
      cutback_factor = 0.4
      dt = 1e-9
      growth_factor = 2
      type = IterationAdaptiveDT
  [../]
  [./TimeStepper0]
      type                     = DT2
      dt                       = 1e-9                        # The initial time step size.
      e_max                    = 1e7                  # Maximum acceptable error.
      e_tol                    = 1e6                  # Target error tolerance.
      max_increase             = 1.5                       # Maximum ratio that the time step can increase.
  [../]
[]

[Outputs]
  #output_linear = true
#  file_base = out
  exodus = true
  csv = true
  console = true
  print_perf_log = true
[]
