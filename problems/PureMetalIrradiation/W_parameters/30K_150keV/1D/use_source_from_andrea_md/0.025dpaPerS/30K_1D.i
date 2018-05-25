#UNITS: um,s,/um^3
# implement grouping method

[GlobalParams]
#set the largest size for vacancy clusters and interstitial clusters. Also defined in blocks to be clearer.

  number_v = 150    #number of vacancy variables i.e. total_groups
  number_single_v = 101  #max size with group size 1
  max_mobile_v = 1

  number_i = 250      #number of interstitial variables, set to 0
  number_single_i = 101  #max size with group size 1
  max_mobile_i = 5

  temperature = 30  #temperature [K]
  SIAMotionDim = 1D
  aux_prefix = rLambda
[]

[Mesh]
  type = GeneratedMesh
  xmin = 0
  xmax = 1 #uniform source for simplicity, no spatical dependence
  dim = 1
  nx = 2
[]

# define defect variables, set variables and boundadry condition as 0 where appropriate
[GVariable]
  [./groups]
#    boundary_value = 0.0
    scaling = 1.0  #important factor, crucial to converge
    bc_type = neumann
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

[Sources]
  [./groups]
    source_i_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100'
    source_i_value = '201650411.47 38205608.66 14438035.36 7238609.25 4237123.3 2735496.17 1889573.62 1371459.99 1033753.73 1231122.26 1037037.24 886696.39 767720.52 671847.52 593385.48 528305.6 473688.81 427376.54 387743.66 353547.03 323822.08 297810.74 274910.34 254636.67 236596.96 220469.83 205990.13 192937.54 181127.66 170405.23 160638.68 151715.94 143541.04 136031.39 129115.6 122731.68 116825.62 111350.13 106263.69 101529.72 97115.88 92993.47 89136.99 85523.68 82133.18 78947.26 75949.52 73125.18 70460.91 67944.64 65565.42 63313.33 61179.32 59155.13 57233.24 55406.76 53669.36 52015.26 50439.13 48936.05 47501.52 46131.35 44821.69 43568.96 42369.87 41221.34 40120.52 39064.76 38051.59 37078.72 36144 35245.42 34381.12 33549.35 32748.46 31976.93 31233.3 30516.23 29824.45 29156.76 28512.03 27889.21 27287.3 26705.36 26142.5 25597.88 25070.7 24560.23 24065.74 23586.56 23122.07 22671.65 22234.73 21810.77 21399.26 20999.69 20611.62 20234.58 19868.17 19511.97'
    source_v_size = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100'
    source_v_value = '420846869.43 52605858.68 15586921.09 6575732.33 3366774.96 1948365.14 1226958.8 821966.54 577293.37 416431.48 354140.82 305446.7 266587.55 235031.4 209020.27 187300.78 168958.63 153313.46 139849.99 128171.82 117969.78 108999.62 101066.27 94012.25 87709.28 82051.94 76952.93 72339.4 68150.17 64333.54 60845.53 57648.58 54710.4 52003.12 49502.59 47187.75 45040.23 43043.88 41184.49 39449.51 37827.79 36309.46 34885.68 33548.57 32291.06 31106.79 29990.05 28935.66 27938.96 26995.7 26102.03 25254.45 24449.76 23685.04 22957.63 22265.06 21605.1 20975.68 20374.88 19800.97 19252.31 18727.41 18224.88 17743.43 17281.87 16839.1 16414.07 16005.84 15613.49 15236.21 14873.2 14523.73 14187.13 13862.76 13550 13248.31 12957.15 12676.01 12404.45 12142.01 11888.28 11642.87 11405.41 11175.55 10952.96 10737.33 10528.36 10325.78 10129.33 9938.74 9753.79 9574.24 9399.89 9230.52 9065.95 8906 8750.47 8599.22 8452.08 8308.9'
    scaling_factor = 2.0 
  [../]
[]

[RecipMeanFreePath]
  [./groups]
    group_constant = group_constant
  [../]
[]

[UserObjects]
  [./material]
    type = GTungsten1D   #definition should be in front of the usage
    i_disl_bias = 1.15
    v_disl_bias = 1.0
    dislocation = 1 #dislocation density 1.0 /um^2
  [../]

  [./group_constant]
    type = GGroup
    material = 'material'
    GroupScheme = RSpace
    dr_coef = 0.5
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
#  petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value =  'hypre    boomeramg  81'
  petsc_options_iname =  '-pc_type -sub_pc_type -ksp_gmres_restart'
  petsc_options_value =  'bjacobi ilu  81'
  #trans_ss_check = true
  #ss_check_tol = 1.0e-14
  l_max_its =  30
  nl_max_its =  40
  nl_abs_tol=  1e-10  #Question: why change to 1e-12 not work!!!
  nl_rel_tol =  1e-7
  l_tol =  1e-8
  num_steps = 500
  start_time = 0
  end_time = 0.558
  #dt = 1.0e-2
  dtmin = 1.0e-10 
  dtmax = 0.005
  active = 'TimeStepper'
  [./TimeStepper]
      cutback_factor = 0.4
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
