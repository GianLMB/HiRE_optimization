#!/bin/csh

setenv debug_status                    normal   # debug or normal or mini
setenv NDIGITS                              7   # Number of digits in the file names
setenv use_qbug                           yes   # yes or no (stops at first step)

setenv Simulation_Type                 serial   # Simulation type: serial, tempering, replica, or saxs
setenv Simulation_Method                   MD   # Simulation method: ART or MD
#setenv Replica_Type                   E_Scale  # if Simulation_Type is replica : T_Exchange or E_Scale 

#setenv chain_file              ichain_RNA.dat

setenv Titration		      .false.   # titration schee on if true
setenv n_steps_tit		         5000   # titration frequency
setenv pH                                 4.0             

setenv single_file                     .true.   # save pdb as a single file

setenv use_xtc                         .true.   # pdb or xtc ?
setenv init_singlefile                 .true.   # a single file for initiation of REMD

setenv Restart_Run                       new    # new, restart or new_velocities

setenv Temperature                   100.0d0    # Tem perature in kelvins
setenv N_Frag                              1    # Number of fragments 
setenv RANDOM_SEED                         0    # set to zero or remove to use the clock

setenv Potential_Scaling_Factor          1.0    # This is use to strengthen the potential when doing MD
#setenv Restrict_Motion               restraint # restricts motion of C_alpha - none (def), constraint,restraint 
#setenv Spring_Constant                  4.0    # Spring constant for the constrained

setenv save_unrotated_structures      .true.    # 

setenv Periodic_Boundary_Condition    .false.   # Periodic Boundary condition 
setenv Box_Length                       25.0    # The box length used for PBC
setenv PDB_center_of_mass             .false.   # writing PDB with respect to the center of mass 

setenv confining_sphere               .false.   # confining sphere
setenv radius_confining_sphere          120     # radius of confining sphere

#setenv Number_Replica                    4     # Number of replica

#setenv Temperatures_Replica      " 150 160 310 315  "  # Temperatures of replica
#setenv Exchange_every_n_steps           1000   # Try an exchange move every n steps

######    Simulated tempering ##################################################################


setenv Number_Tempering                    20   # Number of temperqtures
#setenv Temperatures_tempering   "300.00 308.18 316.57 325.20 334.06 343.16 352.52 362.12 371.99 382.13 392.54 403.24 414.22 425.51 437.11 449.02 461.26 473.82 486.74 500.00"
setenv Temperatures_tempering   "300.00 306.47 313.08 319.83 326.73 333.78 340.98 348.33 355.85 363.52 371.36 379.38 387.56 395.92 404.46 413.18 422.09 431.20 440.50 450.00"
setenv Exchange_every_n_steps             500   # Try an exchange move every n steps

#setenv Free_energy_tempering " 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 "
#setenv Last_accum_energy     " -2500.0 -2000.0 -1500.0 -1000.0 -500.0 0.0 500.0 1000.0 1500.0 2000.0 2500.0 3000.0 3500.0 4000.0 4500.0 5000.0 5500.0 6000.0 6500.0 7000.0 "
#setenv Last_accum_norm       " 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 "
#setenv Last_id                   20

 setenv Free_energy_tempering "  0.000000000000000E+000  -1.63357725235258       -3.59596270065423       -5.98636483099830       -8.88106231155363       -12.2792006173385       -16.0728009792110       -20.1244836214418       -24.3634934886054       -28.7393461815072       -33.2289551201757       -37.8006374047495       -42.4211138927762       -47.1049250461177       -51.8423177360973       -56.6110392112962       -61.4071282355901       -66.2386100109384       -71.0895163221998       -75.9571493938380      "
 setenv Last_accum_energy "   141654.884888799        164130.596353326        190463.551483191        228018.270869671        279943.262631280        358941.921571336        405588.698089831        435380.132574195        473456.876196853        504931.205859926        529211.605467981        567502.107028489        602424.563435654        637082.186810369        660599.579154563        690330.394346595        716266.120692686        741918.619080123        786048.604271161        833483.139441590      "
 setenv Last_accum_norm "   3394.13336557712        3247.63361942175        3038.13277513222        2913.51379364983        2916.16880827039        3177.90045691535        3233.96713451818        3228.83757671123        3313.23819564121        3355.95627817278        3374.37584103128        3494.50770545426        3591.67239660787        3661.52083898103        3688.30872474942        3748.38911885752        3790.01550367311        3810.51898700166        3947.66187211162        4073.56194943977      "
 setenv Last_id           19

######    SAXS-Guided simulation        #########################################################
setenv COMPUTE_SAXS_SERIAL           .false.    # Compute the SAXS force
setenv SAXS_SERIAL_STEP_MIN            1       # SAXS force every <n> steps in minimization 
setenv SAXS_SERIAL_STEP                10      # SAXS force every <n> steps in production

setenv Modulate_SAXS_Serial           .false.   # Add a gaussian periodic modulation factor for SAXS
setenv SAXS_Wave                       30.0    # Wavelength for SAXS modulation
setenv SAXS_Sigma                       3.0    # Inverse sigma for SAXS modulation
setenv SAXS_Energy_Coefficient        200.0    # Energy coefficient used in SAXS-Biased exchange step
setenv SAXS_Profiling_Rate               10    # Recompute the SAXS score every <n> exchange steps (as defined by Exchange_every_n_steps)
setenv SAXS_Norm_Type                     2    # 1 for 1-norm, 2 for 2-norm, 3 for infinity-norm
setenv SAXS_Vect_Max                    0.5    # Maximum value of the scattering vector (default and max is 1)
setenv SAXS_Norm_Range            "0.01 1.0"   # q-range on which to apply the norm calculations UNUSED
setenv SAXS_Calculation_Type          vacuum   # 'vacuum' or 'solution'
setenv SAXS_Refine_Hydration          .false.  # .true. for refined hydration, .false. otherwise (default)
setenv SAXS_Dummy_Water_Radius           3     # Radius of the dummy waters for Hydration algorithm (Only if SAXS_Refine_Hydration=.true.)
setenv SAXS_Dummy_Water_Layers           1     # Number of dummy waters shells/layers (Only if SAXS_Refine_Hydration=.true.)
setenv SAXS_w_shell                    0.037   # Excess electron density for the hydration shell
setenv SAXSPRINT   		      .true.   # Print on files SAXSs (scores) SAXSc (curves)

######  MD MIN ##########################################################################################
setenv MINIMIZATION_TYPE                MIN_damped # MIN_steep or MIN_damped
setenv Max_Num_Iter_Minimization        80000   # Maximum number of iteration (minimization) 
setenv Force_Threshold_Minimization      1.10   # Force threshold for minimization 
setenv Time_Step_DampedMD                2      # Time step for the damped MD minimization (fs)

##====Internal Coordinates
#setenv Internal_Coordinates             .true.
#setenv Backbone_Variables               .true.
#setenv Base_Sugar_Variables             .true.
#setenv Valence_Angle_Variables          .true.

######  MD  #############################################################################################

setenv timeunit                             1.0  # time unit in fs - usually 1.0
setenv timestep                             1.0  # in time units

setenv n_production                         1    # total number of production steps (time units)
#setenv simulation_time                   100    # total time of production (in ns)

setenv n_equilibration                  15000    # number of equilibration steps (time units)
setenv n_rescale_v_equil                   10    # rescale velocity every n steps during equilibration
setenv n_steps_thermalization               0    # Number of temperatures used in thermalization to target
                                                 # temperature (multiple of 1.5)

setenv thermostat                      langevin  # Thermostat (none, berendsen or langevin)
setenv friction_coef                     0.01    # Langevin friction coefficient (frequency)
#setenv berendsen_time                    500    # Berendsen time (time units)

setenv n_correct_center_of_mass          1000    # Correct center of mass every n steps   
setenv n_correct_rotation                1000    # Correct rotation every n steps   

setenv n_statistics                    1000      # Accumulate statistics every n steps
setenv n_stats_average                100000     # Average statistics every n steps

setenv n_save_configs                  1000        # Save configs every n steps (default: 1000)


setenv test_on_LJ                     .false.    # default is false, only used for debug
setenv force_calc_test                .false.    # default is false, only used for debug
setenv chk_ene                        .false.    # default is false, only used for debug

################### Titration ###########################
#sed "s/pH/${pH}/" titration_ff.in > cyl.json

############### Run the simulation ######################################################################

#./batch.sh

#cp cutoffs.dat1 cutoffs.dat
#cp cutoffs.dat-new cutoffs.dat


#system specific parameters
source parametres.csh


# call simulator code
./simulator_test
