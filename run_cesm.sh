# This is a testing bash shell script for running CESM
# component set: FHS94, T85L60 [PK02][MODEL_TOP=2.2e-2 hPa][nlev=60]
# "*" for basic steps
# testing case name: HS94_test

# [PK02] for Stratosphere configuration in Polvani & Kushner (2002)
# https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2001GL014284
# https://www.cesm.ucar.edu/models/simple/change-trefana

######################
# alteration tracing # ==================================================================================================================================================
######################                                                               |                              location                                 |          |
#                                            alteration                              |                        (Directory/Files)                              |   Stage  |                                             |
# =======================================================================================================================================================================
# <A1> Revalue environmental variables in XML files                                  | $CASEROOT                                                             | Step 1.5 |
# <A2> Setting parallel computing configuration (multi-processing multi-threading)   | $CASEROOT                                                             | Step 1.5 |
# <A3> Namelist modification - specify output format.                                | $CASEROOT/user_nl_cam                                                 | Step 2.5 |
# <A4> Namelist modification - specify output variables.                             | $CASEROOT/user_nl_cam                                                 | Step 2.5 |
# <PK1> [PK02] Namelist modification -                                               | $CASEROOT/user_nl_cam                                                 | Step 2.5 |
#              specify namelist parameters to control the use of the                 |                                                                       |          |
#              PK02 stratosphere                                                     |                                                                       |          |
# <PK2> [PK02] Source Code Modification for PK02 Stratosphere configuration          | $CASEROOT/SourceMods/src.cam                                          | Step 2.5 |
# <PK3> [PK02] Adding Namelist definitions for PK02 source code modifications        | $SRCROOT/components/cam/bld/namelist_files/namelist_definition.xml    | Step 2.5 |
# <PK4> [PK02] Regenerate appropriate input data that fit the PK02 configuration.    | $DIN_LOC_ROOT/cam/inic/dabiic/makeic.ncl                              | Step 2.5 |
#              (using uneven grid spacing in vertical direction)                     | $DIN_LOC_ROOT/cam/inic/dabiic/levels.ncl                              | Step 2.5 |
#                                                                                    | $DIN_LOC_ROOT/cam/inic/dabiic/ncdata/HS1994.128x256.L60.nc            | Step 2.5 |
# <PK5> [PK02] Namelist modification - revalue 'ncdata' to the location of           | $CASEROOT/user_nl_cam                                                 | Step 2.5 |
#              expected (new) initial condition file                                 |                                                                       |          |
# <PK6> [PK02] Namelist modification -                                               | $CASEROOT/user_nl_cam                                                 | Step 2.5 |
#              revalue hyperparameters for modified resolution configuration         |                                                                       |          |
# =======================================================================================================================================================================


###################################
# declare environmental variables #
###################################
export CESMDIR=/home/ur12229009/cesm                       # root directory for all CESM releases/tags
export SRCDIR=$CESMDIR/code/cesm2_3_alpha17b               # root directory for CESM code used ($SRCROOT in CESM world)
export DATADIR=/work/ur12229009/cesm_inputdata             # root directory for CESM input data

export CASENAME=HS94_test                                  # Case Name
export CASEDIR=$CESMDIR/cases/$CASENAME                    # Case Directory ($CASEROOT in CESM world)
export COMPSET=FHS94                                       # Component Sets
export GRID=T85z60_T85_mg17                                # Resolution Settings

export MACHINE_NAME=t3                                     # machine where CESM code run
export PROJECT_ID=MST112215                                 # PROJECT value required on $MACHINE
export WHICH_QUEUE=ctest                                   # the machine queue in which to submit the job.


###############################
# * Step 1  - create new case # ========================================================================================
###############################
# After "create_newcase" is issued, the associated directories will be build:
# Case Directory (path: $CASEDIR) - contain 3 scripts needed for the following steps:
#     1. case.setup   (scripts for set up the directories for the following steps)
#     2. case.build   (scripts for compiling the codes)
#     3. case.submit. (script for submitting the case on the machine)
# ======================================================================================================================

# go into scripts directory into the CESM code
cd $SRCDIR/cime/scripts

# create a new case in the directory “cases” in your home directory
./create_newcase --case $CASEDIR --compset $COMPSET --res $GRID --machine $MACHINE_NAME --mpilib intelmpi

# go into the case you just created in the last step
cd $CASEDIR


############################
# Step 1.5 - modified .XML # ===========================================================================================
############################
# reassign environmental variables that are originally assigned by default
# ======================================================================================================================

## directories' path <A1>
./xmlchange DIN_LOC_ROOT=/work/ur12229009/cesm_inputdata/                     # directory of input dataset
./xmlchange CIME_OUTPUT_ROOT=/work/ur12229009/cesm_output/scratch             # directory of scratch data
./xmlchange DOUT_S_ROOT=/work/ur12229009/cesm_output/archive/$CASENAME        # directory of final output data

## characteristics of your simulation <A1>
# simulate directly to the end
# ./xmlchange STOP_OPTION=ndays,STOP_N=1200

# Separate to multi-chunks
./xmlchange STOP_OPTION=ndays,STOP_N=300,RESUBMIT=3        # run the simulation in 4 separate chunks (300 each chunks)


## parallel computing <A2>
# multiprocessing (cores used) (set to -1 for only using cores in one single node)
./xmlchange NTASKS_ATM=56
./xmlchange NTASKS_LND=56
./xmlchange NTASKS_ICE=56
./xmlchange NTASKS_OCN=56
./xmlchange NTASKS_CPL=56
./xmlchange NTASKS_GLC=56
./xmlchange NTASKS_ROF=56

# multithreading
./xmlchange NTHRDS_ATM=1
./xmlchange NTHRDS_LND=1
./xmlchange NTHRDS_ICE=1
./xmlchange NTHRDS_OCN=1
./xmlchange NTHRDS_CPL=1
./xmlchange NTHRDS_GLC=1
./xmlchange NTHRDS_ROF=1


## set project and queue for submitting jobs on the machine (can also added after Step 3.) <A1>
./xmlchange PROJECT_REQUIRED=TRUE
./xmlchange PROJECT=$PROJECT_ID
./xmlchange JOB_QUEUE=$WHICH_QUEUE

##############################
# * Step 2 - set up the case # =========================================================================================
##############################
# creates more files and directories needed to run the model.
# The following directory will be created:
#    Build/Run Directory (path: specify by .xml files under $CASEDIR) -
#       where the CESM code are actually compiled (bld/) and executed (run/),
#       this is also the place the outputting data temporarily stored (run/)
# ======================================================================================================================

# setup your case
./case.setup


#################################
# Step 2.5 - More Modifications # ======================================================================================
#################################
# 1. Modify namelist to specify model's controls. ($CASEROOT/user_nl_*)
# 2. Modify source code if needed. (under $CASEROOT/SourceMods/)
# 3. Modify namelist definitions if needed. ($SRCROOT/components/cam/bld/namelist_files/namelist_definition.xml)
# 4. Select/Regenerate appropriate input data that fit the modified model configuration. ($DIN_LOC_ROOT)
#    [remember to reset the location (path) to which 'ncdata' points (i.e. reset 'myfilepath' in 'user_nl_cam' file),
#    so that the model inputs the correct initial condition file.]
# ======================================================================================================================

## Adding Namelist: specify output frequency and numbers of output samples per file <A3>
cat <<END_OF_INSERT >> user_nl_cam

! specify output frequency and numbers of output samples per file
! assign values in order of h0, h1, and h2, ... files
nhtfrq = 0, -24
mfilt = 1, 365

END_OF_INSERT


## Adding Namelist: adding output variables <A4>
cat <<END_OF_INSERT >> user_nl_cam

! adding output variables
! 'fincl<x>' controls the fields in h<x+1>
fincl1 = 'U', 'V', 'OMEGA', 'T', 'Z3', 'VZ', 'VT', 'VU', 'OMEGAV', 'OMEGAT', 'OMEGAU', 'PS', 'PSL', 'T010', 'U010', 'QRS', 'DTCORE'
fincl2 = 'U', 'V', 'OMEGA', 'T', 'Z3', 'VZ', 'VT', 'VU', 'OMEGAV', 'OMEGAT', 'OMEGAU', 'PS', 'PSL', 'T010', 'U010', 'QRS', 'DTCORE'

END_OF_INSERT


## Adding Namelist: specify namelist parameters to control the use of the PK02 stratosphere <PK1>
#  namelist -
#    pkstrat: set to true to use the Polvani and Kushner stratosphere (Polvani and Kushner, 2002)
#     vgamma: The gamma value (K/Km) in the Polvani and Kushner stratosphere (Polvani and Kushner, 2002)
#             that represents the lapse rate in the polar stratosphere.  Default=2.

cat <<END_OF_INSERT >> user_nl_cam

! PKSTRAT
! Add namelist parameters to control the use of the PK02 relaxation temperature profile
&pkstrat_nl
  pkstrat=.True.
  vgamma=2.
/
! END-PKSTRAT

END_OF_INSERT


## Source Code Modification for PK02 Stratosphere configuration <PK2>
cp $CESMDIR/code/PK02_configuration/pkstrat_CESM2_2_3/SourceMods/src.cam/held_suarez_1994.F90 $CASEDIR/SourceMods/src.cam/held_suarez_1994.F90
cp $CESMDIR/code/PK02_configuration/pkstrat_CESM2_2_3/SourceMods/src.cam/held_suarez_cam.F90 $CASEDIR/SourceMods/src.cam/held_suarez_cam.F90
cp $CESMDIR/code/PK02_configuration/pkstrat_CESM2_2_3/SourceMods/src.cam/runtime_opts.F90 $CASEDIR/SourceMods/src.cam/runtime_opts.F90


## Adding Namelist definitions for PK02 source code modifications <PK3>
cp $CESMDIR/code/PK02_configuration/pkstrat_CESM2_2_3/namelist_definitions/namelist_definition.xml $SRCDIR/components/cam/bld/namelist_files/namelist_definition.xml
# (since the modification is made right in $SRCROOT, it is only needed for the first-time running)


## Regenerate appropriate input data that fit the PK02 configuration. <PK4>

cp $CESMDIR/code/PK02_configuration/dabiic/* $DATADIR/atm/cam/inic/dabiic/
# (since the makeic.ncl and its related files is put in $DIN_LOC_ROOT, and will remain existed globally,
# this command is only needed for the first-time running)

nlev=60             # numbers of vertical layers

# replace context in makeic.ncl
sed -i "66c\nlev = $nlev" $DATADIR/atm/cam/inic/dabiic/makeic.ncl

# execute makeic.ncl - generate input data
cd $DATADIR/atm/cam/inic/dabiic/
module load coda/tools/ncl
ncl $DATADIR/atm/cam/inic/dabiic/makeic.ncl
cd $CASEDIR

# reassign number of levels to fit the actual vertical resolution
# (if the resolution is different from the default of $COMPSET)
./xmlchange --file env_build.xml --id CAM_CONFIG_OPTS --val "--phys held_suarez --nlev=$nlev"


## Adding Namelist: revalue 'ncdata' to the location of expected (new) initial condition file <PK5>
# remember to check ncdata path
cat <<END_OF_INSERT >> user_nl_cam

! input data (I.C.) for modified vertical resolution configuration
! ncdata='myfile'   ! 'myfile' for location of expected initial condition
! for T85L60
ncdata='/work/ur12229009/cesm_inputdata/atm/cam/inic/dabiic/ncdata/HS1994.128x256.L60.nc'

END_OF_INSERT


## Adding Namelist: revalue hyperparameters for modified resolution configuration <PK6>
cat <<END_OF_INSERT >> user_nl_cam

! hyperparameters for modified resolution configuration
! usually automatically set up as case is created,
! and is consistent to the horizontal resolution assigned when creating the case
! for T85
eul_hdif_order = 4
eul_hdif_coef = 7.14D14

END_OF_INSERT


## generate the namelists
./preview_namelists

########################################
# * Step 3 - build (compile) the model # ===============================================================================
########################################
# Compile the model that has already set up.
# It will write files in the build directory
# ======================================================================================================================

# build the final CESM model executable
./case.build

# if rebuild:
# ./case.build --clean-all   # cleaning bld directory
# ./case.build


##############################
# * Step 4 - submit the jobs # =========================================================================================
##############################
# After issuing the command, the system will:
# 1. Checking archive and run options
# 2. Checking in namelists that need to be rebuilt
# 3. Check that all of the necessary input data for running is available (downloading whatever is missing)
# 4. Submitting the "case.run" script to the HPC batch job scheduler
# 5. As "case.run" executes, it will perform (run) the actual simulation.
# ======================================================================================================================

# submits the experiment to start it running
./case.submit

#######################
# Output data storage #
#######################
# When the run is successfully finished, case.submit will submit the "case.st_archive" script to archive the model output.
# The data will then moved from the Run Directory to the Archive Directory.
