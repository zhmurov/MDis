/*
 * parameters.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>

#define PARAMETER_LENGTH						100
#define PARAMETER_STRING_UNDEFINED				"NONE"

/*
 * Configuration file parameters names
 */
#define PARAMETER_DEVICE						"device"

#define PARAMETER_REDIRECT_STDOUT               "stdout"
#define PARAMETER_REDIRECT_STDERR               "stderr"


#define PARAMETER_FF_TYPE						"paratype"
#define PARAMETER_FF_FILE						"parameters"
#define PARAMETER_TOP_FILE						"topologies"

#define PARAMETER_COORD_FILE					"coordinates"
#define PARAMETER_VEL_FILE						"velocities"

#define PARAMETER_TOPOLOGY_FILE					"structure"

#define PARAMETER_FIRSTSTEP						"firststep"
#define PARAMETER_NUMSTEPS						"numsteps"
#define PARAMETER_TOTALTIME						"totaltime"
#define PARAMETER_RUN							"run"
#define PARAMETER_FIRSTRUN						"firstrun"
#define PARAMETER_RUNNUM						"runnum"
#define PARAMETER_NSIM							"simmetrynumber"

#define PARAMETER_INTEGRATOR					"integrator"
#define PARAMETER_VALUE_INTEGRATOR_LEAPFROG		"leap-frog"
#define PARAMETER_VALUE_INTEGRATOR_STEEPEST_DESCENT	"sd"
#define PARAMETER_VALUE_INTEGRATOR_VELOCITY_VERLET	"velocity-verlet"
#define PARAMETER_TEMPERATURE					"temperature"
#define PARAMETER_TIMESTEP						"timestep"

#define PARAMETER_LANGEVIN_ON					"langevin"
#define PARAMETER_DAMPING						"damping"
#define PARAMETER_RSEED							"seed"

#define PARAMETER_RIGIDBONDS					"rigidbonds"
#define PARAMETER_VALUE_RIGIDBONDS_NONE			"none"
#define PARAMETER_VALUE_RIGIDBONDS_HYDROGEN		"bonh"
#define PARAMETER_VALUE_RIGIDBONDS_ALL			"bond"

#define PARAMETER_RIGIDTOL						"rigidtol"


#define PARAMETER_MIN_CHANGE_TIMESTEP			"sdChangeTimestep"
#define PARAMETER_MIN_MAXFORCE					"sdMaxForce"
#define PARAMETER_MIN_UPDATE_TIMESTEP_FREQ		"sdUpdateFreq"
#define PARAMETER_MIN_CHANGE_TIMESTEP_AT		"sdVelThreshold"
#define PARAMETER_MIN_TIMESTEP_FACTOR			"sdTSFactor"
#define PARAMETER_MIN_FINAL_TIMESTEP			"sdFinalTimestep"

#define PARAMETER_HEATING						"heating"
#define PARAMETER_INITIAL_TEMPERATURE			"initialtemp"
#define PARAMETER_FINAL_TEMPERATURE				"finaltemp"
#define PARAMETER_TEMPERATURE_UPDATE_FREQ		"tempfreq"
#define PARAMETER_TEMPERATURE_INCREMENT			"tempincr"

#define PARAMETER_PULLING						"pulling"
#define PARAMETER_PULLING_PROTOCOL				"pullProtocol"
#define PARAMETER_VALUE_PULLING_PROTOCOL_FCONST	"fconst"
#define PARAMETER_VALUE_PULLING_PROTOCOL_FRAMP	"smd"
#define PARAMETER_PULLING_REFFILE				"pullFile"
#define PARAMETER_PULLING_SMD_UPDATE_FREQ		"smdUpdateFreq"
#define PARAMETER_PULLING_SMD_FILE				"smdOutputFile"
#define PARAMETER_PULLING_SMD_RESTART			"smdRestart"
#define PARAMETER_PULLING_SMD_INITIAL_TIP_FILE	"initialTipPositionFile"
#define PARAMETER_PULLING_SMD_INITIAL_TIP_VECTOR        "initialTipPosition"
#define PARAMETER_PULLING_SMD_RVEL              "smdRescaleVelocity"
#define PARAMETER_PULLING_SMD_RKS               "smdRescaleKS"

#define PARAMETER_PULLING_PLANE					"pPulling"
#define PARAMETER_PULLING_PLANE_PROTOCOL		"pPullprotocol"
#define PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FCONST	"pfconst"
#define PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FRAMP	"psmd"
#define PARAMETER_PULLING_PLANE_POINT			"pPoint"
#define PARAMETER_PULLING_PLANE_NORM			"pNorm"
#define PARAMETER_PULLING_PLANE_PULLFORCE		"pForce"
#define PARAMETER_PULLING_PLANE_PULLSPEED		"pSpeed"
#define PARAMETER_PULLING_PLANE_PULLSPRING		"pSpring"
#define PARAMETER_PULLING_PLANE_MASS			"pMass"
#define	PARAMETER_PULLING_PLANE_FIXSPRING		"pFixspring"
#define PARAMETER_PULLING_PLANE_REFFILE			"pPullfile"
#define PARAMETER_PULLING_PLANE_UPDATE_FREQ		"pUpdateFreq"
#define PARAMETER_PULLING_PLANE_OUTPUT_FILE		"pOutputFile"
#define PARAMETER_PULLING_PLANE_OUTPUT_FREQ		"pOutputFreq"
//#define PARAMETER_PULLING_PLANE_CUTOFF		"pCutoff"

#define PARAMETER_PUSHING_PLANE					"pPushing"
#define PARAMETER_PUSHING_PLANE_PROTOCOL		"pPushprotocol"
#define PARAMETER_PUSHING_PLANE_COUNT			"pPushCount"
#define PARAMETER_PUSHING_PLANE_NORM_FIX		"pNormFix"
#define PARAMETER_PUSHING_PLANE_POSITION_FIX		"pPositionFix"
#define PARAMETER_PUSHING_PLANE_NORM_PUSH		"pNormPush"
#define PARAMETER_PUSHING_PLANE_POSITION_PUSH		"pPositionPush"
#define PARAMETER_PUSHING_PLANE_SIGMA			"pSigma"
#define PARAMETER_PUSHING_PLANE_EPSILON			"pEpsilon"
#define PARAMETER_PUSHING_PLANE_PUSHSPEED		"pSpeed"
#define PARAMETER_PUSHING_PLANE_PUSHSPRING		"pSpring"
#define PARAMETER_PUSHING_PLANE_PUSHDAMPING		"pDamping"
#define PARAMETER_PUSHING_PLANE_REF_FILE		"pRefFile"
#define PARAMETER_PUSHING_PLANE_UPDATE_FREQ		"pUpdateFreq"
#define PARAMETER_PUSHING_PLANE_OUTPUT_FILE		"pOutputFile"
#define PARAMETER_PUSHING_PLANE_OUTPUT_FREQ		"pOutputFreq"

#define PARAMETER_DRUM_POTENTIAL			"Drum"
#define PARAMETER_DRUM_POTENTIAL_SIGMA			"dSigma"
#define PARAMETER_DRUM_POTENTIAL_EPSILON		"dEpsilon"
#define PARAMETER_DRUM_POTENTIAL_RADIUS			"dRadius"
#define PARAMETER_DRUM_POTENTIAL_POINT			"dPoint"
#define PARAMETER_DRUM_POTENTIAL_NORM			"dNorm"
#define PARAMETER_DRUM_POTENTIAL_XYZ_OUTPUT_FILE	"dXyzOutput"


#define PARAMETER_POSSIBLEPAIRS_FREQ			"possiblepairsfreq"
#define PARAMETER_POSSIBLEPAIRS_CUTOFF			"possiblepairscutoff"

#define PARAMETER_PAIRS_FREQ					"pairsfreq"
#define PARAMETER_PAIRS_CUTOFF_LJ				"pairscutofflj"
#define PARAMETER_PAIRS_CUTOFF_COULOMB			"pairscutoffcoulomb"
#define PARAMETER_PAIRS_CUTOFF_NB				"pairscutoff"

#define PARAMETER_MAX_PAIRS_12					"maxpairs12"
#define PARAMETER_MAX_PAIRS_13					"maxpairs13"
#define PARAMETER_MAX_PAIRS_14					"maxpairs14"
#define PARAMETER_MAX_PAIRS_EXCLUDED			"maxpairsExcluded"
#define PARAMETER_MAX_PAIRS_LJ					"maxpairslj"
#define PARAMETER_MAX_PAIRS_COULOMB				"maxpairscoulomb"
#define PARAMETER_MAX_PAIRS_NB					"maxpairsnb"
#define PARAMETER_MAX_POSSIBLEPAIRS				"maxpossiblepairs"
#define PARAMETER_MAX_PAIRS_GB					"maxpairsgb"
#define PARAMETER_PAIRSLISTS_EXTENSION			"pairslistextension"

#define PARAMETER_NB_TYPE						"nbtype"
#define PARAMETER_VALUE_NB_TYPE_RDIESHIFT		"rdieshift"
#define PARAMETER_VALUE_NB_TYPE_CDIESWITCH		"cdieswitch"

#define PARAMETER_LJ_CUTOFF						"ljcutoff"
#define PARAMETER_LJ_SWITCH						"ljswitch"
#define PARAMETER_COULOMB_CUTOFF				"coulombcutoff"

#define PARAMETER_NB_CUTOFF						"nbcutoff"
#define PARAMETER_NB_SWITCH						"nbswitch"

#define PARAMETER_14_TYPE						"1-4"
#define PARAMETER_VALUE_14_EXCLUDE				"exclude"
#define PARAMETER_VALUE_14_SCALE				"scale"
#define PARAMETER_SCALE_14_FACTOR				"scale14"

#define PARAMETER_SASA_ON						"sasa"
#define PARAMETER_SASA_PAIRS_FREQ				"sasapairsfreq"
#define PARAMETER_SASA_PARAMETERS_FILE			"sasaparameters"
#define PARAMETER_SASA_PAIRSLIST_CUTOFF			"sasapairslistcutoff"
#define PARAMETER_SASA_RPROBE					"sasarprobe"
#define PARAMETER_SASA_PIJ_COV					"sasapijcov"
#define PARAMETER_SASA_PIJ_NB					"sasapijnb"
#define PARAMETER_MAX_PAIRS_SASA				"maxpairssasa"
#define PARAMETER_MAX_PAIRLIST_SASA				"maxpairlistsasa"

#define PARAMETER_GENBORN_ON					"genborn"
#define PARAMETER_GENBORN_PARAMETERS_FILE		"genbornparameters"
#define PARAMETER_GENBORN_EIN					"ein"
#define PARAMETER_GENBORN_EOUT					"eout"
#define PARAMETER_GENBORN_CUTALPHA				"cutalpha"

#define PARAMETER_GBSW_ON						"gbsw"
#define PARAMETER_GBSW_ANG_QUADR_FILENAME		"angularquadrfile"
#define PARAMETER_GBSW_RAD_QUADR_POINTCOUNT		"radialquadrn"
#define PARAMETER_GBSW_CUTOFF					"gbswcutoff"
#define PARAMETER_GBSW_SMOOTH_LENGTH			"gbswsmooth"
#define PARAMETER_GBSW_A0						"gbswa0"
#define PARAMETER_GBSW_A1						"gbswa1"
#define PARAMETER_GBSW_PBRADII_FILENAME			"gbswradiifile"
#define PARAMETER_MESH_FREQ						"meshfreq"
#define PARAMETER_MESH_MARGIN					"meshmargin"
#define PARAMETER_MESH_STEP						"meshstep"
#define PARAMETER_MESH_EXTENSION				"meshextension"

#define PARAMETER_GB_ON							"gb"
#define PARAMETER_GB_MODEL						"gbmodel"
#define PARAMETER_GB_PARAMETERS_FILE			"gbparameters"
#define PARAMETER_GB_DIELECTRIC_OFFSET			"gbdielectricoffset"
#define PARAMETER_GB_IONS						"gbions"
#define PARAMETER_GB_OBC_ALPHA					"gbalpha"
#define PARAMETER_GB_OBC_BETA					"gbbeta"
#define PARAMETER_GB_OBC_GAMMA					"gbgamma"
#define PARAMETER_SA_ON							"sa"
#define PARAMETER_SA_SIGMA						"sasigma"
#define PARAMETER_SA_RPROBE						"sarprobe"

#define PARAMETER_REPULSIVE_BOUNDARY_ON			"repulsiveboundary"
#define PARAMETER_REPULSIVE_BOUNDARY_FOR_MD		"boundary_for_md"
#define PARAMETER_REPULSIVE_BOUNDARY_GEOMETRY	"boundarygeometry" //redundant 3-01-2010
#define PARAMETER_REPULSIVE_BOUNDARY_RADIUS		"boundaryradius"
#define PARAMETER_REPULSIVE_BOUNDARY_LOCATION 	"boundarylocation"
#define PARAMETER_REPULSIVE_BOUNDARY_EPSILON 	"boundaryepsilon"
#define PARAMETER_REPULSIVE_BOUNDARY_SIGMA	 	"boundarysigma"
#define PARAMETER_REPULSIVE_BOUNDARY_PDB		"boundaryfile"

#define PARAMETER_PBC_ON                        "pbc"
#define PARAMETER_PBC_LX                        "pbc_lx"
#define PARAMETER_PBC_LY                        "pbc_ly"
#define PARAMETER_PBC_LZ                        "pbc_lz"

#define PARAMETER_ENERGYOUTPUT_FILE				"outputname"
#define PARAMETER_ENERGYOUTPUT_FREQ				"outputfreq"
#define PARAMETER_ENERGYOUTPUT_LEGEND_FILE		"outputlegend"
#define PARAMETER_DCDOUTPUT_FILE				"dcdfile"
#define PARAMETER_DCDOUTPUT_FREQ				"dcdfreq"

#define PARAMETER_RESTART_COORD_OUTPUT_FILE		"restartcoord"
#define PARAMETER_RESTART_VEL_OUTPUT_FILE		"restartvel"
#define PARAMETER_FINAL_COORD_OUTPUT_FILE		"finalcoord"
#define PARAMETER_FINAL_VEL_OUTPUT_FILE			"finalvel"
#define PARAMETER_RESTARTOUTPUT_FREQ			"restartfreq"
#define PARAMETER_RESTART_KEY			        "restartkey"
#define PARAMETER_RESTART_COUNT	                "restartcount"

#define PARAMETER_PARAMETRIZE_SOP				"createSOP"

#define PARAMETER_RB_ENABLED          "RB_Enabled"
#define PARAMETER_RB_DEBUG            "RB_Debug"
#define PARAMETER_RB_CENTRALIZE       "RB_Centralize"
#define PARAMETER_RB_FREQ             "RB_Freq"
#define PARAMETER_RB_TIME_STEP        "RB_Timestep"
#define PARAMETER_RB_CUTOFF_DISTANCE  "RB_CutoffDistance"
#define PARAMETER_RB_CUTOFF_BARRIER   "RB_CutoffBarrier"
#define PARAMETER_RB_GROUPS_PDB       "RB_GroupsPDB"
#define PARAMETER_RB_RADII_TYPE       "RB_Radii_Type" // possible values: sup, avgsq, gyr
#define PARAMETER_RB_OUTPUT_FREQ      "RB_outputfreq"
#define PARAMETER_RB_MAXSTEP_MULT      "RB_MaxStepMultiplier"


#define PARAMETER_REMD_ENABLED           "REMD_Enabled"
#define PARAMETER_REMD_DEBUG             "REMD_Debug"
#define PARAMETER_REMD_DISABLE_EXCHANGES "REMD_DisableExchanges"
#define PARAMETER_REMD_HEATSTEPS         "REMD_HeatSteps"
#define PARAMETER_REMD_FREQ              "REMD_Freq"
#define PARAMETER_REMD_MIN_TEMPERATURE   "REMD_MinTemperature"
#define PARAMETER_REMD_MAX_TEMPERATURE   "REMD_MaxTemperature"
#define PARAMETER_REMD_EXCHANGESFILE         "REMD_ExchangesFile"

#define PARAMETER_CONSTRAINTS_ENABLED		"constraints"
#define PARAMETER_CONSTRAINTS_REF_FILENAME	"constraintsfile"
#define PARAMETER_CONSTRAINTS_KS			"constraintsks"
#define PARAMETER_CONSTRAINTS_GROUP			"constraintsgroup"
#define PARAMETER_CONSTRAINTS_CUTOFF		"constraintscutoff"

#define PARAMETER_VALUE_CONSTRAINTS_GROUP_IN	"in"
#define PARAMETER_VALUE_CONSTRAINTS_GROUP_BW	"bw"
#define PARAMETER_VALUE_CONSTRAINTS_GROUP_NONE 	"none"

#define PARAMETER_SC_ON 						"structconstr"

#define DEFAULT_LANGEVIN_ON						1
#define DEFAULT_TEMPERATURE						300.0f
#define DEFAULT_DAMPING							50.0f
#define DEFAULT_HEATING							0
#define DEFAULT_INITIAL_TEMPERATURE				0.0f
#define DEFAULT_FINAL_TEMPERATURE				300.0f
#define DEFAULT_MIN_MAXFORCE					10000.0f

#define DEFAULT_PULLING							0

#define DEFAULT_PULLING_PLANE					0
#define DEFAULT_PUSHING_PLANE					0
#define DEFAULT_DRUM_POTENTIAL					0

#define DEFAULT_PAIRS_CUTOFF_LJ					0.8f
#define DEFAULT_PAIRS_CUTOFF_COULOMB			0.8f
#define DEFAULT_PAIRS_CUTOFF_NB					0.8f
#define DEFAULT_LJ_CUTOFF						0.75f
#define DEFAULT_LJ_SWITCH						0.65f
#define DEFAULT_COULOMB_CUTOFF					0.75f
#define DEFAULT_POSSIBLEPAIRS_CUTOFF			1.6f
#define DEFAULT_PAIRSLISTS_EXTENSION			256

#define DEFAULT_SASA_RPROBE						0.14f
#define DEFAULT_SASA_PIJ_COV					0.8875f
#define DEFAULT_SASA_PIJ_NB 					0.3516f

#define DEFAULT_GENBORN_EIN						1.0f
#define DEFAULT_GENBORN_EOUT					80.0f
#define DEFAULT_GENBORN_CUTALPHA				1.0e5f

#define DEFAULT_GB_MODEL						"obc"
#define DEFAULT_GB_DIELECTRIC_OFFSET			0.009f
#define DEFAULT_GB_OBC_ALPHA					1.0f
#define DEFAULT_GB_OBC_BETA						0.8f
#define DEFAULT_GB_OBC_GAMMA					4.85f

#define DEFAULT_SA_ON							1
#define DEFAULT_SA_SIGMA						2.259357506f
#define DEFAULT_SA_RPROBE						0.14f

#define DEFAULT_GBSW_CUTOFF						2.0f
#define DEFAULT_GBSW_SMOOTH_LENGTH				0.01f
#define DEFAULT_GBSW_A0							-0.081f
#define DEFAULT_GBSW_A1							1.6f
#define DEFAULT_GBSW_RAD_QUADR_POINTCOUNT		0

#define DEFAULT_MESH_MARGIN						2.0f
#define DEFAULT_MESH_STEP						0.1f
#define DEFAULT_MESH_EXTENSION					1000000

#define DEFAULT_RESTART_COUNT	            0
/*
 * Output
 */
#define ENERGY_OUTPUT_WIDTH					15
#define ENERGY_OUTPUT_TITLE					"ENERGY:"

#define ENERGY_OUTPUT_NAME_HARMONIC			"Harmonic"
#define ENERGY_OUTPUT_NAME_UREY_BRADLEY		"Urey-Bradley"
#define ENERGY_OUTPUT_NAME_ANGLE			"Angles"
#define ENERGY_OUTPUT_NAME_DIHEDRAL 		"Dihedrals"
#define ENERGY_OUTPUT_NAME_IMPROPER			"Impropers"

#define ENERGY_OUTPUT_NAME_LENNARD_JONES	"VdW"
#define ENERGY_OUTPUT_NAME_ELECTROSTATIC	"Coulomb"

#define ENERGY_OUTPUT_NAME_SASA				"SASA"
#define ENERGY_OUTPUT_NAME_GEN_BORN			"GenBorn"
#define ENERGY_OUTPUT_NAME_GB				"GB"
#define ENERGY_OUTPUT_NAME_SA				"SA"
#define ENERGY_OUTPUT_NAME_SC				"StrConstr"

#define ENERGY_OUTPUT_NAME_TEMPERATURE		"Temperature"
#define ENERGY_OUTPUT_NAME_POTENTIAL		"TOTAL"

#define ENERGY_OUTPUT_NAME_TRAJECTORY		"Trajectory"
#define ENERGY_OUTPUT_NAME_TIMESTEP			"Step"
#define ENERGY_OUTPUT_NAME_TIME				"Time"

#define PULLING_OUTPUT_WIDTH				16
#define PULLING_OUTPUT_TITLE				"SMD:"

#define PULLING_OUTPUT_NAME_TRAJECTORY		"Trajectory"
#define PULLING_OUTPUT_NAME_ATOM			"Atom"
#define PULLING_OUTPUT_NAME_STEP			"Step"

#define PULLING_OUTPUT_NAME_FEXT			"Total Force"
#define PULLING_OUTPUT_NAME_FX				"Fx"
#define PULLING_OUTPUT_NAME_FY				"Fy"
#define PULLING_OUTPUT_NAME_FZ				"Fz"

#define PULLING_OUTPUT_NAME_R				"R"
#define PULLING_OUTPUT_NAME_RX				"Rx"
#define PULLING_OUTPUT_NAME_RY				"Ry"
#define PULLING_OUTPUT_NAME_RZ				"Rz"

#define PARAMETER_MPI_RANK                  "mpi_rank"
#define PARAMETER_MPI_DEV_AUTO              "mpi_device_auto"
#define PARAMETER_MPI_DEVPERNODE            "mpi_dpn"
#define PARAMETER_MPI_DEVICE(i)             parameter_mpi_device(i)
#define PARAMETER_MPI_DEV_CUR               "mpi_device_current"
#define PARAMETER_MPI_FIRSTRUN              "mpi_firstrun_auto"

#define PARAMETER_UMBRELLA                  "umbrella"
#define PARAMETER_UMBRELLA_DEBUG            "umbrella_debug"
#define PARAMETER_UMBRELLA_OUTFILE          "umbrella_file"
#define PARAMETER_UMBRELLA_STAGE            "umbrella_stage"
#define PARAMETER_UMBRELLA_METHOD           "umbrella_method"
#define PARAMETER_UMBRELLA_WINSTEP          "umbrella_win_step"
#define PARAMETER_UMBRELLA_KS               "umbrella_ks"
#define PARAMETER_UMBRELLA_ORTHOGONAL_KS    "umbrella_ks_o"
#define PARAMETER_UMBRELLA_RESTRAIN_KS      "umbrella_ks_r"
#define PARAMETER_UMBRELLA_FREQ             "umbrella_freq"
#define PARAMETER_UMBRELLA_COM_CHAIN_FIXED  "umbrella_com_fixed"
#define PARAMETER_UMBRELLA_COM_FIXED_ALL    "umbrella_com_fixall"
#define PARAMETER_UMBRELLA_COM_FIXED_LIST   "umbrella_com_fixedlist"
#define PARAMETER_UMBRELLA_COM_FIXED_MASK   "umbrella_com_fmask"
#define PARAMETER_UMBRELLA_COM_CHAIN_PULLED "umbrella_com_pulled"
#define PARAMETER_UMBRELLA_COM_PULLED_MASK  "umbrella_com_pmask"
#define PARAMETER_UMBRELLA_COM_USEGPU       "umbrella_com_usegpu"
#define PARAMETER_UMBRELLA_COM_FIX_COM      "umbrella_com_fix_com"
#define PARAMETER_UMBRELLA_COM_LOAD_VECTOR  "umbrella_com_load_vector"
#define PARAMETER_UMBRELLA_COM_COORD_ZERO   "umbrella_com_coord_zero"
#define PARAMETER_UMBRELLA_DIRECTION        "umbrella_direction"    
#define PARAMETER_UMBRELLA_VELOCITY         "umbrella_velocity"
#define PARAMETER_UMBRELLA_FILE_COM_TEMP    "umbrella_file_temp"
#define PARAMETER_UMBRELLA_FILE_FRAME       "umbrella_file_frame"
#define PARAMETER_VALUE_UMBRELLA_STAGE_PREPARE      "preparation"
#define PARAMETER_VALUE_UMBRELLA_STAGE_PRODUCTION   "production"
#define PARAMETER_VALUE_UMBRELLA_METHOD_COM         "com"
#define PARAMETER_VALUE_UMBRELLA_METHOD_ROT         "rot"


typedef struct {

	int device;

	int firstrun;
	int Ntr;
	float initialTime;
	float totalTime;
	long long firststep;
	long long int numsteps;
	int rseed;

	// Files
	char topologyFilename[PARAMETER_LENGTH];
	char coordFilename[PARAMETER_LENGTH];
	char velFilename[PARAMETER_LENGTH];
	char forceFieldFilename[PARAMETER_LENGTH];
	char topologiesFilename[PARAMETER_LENGTH];
	char parametersType[PARAMETER_LENGTH];
	char dcdFilename[PARAMETER_LENGTH];
	char restartKey[PARAMETER_LENGTH];

} Parameters;

void parseParametersFileNAMD(char* filename, Parameters* parameters);
FILE* parametersPrepareRestart(Parameters* parameters);
char *parameter_mpi_device(int i);
