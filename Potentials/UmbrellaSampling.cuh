/*
 * UmbrellaSampling.cuh
 *
 *  Created on: Jul 21, 2011
 *      Author: alekseenko
 *
 *  There is a lot of info about this method
 *  I would personally recommend Kästner, 2011 (doi:10.1002/wcms.66) and Buch et al., 2010 (doi:10.1021/ci900455r)
 *
 */

#pragma once


namespace umbrella_sampling {

int blockSize;
int blockCount;
int blockCountTr;

#define UMBRELLA_ATOM_GROUP_OTHER 0
#define UMBRELLA_ATOM_GROUP_FIXED 1
#define UMBRELLA_ATOM_GROUP_PULLED 2

#define UMBRELLA_METHOD_COM 1
#define UMBRELLA_METHOD_ROT 2

#define UMBRELLA_STAGE_PREPARATION 0
#define UMBRELLA_STAGE_PRODUCTION 1


typedef struct {
	int* atomGroups; //[N] What this mean depends on choice of RCs
	float ks; // Spring constant. We use harmonic potential
	float ks_ortho; // Spring constant for perpendicular movement
	float ks_rest; // K-spring for restraining fixed atoms
	float4 d; // Direction of pulling, normalized
	float v; // Velocity of pulling
} Data;

typedef struct {
	float4 *rcs; //[Ntr] Reaction coordinates (RCs). For CoM sampling they are CoM_para, CoM_orth, CoM_par_norm.
	float2 *energy; //[Ntr] Umbrella energy, parallel and orthogonal
	int method; // What Reaction Coordinates (RCs) to use. Now only 'CoM' is supported
	//    int win_count; // Number of windows
	float win_step; // Distance between windows
	int stage; // Whether we are preparing initial positions or actually sampling the energy landscape
	std::string energyfile; // File to save energies (and RCs) to
} SamplingParameters;

typedef struct {
	float4 *positions; //[Ntr] Desired positions of pulled chain
	float4 *d_dfp, *h_dfp; //[Ntr] Forces acting on pulled chain
	float4 initial_pos; // Coordinates of CoM of pulled chain at t=0
	float  rc_zero; // The adjusstment to trajectories desired RCs during pulling
	float4 initial_fixpos; // Coordinates of CoM of fixed chain at t=0, necessary for alt_fix
	float4 *d_initial_coords; // Initial coords of the system to restrain ligand to
	bool use_gpu; // Use GPU for CoM calculations
	bool fix_com; // Instead of individually fixing FIXED atoms, fix their com
	float4 *d_com, *h_com;
} CoMParameters;

typedef struct {
	float4 *d_initial_coords; // Initial coords of the system to restrain ligand to
	float4 center; // Rotational center
	// float4 axis; // Rotational center
	float *d_phi0;
	float3 *d_energy, *h_energy; // (0.5*k*phi^2, 0.5*k_o*dr^2, phi);
	int n_pulled;
} RotParameters;

Data h_data, d_data;
CoMParameters com_params;
RotParameters rot_params;
SamplingParameters sampling_params;

Potential potential;
Updater updater;
bool debug;

void create();
void init();
void initCoM();
//void initRot();
inline void compute();
void destroy();

void updateUmbrellaCoM();
void updateUmbrellaRot();
void updaterDestroy();
void computeUmbrellaCoM(int dry);
//void computeUmbrellaRot(int dry);
// Due to using pointers to computeX functions of type void(*)(), we can't simply specify default values :(
void computeUmbrellaCoM() {computeUmbrellaCoM(0);}
//void computeUmbrellaRot() {computeUmbrellaRot(0);}

} // namespace umbrella_sampling

/*
 
                    MM  :::::::::....::::::::  MM
                   MMMM  :::::::::::::::::::  MMMM
                  MMMMMM  :::::::::::::::::  MMMMMMM
                MMMMMMMMM  :::::::::::::::  MMMMMMMMM
              MMMMMMMMMMMM  :::::::::::::  MMMMMMMMMMMM
            MMMMMMMMMMMMMM   :::::::::::   MMMMMMMMMMMMMM
         MMMMMMMMMMMMMMMMMM   ::::::::::  MMMMMMMMMMMMMMMMMM
      MMMMMMMMMMMMMMMMMMMMMM   ::::::::  MMMMMMMMMMMMMMMMMMMMM
         MMMMMMMMMMMMMMMMMMMM   ::::::  MMMMMMMMMMMMMMMMNMMMMM 
      ....    MMMMMMMMMMMMMMM   :::::  MMMMMMMMMMMMMMMMMM   ..::
       :::::..     MMMMMMMMMMM   :::   MMMMMMMMMMMMMM  ..::::::
       :::::::::...   MMMMMMMMM   :   MMMNMMMMM   ...::::::::::
       ::::::::::::::...   MMMMM     MMMMM   ...::::::::::::::
       :::::::::::::::::::..   M     M     :::::::::::::::::::
       :::::::::::::::::::::.           .:::::::::::::::::::::
       :::::::::::::::::      MM     M        ::::::::::::::::
       ::::::::::::      MMMMMMM  .  MMMMMMM      ::::::::::::.
       ::::::::      MMMMMMMMMM  .:.  MMMMMMMMMM     ::::::::::
      .:::      MMMMMMMMMMMMMM  .:::.  MMMMMMMMMMMMM      ::::::
            MMMMMMMMMMMMMMMMMM  :::::.  MMMMMMMMMMMMMMMMM     
       MMMMMMMMMMMMMMMMMMMMMM  :::::::.  MMMMMMMMMMMMMMMMMMMMMM
        MMMMMMMMMMMMMMMMMMMM  :::::::::  MMMMMMMMMMMMMMMMMMMMMM
          MMMMMMMMMMMMMMMMM   ::::::::::  MMMMMMMMMMMMMMMMM 
             MMMMMMMMMMMMM   :::::::::::.  MMMMMMMMMMMMM 
               MMMMMMMMMM  .:::::::::::::.  MMMMMMMMMM
                 MMMMMMMM .:::::::::::::::.  MMMMMMM
                  MMMMMM  :::::::::::::::::.  MMMMM
                    MMM  :::::::::::::::::::   MMM
                     M  ::::....     ...:::::  MM

*/
