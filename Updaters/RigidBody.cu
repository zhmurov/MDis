#include "RigidBody.cuh"
#include <set>
#include <iterator>
#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Potentials/RepulsiveBoundaryPotential.cuh"
#include "CoordinatesOutputManagerDCD.cuh"
#include "EnergyOutputManager.cuh"
#include "../Util/Log.h"
#include "../Util/wrapper.h"
#include "../Util/mystl.h"

namespace rigid_body {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<rigid_body> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create() {
	debug = (bool)getYesNoParameter(PARAMETER_RB_DEBUG, 0);
	if (!getYesNoParameter(PARAMETER_RB_ENABLED, 0))
		return;
	LOG << "create";
	updater.update = update;
	updater.destroy = destroy;
	updater.frequency = getIntegerParameter(PARAMETER_RB_FREQ);
	sprintf(updater.name, "Brownian Rigid Body");
	updaters[updatersCount] = &updater;
	updatersCount ++;
    // read radii_type from config
    char radii_type_s[256];
    getParameter(radii_type_s, PARAMETER_RB_RADII_TYPE, "sup");
    radii_type = -1;
    if ( !strcmp(radii_type_s, "sup"  ) )  radii_type = 0;
    if ( !strcmp(radii_type_s, "avgsq") )  radii_type = 1;
    if ( !strcmp(radii_type_s, "gyr"  ) )  radii_type = 2;
    if (radii_type == -1) {
        LOG << "WARNING: Unknown radii type '" << radii_type_s << "'. Assuming 'sup'";
        radii_type = 0;
    }

    maxstep_mult = getIntegerParameter(PARAMETER_RB_MAXSTEP_MULT, 100000);
	init();
	LOG << "create done";
}

inline std::string getGroup(bool usePDB, PDB& pdb, int index) {
	if (usePDB)
		return any2str(pdb.atoms[index].occupancy);
	return topology.atoms[index].segment;
}

inline void writeRange(std::ostream& ranges, int& first, int second) {
	if (first == -1)
		return;
	if (first == second)
		ranges << " " << first;
	else
		ranges << " " << first << "-" << second;
	first = -1;
}

void init() {
	rseed = getIntegerParameter(PARAMETER_RSEED);
	groups.resize(gsystem.N);
    coordinateSystem.resize(parameters.Ntr, make_float4(0.0f, 0.0f, 0.0f, 0.0f));

	std::set<std::string> groupsSet;
	PDB groupsPDB;
	char groupsPDBFilename[256];
	getMaskedParameter(groupsPDBFilename, PARAMETER_RB_GROUPS_PDB, "-");
	bool usePDB = std::string(groupsPDBFilename) != "-";
	if (usePDB) {
		// Read PDB
		readPDB(groupsPDBFilename, &groupsPDB);
		if (groupsPDB.atomCount != gsystem.N) {
			DIE("Groups file contain wrong number of lines. It has %d atoms and %d atoms are required!", groupsPDB.atomCount, gsystem.N);
		}
		// Find all flags
	}
	for (int i = 0; i < gsystem.N; i++)
		groupsSet.insert(getGroup(usePDB, groupsPDB, i));
	bodyCount = groupsSet.size();
	body.resize(bodyCount);

	for (std::set<std::string>::iterator i = groupsSet.begin(), e = groupsSet.end(); i != e; ++i) {
		std::stringstream ranges;
		int first = -1, second;
		for (int id = 0; id < gsystem.N; id++) {
			if (getGroup(usePDB, groupsPDB, id) == *i) {
				if (first == -1)
					first = second = id;
				else
					second = id;
			} else {
				writeRange(ranges, first, second);
			}
		}
		writeRange(ranges, first, second);
		LOG << "Found body with " << (usePDB? "\"occupancy\" ": "segment ") << *i << ". AtomIds:" << ranges.str();
	}

	// Fill groups and find out bodies' masses
	for (int j = 0; j < bodyCount; j++)
		body[j].mass = 0.0f;
	for (int i = 0; i < gsystem.N; i++) {
		int j = groups[i] = std::distance(groupsSet.begin(), groupsSet.find(getGroup(usePDB, groupsPDB, i)));
		body[j].mass += topology.atoms[i].mass;
	}

	// Make a copy of groups on device (needed for write back LD->MD procedure)
	cudaMalloc((void**)&d_groups, sizeof(int) * groups.size());
	cudaMemcpy(d_groups, &groups[0], sizeof(int) * groups.size(), cudaMemcpyHostToDevice);

	// Calculate viscosity of (assume) water using empirical equations
	T = getFloatParameter(PARAMETER_TEMPERATURE, DEFAULT_TEMPERATURE);
	timeStep = getFloatParameter(PARAMETER_RB_TIME_STEP, 100);
	cutoffDistance = getFloatParameter(PARAMETER_RB_CUTOFF_DISTANCE, 1.5f);
	cutoffBarrier = getFloatParameter(PARAMETER_RB_CUTOFF_BARRIER, 0.2f);
	viscosity = calculateViscosity(T);
}

void destroy() {}

void readTrajectory(int traj) {
	// Read rigid bodies and find centers of mass
	for (int j = 0; j < bodyCount; j++)
		body[j].r = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	for (int i = 0, ti = traj * gsystem.N; i < gsystem.N; i++, ti++) {
		int j = groups[i];
		body[j].r.x += gsystem.h_coord[ti].x*topology.atoms[i].mass;
		body[j].r.y += gsystem.h_coord[ti].y*topology.atoms[i].mass;
		body[j].r.z += gsystem.h_coord[ti].z*topology.atoms[i].mass;
	}
	for (int j = 0; j < bodyCount; j++) {
		body[j].r.x /= body[j].mass;
		body[j].r.y /= body[j].mass;
		body[j].r.z /= body[j].mass;
        body[j].r.w = 0.0;
		body[j].r0 = body[j].r;
		body[j].q = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		body[j].R = 0.0f;
	}
	// Find bodies' radii
    for (int i = 0, ti = traj * gsystem.N; i < gsystem.N; i++, ti++) {
		int j = groups[i];
        float rt = 0.0f;
        switch (radii_type) {
            case 0:
        		body[j].R = std::max(body[j].R, findDistance(body[j].r, gsystem.h_coord[ti]));
                break;
            case 1:
                rt = findDistance(body[j].r, gsystem.h_coord[ti]);
                body[j].R += rt*rt;
                body[j].r.w += 1;
                break;
            case 2:
                for (int k = i + 1, tk = traj * gsystem.N + i + 1; k < gsystem.N; k++, tk++) {
                    if (groups[i] != groups[k]) continue;
                    rt = findDistance(gsystem.h_coord[ti], gsystem.h_coord[tk]);
                    body[j].R += rt*rt;
                }
                body[j].r.w += 1;
                break;
        }
	}
    if (radii_type == 1)
        for (int j = 0; j < bodyCount; j++)
            body[j].R = sqrtf(body[j].R / body[j].r.w);
    if (radii_type == 2)
        for (int j = 0; j < bodyCount; j++)
            body[j].R = sqrtf(body[j].R / ( body[j].r.w * body[j].r.w ));
	// Update parameters of bodies which depend on body radius
	for (int j = 0; j < bodyCount; j++)
		body[j].updateParameters();
	// Log parameters found
	for (int j = 0; j < bodyCount; j++) {
		LOG << "TRAJ" << traj << " BODY" << j
			  << " R=" << body[j].R
			  << " ksi_trans=" << body[j].ksi_trans
			  << " var_trans=" << body[j].var_trans
			  << " ksi_rot=" << body[j].ksi_rot
			  << " var_rot=" << body[j].var_rot;
	}
}

__global__ void moveBodies(unsigned size, int* groups, float4* curPosition, float4* oldPosition, Matrix<float, 3, 3>* rotation, float4* coord) {
	unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= size)
		return;
	int j = groups[i];
	Vector<float, 3>& r = AsVector<3>(coord[i]);
	Vector<float, 3>& r0 = AsVector<3>(oldPosition[j]);
	Vector<float, 3>& r1 = AsVector<3>(curPosition[j]);
	r = r1 + rotation[j] * (r - r0);
}

void centralizeTrajectory(int traj) {
	// Stop if centralization disabled
	if (!getYesNoParameter(PARAMETER_RB_CENTRALIZE, 0))
		return;

	LOG << "TRAJ" << traj << " centralize";

	// Find a center
	float4 center = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	for (int j = 0; j < bodyCount; j++)
		center += body[j].r;
	center /= bodyCount;

	// Set new coordinate system to the center
	coordinateSystem[traj] += center;
	LOG << "TRAJ" << traj << " centralized. coordinate system center in now at point " << coordinateSystem[traj];

	// Adjust bodies' position according to changes above
	for (int j = 0; j < bodyCount; j++) {
		body[j].r -= center;
    LOG << "TRAJ" << traj << " BODY" << j << " center at " << body[j].r;
  }
}

void writeTrajectory(int traj) {
	static std::vector< float4 > h_curPosition(bodyCount);
	static std::vector< float4 > h_oldPosition(bodyCount);
	static std::vector< Matrix<float, 3, 3> > h_rotation(bodyCount);
	for (int j = 0; j < bodyCount; j++) {
		h_curPosition[j] = body[j].r;
		h_oldPosition[j] = body[j].r0;
		quaternion::toMatrix(body[j].q, h_rotation[j]);
		// Reset previous location to current location
		body[j].r0 = body[j].r;
		body[j].q = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	}

	static float4* d_curPosition = 0;
	static float4* d_oldPosition = 0;
	static Matrix<float, 3, 3>* d_rotation = 0;
	if (!d_curPosition) {
		cudaMalloc((void**)&d_curPosition, sizeof(float4) * bodyCount);
		cudaMalloc((void**)&d_oldPosition, sizeof(float4) * bodyCount);
		cudaMalloc((void**)&d_rotation, sizeof(Matrix<float, 3, 3>) * bodyCount);
	}

	cudaMemcpy(d_curPosition, &h_curPosition[0], sizeof(float4) * bodyCount, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldPosition, &h_oldPosition[0], sizeof(float4) * bodyCount, cudaMemcpyHostToDevice);
	cudaMemcpy(d_rotation, &h_rotation[0], sizeof(Matrix<float, 3, 3>) * bodyCount, cudaMemcpyHostToDevice);

	CUDA_SIMPLE_INVOKE(moveBodies, gsystem.N)
			(gsystem.N, d_groups, d_curPosition, d_oldPosition, d_rotation, gsystem.d_coord + traj * gsystem.N);
}

void calcUnboundEnergy() {
    copyCoordinatesFromGPU();
    float4 old_coords[parameters.Ntr*gsystem.N];
    // back up coordinates
    memcpy(old_coords, gsystem.h_coord, sizeof(float4)*parameters.Ntr*gsystem.N);
    for (int traj = 0; traj < parameters.Ntr; traj++) {
        readTrajectory(traj);
        float x = 0;
        float d = cutoffDistance + cutoffBarrier + 5.0f;
        body[0].r = make_float4(x,0,0,0);
        for (int i = 1; i < bodyCount; ++i ) {
            x += body[i-1].R + d + body[i].R;
            body[i].r = make_float4(x,0,0,0);
        }
        for (int i = 0; i < bodyCount; ++i ) {
            body[i].q = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
        }
        writeTrajectory(traj);
    }
    // TODO: calculate energy
    // restore coordinates
    memcpy(gsystem.h_coord, old_coords, sizeof(float4)*parameters.Ntr*gsystem.N);
    for (int traj = 0; traj < parameters.Ntr; traj++) {
        copyCoordinatesToGPU(traj, 0);
    }
}

void update() {
	LOG << "update";
    static int dcdOutputFreq = 0;
    int logFreq = 100;
	if (!dcdOutputFreq) {
		if (updater.frequency % possiblePairsListUpdater.frequency != 0) {
			DIE("rigid body updater frequency (%d) is not a multiple of pairs list updater frequency (%d)", updater.frequency,  possiblePairsListUpdater.frequency);
		}
		if (energy_output::updater.frequency != coordinates_output_dcd::updater.frequency) {
			DIE("energy_output::updater.frequency != coordinates_output_dcd::updater.frequency");
		}
		dcdOutputFreq = (int)ceil(coordinates_output_dcd::updater.frequency * getFloatParameter(PARAMETER_TIMESTEP) / timeStep);
		if (dcdOutputFreq < 100)
			dcdOutputFreq = 100;
		LOG << "dcdOutputFreq = " << dcdOutputFreq;
	}
    const int rb_dcdOutputFreq = getIntegerParameter(PARAMETER_RB_OUTPUT_FREQ, dcdOutputFreq);
	copyCoordinatesFromGPU();
	for (int traj = 0; traj < parameters.Ntr; traj++) {
		LOG << "TRAJ" << traj << " begin";
		readTrajectory(traj);

		// Main Langevin Dynamics cycle
		int stepCount = 0;
		float startTime = trajectoryTime[traj];
		long long maxSteps = maxstep_mult*ceil((parameters.totalTime - startTime) / timeStep);
		for (; stepCount < maxSteps; stepCount++) {
			// Check distances
			bool finish = false;
			float minDistance = 1e10f;
			float range = cutoffDistance + (stepCount == 0? 0: cutoffBarrier); // adjust range for the first time to avoid blinking MD<->LD
			for (int i = 0; !finish && i < bodyCount - 1; i++) {
				for (int j = i + 1; !finish && j < bodyCount; j++) {
					float distance = findDistance(body[i].r, body[j].r);
                    //LOG << "Distance between bodies (" << i << "," << j << ") = " << distance;
					minDistance = std::min(minDistance, distance);
					if (distance < body[i].R + body[j].R + range)
						finish = true;
				}
			}

			// Log current state
			if (stepCount % logFreq == 0)
				LOG << "TRAJ" << traj << " step = " << stepCount << " maxSteps = " << maxSteps << " minDistance = " << minDistance;

			// Finish (go to MD) if we have at least one pair of bodies closer than cutoff range
			if (finish)
				break;

			// Move bodies
			for (int j = 0; j < bodyCount; j++)
				body[j].step(traj);

			// Write DCD files if needed
			if (stepCount % rb_dcdOutputFreq == 0) {
				trajectoryTime[traj] = startTime + stepCount * timeStep;
				writeTrajectory(traj);
				copyCoordinatesFromGPU(true);
				coordinates_output_dcd::saveTrajectory(traj);
				energy_output::saveRigidBody(traj);
			}
		}

		// Dont write back if we dont change anything
		if (stepCount == 0) {
			LOG << "TRAJ" << traj << " nothing to do";
			continue;
		}

		// Save simulation results
		trajectoryTime[traj] = startTime + stepCount * timeStep;
		centralizeTrajectory(traj);
		writeTrajectory(traj);
		LOG << "TRAJ" << traj << " " << stepCount << " steps done";

		// TODO: somehow tell to MD that we should stop if timeLimit have been reached
	}
}

void RigidBody::addRepulsiveBoundaryForce(int traj) {
	float4 coord = r;
	float4 sphere = repulsiveBoundaryData.geometry;
    sphere += coordinateSystem[traj];
	coord.x -= sphere.x;
	coord.y -= sphere.y;
	coord.z -= sphere.z;
	coord.w = sqrtf(coord.x*coord.x + coord.y*coord.y + coord.z*coord.z);
	sphere.w -= coord.w;
	sphere.w -= repulsiveBoundaryData.sigma;
	if (sphere.w < 0) {
		float mult = 2.0f*repulsiveBoundaryData.epsilon*sphere.w/coord.w;
		f.x += mult*coord.x;
		f.y += mult*coord.y;
		f.z += mult*coord.z;
	}
}

void RigidBody::step(int traj) {
	zeroForceAndTorque();
	addRandomForce();
    static int bounded = getYesNoParameter(PARAMETER_REPULSIVE_BOUNDARY_ON, 0);
    if (bounded)
        addRepulsiveBoundaryForce(traj);
	addRandomTorque();
	integrate();
}

#undef LOG

} // namespace rigid_body
