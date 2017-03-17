#pragma once

#include <vector>
#include "../Util/LinearAlgebra.h"
#include "../Util/Constants.h"
#include "../Util/Quaternion.h"
#include "../Util/ran2.h"
#include "../Core/gsystem.h"

namespace rigid_body {
	struct RigidBody {
		float mass;
		float4 f;
		float4 M;
		float4 r;
		float4 r0;
		float4 q;
		float ksi_trans;
		float ksi_rot;
		float var_trans;
		float var_rot;
		float R;

		void updateParameters();
		void zeroForceAndTorque();
		void addRandomForce();
		void addRandomTorque();
		void addRepulsiveBoundaryForce(int traj);
		void integrate();
		void step(int traj);
		void reset(float4 r_, float mass_, float R_);
	};

	std::vector<float4> coordinateSystem;
	std::vector<RigidBody> body; // Container for all rigid bodies
	size_t bodyCount;      // Number of bodies (the same for each trajectory)
	float cutoffDistance;  // If distance between any pair of bodies' surfaces greater than cutoffDistance, we continue MD
	float cutoffBarrier;   // If distance between any pair of bodies' surfaces less than cutoffDistance + cutoffBarrier, we enter LD
	float timeStep;        // LD time step

	std::vector<int> groups;
	int* d_groups;
	float viscosity;
	int rseed;
	float T;
    int radii_type, maxstep_mult;
	bool debug;

	Updater updater;

	void create();
	void init();
	void destroy();
	void update();

	inline float calculateViscosity(float T) {
		return unitssystem::Convert(unitssystem::si, unitssystem::gromacs,
			(1.78e-6f / (1 + 3.37e-2f * (T - 273) + 2.21e-4f * (T - 273) * (T - 273))) * // kinematic viscosity
			(9.957e2f / (0.984f + 0.483e-3f * (T - 273))), // density of water
			_M - _L - _T); // kg*m^-1*s^-1;
	}

	void RigidBody::updateParameters() {
		ksi_trans = 6 * M_PI * viscosity * R; //mass*Kb_MD*T/D;
		ksi_rot = 8 * M_PI * viscosity * R * R * R; //mass*Kb_MD*T/D;
		var_trans = sqrtf(2.0f*ksi_trans*constant::Boltzmann(unitssystem::gromacs)*T/timeStep);
		var_rot = sqrtf(2.0f*ksi_rot*constant::Boltzmann(unitssystem::gromacs)*T/timeStep);
	}

	void RigidBody::zeroForceAndTorque() {
		f.x = f.y = f.z = M.x = M.y = M.z = 0.0f;
	}

	void RigidBody::addRandomForce() {
		f.x += var_trans * ran2::gasdev(&rseed);
		f.y += var_trans * ran2::gasdev(&rseed);
		f.z += var_trans * ran2::gasdev(&rseed);
	}

	void RigidBody::addRandomTorque() {
		M.x += var_rot * ran2::gasdev(&rseed);
		M.y += var_rot * ran2::gasdev(&rseed);
		M.z += var_rot * ran2::gasdev(&rseed);
	}

	void RigidBody::integrate() {
		r += f * (timeStep / ksi_trans);
		q = quaternion::mul(quaternion::makeRotationAround(M * (timeStep / ksi_rot)), q);
	}

	void RigidBody::reset(float4 r_, float mass_, float R_) {
		r0 = r = r_;
		mass = mass_;
		R = R_;
		q = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		updateParameters();
	}
}
