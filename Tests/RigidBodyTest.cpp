#include <map>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "Updaters/RigidBody.cuh"

int runsCount = 10000;
int stepsCount = 10000;
float histStep = 0.1f;

namespace rigid_body {

struct Stats {
	float4 sum;
	float4 square_sum;
	int samples;
	typedef std::map<float, int> Hist;
	Hist x_hist;
	Hist y_hist;
	Hist z_hist;

	Stats()
		: sum(make_float4(0.0f, 0.0f, 0.0f, 0.0f))
		, square_sum(make_float4(0.0f, 0.0f, 0.0f, 0.0f))
		, samples(0)
	{}

	void addHist(std::map<float, int>& hist, float v) {
		hist[ceil((v - histStep/2.0f)/histStep)*histStep]++;
	}

	void add(float4 sample) {
		addHist(x_hist, sample.x);
		addHist(y_hist, sample.y);
		addHist(z_hist, sample.z);
		sum += sample;
		sample.x *= sample.x;
		sample.y *= sample.y;
		sample.z *= sample.z;
		square_sum += sample;
		samples++;
	}

	float4 getAverage() const {
		return sum * (1/(float)samples);
	}

	float4 getSquareAverage() const {
		return square_sum * (1/(float)samples);
	}

	float4 getDispersion() const {
		float4 avg = getAverage();
		avg.x *= avg.x;
		avg.y *= avg.y;
		avg.z *= avg.z;
		return getSquareAverage() - avg;
	}
};

struct AllStats {
	Stats position;
	Stats rotation;

	void add(float4 r, float4 q) {
		position.add(r);
		rotation.add(quaternion::toRotationVector(q));
	}
};

std::ostream& operator<<(std::ostream& o, const Stats::Hist& hist) {
	for (Stats::Hist::const_iterator i = hist.begin(), e = hist.end(); i != e; ++i)
		o << i->first << ";" << i->second << "\n";
	return o;
}

std::ostream& operator<<(std::ostream& o, const Stats& s) {
	return o << "E(x) = " << s.getAverage().x << "\t"
				"E(y) = " << s.getAverage().y << "\t"
				"E(z) = " << s.getAverage().z << "\t"
				"D(x) = " << s.getDispersion().x << "\t"
				"D(y) = " << s.getDispersion().y << "\t"
				"D(z) = " << s.getDispersion().z << "\n";
}

void saveToFile(const std::string& filename, const Stats::Hist& hist) {
	std::ofstream o(filename.c_str());
	o << hist;
}

void saveToFile(const std::string& prefix, int step, const Stats::Hist& hist) {
	std::stringstream ss;
	ss << prefix << "." << std::setw(7) << std::setfill('0') << step << ".txt";
	saveToFile(ss.str(), hist);
}

float theoreticalDispersionCoef(float T, float ksi) {
	float D = constant::Boltzmann(unitssystem::gromacs)*T/ksi;
	return 2.0f*D;
}

float theoreticalDispersion(float t, float T, float ksi) {
	return theoreticalDispersionCoef(T, ksi)*t;
}

void test() {
	// Global Parameters
	timeStep = 10;
	T = 300.0f;
	viscosity = calculateViscosity(T);
	rseed = time(0);

	// Resulting data storage
	typedef std::map<int, AllStats> Data;
	Data data;

	// Simulation
	RigidBody rb;
	int progress = 0;
	printf("Rigid Body Test Progress:   0%%\n");
	for (int run = 0; run < runsCount; run++) {
		rb.reset(make_float4(0.0f, 0.0f, 0.0f, 0.0f), 1.0f, 1.0f);
		for (int step = 0; step < stepsCount; step++) {
			if ((step & 0x0f) == 0) // each 16 steps
				data[step].add(rb.r, rb.q);
			rb.zeroForceAndTorque();
			rb.addRandomForce();
			rb.addRandomTorque();
			rb.integrate();
		}
		// Output progress
		int curProgress = int(100*(run+1)/runsCount);
		if (progress != curProgress) {
			progress = curProgress;
			printf("Rigid Body Test Progress: %3d%%\n", progress);
		}
	}
	printf("\n");

	{ // Analyze position
		std::cout << ">>>>>>> TRANSITION <<<<<<<" << std::endl;
		std::ofstream Ex("average_position_x.txt");
		std::ofstream Ey("average_position_y.txt");
		std::ofstream Ez("average_position_z.txt");
		std::ofstream Dx("dispersion_position_x.txt");
		std::ofstream Dy("dispersion_position_y.txt");
		std::ofstream Dz("dispersion_position_z.txt");
		std::ofstream thD("dispersion_position_theoretical.txt");
		for (Data::iterator i = data.begin(), e = data.end(); i != e; ++i) {
			std::cout << "step = " << i->first << "\t" << i->second.position;
			Ex << i->first << ';' << i->second.position.getAverage().x << '\n';
			Ey << i->first << ';' << i->second.position.getAverage().y << '\n';
			Ez << i->first << ';' << i->second.position.getAverage().z << '\n';
			Dx << i->first << ';' << i->second.position.getDispersion().x << '\n';
			Dy << i->first << ';' << i->second.position.getDispersion().y << '\n';
			Dz << i->first << ';' << i->second.position.getDispersion().z << '\n';
			thD << i->first << ';' << theoreticalDispersion(i->first*timeStep, T, rb.ksi_trans) << '\n';
			saveToFile("hist_position_x", i->first, i->second.position.x_hist);
			saveToFile("hist_position_y", i->first, i->second.position.y_hist);
			saveToFile("hist_position_z", i->first, i->second.position.z_hist);
		}
	}

	{ // Analyze rotation
		std::cout << ">>>>>>> ROTATION <<<<<<<" << std::endl;
		std::ofstream Ex("average_rotation_x.txt");
		std::ofstream Ey("average_rotation_y.txt");
		std::ofstream Ez("average_rotation_z.txt");
		std::ofstream Dx("dispersion_rotation_x.txt");
		std::ofstream Dy("dispersion_rotation_y.txt");
		std::ofstream Dz("dispersion_rotation_z.txt");
		std::ofstream thD("dispersion_rotation_theoretical.txt");
		for (Data::iterator i = data.begin(), e = data.end(); i != e; ++i) {
			std::cout << "step = " << i->first << "\t" << i->second.rotation;
			Ex << i->first << ';' << i->second.rotation.getAverage().x << '\n';
			Ey << i->first << ';' << i->second.rotation.getAverage().y << '\n';
			Ez << i->first << ';' << i->second.rotation.getAverage().z << '\n';
			Dx << i->first << ';' << i->second.rotation.getDispersion().x << '\n';
			Dy << i->first << ';' << i->second.rotation.getDispersion().y << '\n';
			Dz << i->first << ';' << i->second.rotation.getDispersion().z << '\n';
			thD << i->first << ';' << theoreticalDispersion(i->first*timeStep, T, rb.ksi_rot) << '\n';
			saveToFile("hist_rotation_x", i->first, i->second.rotation.x_hist);
			saveToFile("hist_rotation_y", i->first, i->second.rotation.y_hist);
			saveToFile("hist_rotation_z", i->first, i->second.rotation.z_hist);
		}
	}

	float th_k_trans = theoreticalDispersionCoef(T, rb.ksi_trans);
	float th_k_rot = theoreticalDispersionCoef(T, rb.ksi_rot);
	float4 ex_k_trans = data.rbegin()->second.position.getDispersion() * (1.0f/data.rbegin()->first/timeStep);
	float4 ex_k_rot = (++data.begin())->second.rotation.getDispersion() * (1.0f/(++data.begin())->first/timeStep);
	float4 err_trans = ex_k_trans*(100.0f/th_k_trans) - make_float4(100.0f, 100.0f, 100.0f, 0.0f);
	float4 err_rot = ex_k_rot*(100.0f/th_k_rot) - make_float4(100.0f, 100.0f, 100.0f, 0.0f);
	std::cout << ">>>>>>> COMPARE RESULTS <<<<<<" << std::endl;
	std::cout << "TRANSITION DISPERSION COEF:"
				 "\ttheoretical=" << th_k_trans <<
				 "\tx_sim=" << ex_k_trans.x <<
				 "\ty_sim=" << ex_k_trans.y <<
				 "\tz_sim=" << ex_k_trans.z <<
				 std::endl;
	std::cout << "ROTATION DISPERSION COEF:"
				 "\ttheoretical=" << th_k_rot <<
				 "\tx_sim=" << ex_k_rot.x <<
				 "\ty_sim=" << ex_k_rot.y <<
				 "\tz_sim=" << ex_k_rot.z <<
				 std::endl;
	std::cout << "TRANSITION ERROR [%]:"
				 "\tx_err=" << fabs(err_trans.x) <<
				 "\ty_err=" << fabs(err_trans.y) <<
				 "\tz_err=" << fabs(err_trans.z) <<
				 std::endl;
	std::cout << "ROTATION ERROR [%]:"
				 "\tx_err=" << fabs(err_rot.x) <<
				 "\ty_err=" << fabs(err_rot.y) <<
				 "\tz_err=" << fabs(err_rot.z) <<
				 std::endl;

}

} // namespace rigid_body

int main(int ac, char** av) {
	if (ac == 2 && std::string(av[1]) == "--help") {
		std::cerr << "usage: " << av[0] <<
					 " [runsCount, Def.Val: " << runsCount << "]"
					 " [stepsCount, Def.Val: " << stepsCount << "]"
					 " [histStep, Def.Val: " << histStep << "]"
				  << std::endl;
		exit(1);
	}
	if (ac > 1)
		runsCount = atoi(av[1]);
	if (ac > 2)
		stepsCount = atoi(av[2]);
	if (ac > 3)
		histStep = atof(av[3]);
	rigid_body::test();
}
