// SPDX-License-Identifier: MIT
//
// MultiBolt - a multi-term Boltzmann equation solver for low temperature plasma
// 
// Copyright 2021 Max Flynn, Jacob Stephens
// 
// Licensed under the MIT License.
// A copy of the license should be distributed with this software (see LICENSE).
// The license text may also be found at: https://opensource.org/licenses/MIT
// --------------------------------------------------------------------------------- 

// a minimal file for playing around in MultiBolt

#include <iostream>
#include <multibolt>


int main() {

	//const int N_xsec = 100;
	const int kmax = 17; // N vib
	const int N_spec = 59; // Changes to be number of de's
	const int N_trial = 55;
	double wtime_total[N_trial] = { };
	double wtime_spsolve[N_trial] = { };
	double wtime_grid[N_trial] = { };
	double wtime_Abuild[N_trial] = { };
	double ebar[N_trial] = { };
	int N_spsolve[N_trial] = { };
	int N_grid[N_trial] = { };

	double wtime_total2[kmax] = { };
	double wtime_spsolve2[kmax] = { };
	double wtime_grid2[kmax] = { };
	double wtime_Abuild2[kmax] = { };
	int N_spsolve2[kmax] = { };
	int N_grid2[kmax] = { };

	// imax is now the number of vibration species
	int imax[kmax] = {1, 1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59};

	double EN0 = 10.0;
	int Nu_grid = 631;

	double T_array[N_trial] = {};
	T_array[0] = 300.0;
	for (int i = 1; i < N_trial; ++i) {
		T_array[i] = T_array[i - 1] + 50.0;
	}

	// stdout mode
	//mb::mbSpeaker.printmode_normal_statements();
	//mb::mbSpeaker.printmode_debug_statements();
	mb::mbSpeaker.printmode_no_statements();

	// Species Data
	double de[N_spec] = { 0.000000e0, 2.88600e-1, 5.73700e-1, 8.55300e-1, \
						  1.133500e0, 1.408200e0, 1.679400e0, 1.947000e0, \
						  2.211100e0, 2.471700e0, 2.728700e0, 2.982100e0, \
						  3.232000e0, 3.478200e0, 3.720800e0, 3.720800e0, \
						  4.195100e0, 4.426800e0, 4.654700e0, 4.878900e0, \
						  5.099300e0, 5.315900e0, 5.528600e0, 5.737500e0, \
						  5.942300e0, 6.143200e0, 6.340000e0, 6.532700e0, \
						  6.721100e0, 6.905200e0, 7.085000e0, 7.260200e0, \
						  7.430900e0, 7.596800e0, 7.757900e0, 7.914000e0, \
						  8.065000e0, 8.210800e0, 8.351100e0, 8.485700e0, \
						  8.614600e0, 8.737400e0, 8.853900e0, 8.964000e0, \
						  9.067500e0, 9.163900e0, 9.253300e0, 9.335300e0, \
						  9.409800e0, 9.476700e0, 9.535900e0, 9.587500e0, \
						  9.631400e0, 9.667700e0, 9.696500e0, 9.718100e0, \
						  9.733300e0, 9.743300e0, 9.750100e0 };

	// BoltzmannParameters Definition //

	mb::BoltzmannParameters p = mb::BoltzmannParameters();

	// Set Parameters
	p.model = mb::ModelCode::HD;
	p.N_terms = 2;
	p.Nu = Nu_grid;
	p.iter_min = 3;		// Force continue if below this number of iterations
	p.iter_max = 70;	// Force break if you haven't converged within this many iterations to avoid wasting time
	p.weight_f0 = 1.0;	// Iteration weight. N2 and Ar behave well, 1.0 will suffice.
	p.conv_err = 1.0e-8;	// relative error
	p.USE_ENERGY_REMAP = false; // Look for a new grid if necessary
	p.remap_grid_trial_max = 60;				// try a new grid no more than this many times
	p.remap_target_order_span = 10.0;			// how many orders would you prefer your EEDF to span (in eV^-3/2)
	p.USE_EV_MAX_GUESS = false;
	p.initial_eV_max = 100;
	p.p_Torr = 7.60000322e+02;
	p.DONT_ENFORCE_SUM = false;

	// Create Base Library
	std::string fname = { "../../cross-sections-unix/Laporta_N2_vib_set_new2.txt" };
	lib::Library lib0 = lib::Library(fname);
	lib0.assign_scattering_by_type(lib::isotropic_scattering(), lib::CollisionCode::elastic);
	lib0.assign_scattering_by_type(lib::isotropic_scattering(), lib::CollisionCode::excitation);
	lib0.assign_scattering_by_type(lib::isotropic_scattering(), lib::CollisionCode::ionization);
	lib0.assign_scattering_by_type(lib::isotropic_scattering(), lib::CollisionCode::superelastic);

	// Loop Maximum sets
	for (int k = 0; k < kmax; ++k) {

		// Set Number of retained vibrational states
		int ispec = imax[k];

		// copy library
		lib::Library lib = lib0;

		// Erase Species - Add Cross Sections C:\Users\lynch\source\repos\MultiBolt\cross-sections\Biagi_Kr.txt ""
		for (int i = lib0.allspecies.size()-1; i >= ispec; --i) {
			lib.erase_by_index(i);
		}

		// Loop Temperature Trials
		for (int it = 0; it < N_trial; ++it) {

			// Update Parameters
			p.T_K = T_array[it];
			p.EN_Td = EN0 * T_array[it] / 300.0;

			// Get Normalized Fractions of root species
			double sfrac[100] = { };
			double ssum = 0.0;
			for (int i = 0; i < ispec; ++i) {
				sfrac[i] = exp(-de[i] / (T_array[it] * 8.617333262145179e-05));
				ssum += sfrac[i];
			}
			for (int i = 0; i < ispec; ++i) {
				sfrac[i] /= ssum;
				std::string sname = { "N2_v" + std::to_string(i) };
				lib.assign_fracs_by_names(sname, sfrac[i]);
				//mb::extra_normal_statement(std::to_string(i) + ",  " + sname + ",  " + mb::mb_format(sfrac[i]) );
			}

			// Run Solver
			mb::BoltzmannSolver solver = mb::BoltzmannSolver(p, lib);

			// Store Times
			N_spsolve[it] = solver.N_spsolve;
			N_grid[it] = solver.N_grid;
			wtime_total[it] = solver.wtime_total;
			wtime_spsolve[it] = solver.wtime_spsolve;
			wtime_grid[it] = solver.wtime_grid;
			wtime_Abuild[it] = solver.wtime_Abuild;
			ebar[it] = solver.out.avg_en;
		}

		// Print Times in 28 trials
		for (int it = 0; it < N_trial; ++it) {
			mb::extra_normal_statement(std::to_string(k + 1) + ", " \
				+ std::to_string(it + 1) + ", " \
				+ std::to_string(N_spsolve[it]) + ", " \
				+ std::to_string(N_grid[it]) + ", " \
				+ mb::mb_format(wtime_total[it]) + ", " \
				+ mb::mb_format(wtime_spsolve[it]) + ", " \
				+ mb::mb_format(wtime_grid[it]) + ", " \
				+ mb::mb_format(wtime_Abuild[it]) + ", " \
				+ mb::mb_format(ebar[it]));
			wtime_total2[k] += wtime_total[it];
			wtime_spsolve2[k] += wtime_spsolve[it];
			wtime_grid2[k] += wtime_grid[it];
			wtime_Abuild2[k] += wtime_Abuild[it];
			N_spsolve2[k] += N_spsolve[it];
			N_grid2[k] += N_grid[it];
		}

	}

	mb::extra_normal_statement("SUMMARY:");
	mb::extra_normal_statement(std::to_string(N_trial));
	for (int k = 0; k < kmax; ++k) {
		mb::extra_normal_statement(std::to_string(k + 1) + ", " \
			+ std::to_string(N_spsolve2[k]) + ", " \
			+ std::to_string(N_grid2[k]) + ", " \
			+ mb::mb_format(wtime_total2[k]) + ", " \
			+ mb::mb_format(wtime_spsolve2[k]) + ", " \
			+ mb::mb_format(wtime_grid2[k]) + ", " \
			+ mb::mb_format(wtime_Abuild2[k]) + ", " \
			+ mb::mb_format((wtime_total2[k] - wtime_spsolve2[k] \
				- wtime_grid2[k] - wtime_Abuild2[k]) / wtime_total2[k])
		);
	}

	return 0;

}
