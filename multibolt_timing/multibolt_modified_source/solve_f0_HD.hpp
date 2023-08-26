// SPDX-License-Identifier: MIT
//
// MultiBolt - a multi-term Boltzmann equation solver for low temperature plasma
// 
// Copyright 2021-2022 Max Flynn, Jacob Stephens
// 
// Licensed under the MIT License.
// A copy of the license should be distributed with this software (see LICENSE).
// The license text may also be found at: https://opensource.org/licenses/MIT
// --------------------------------------------------------------------------------- 

#pragma once

#include "multibolt"
#include <chrono>

// solve governing equation for f0 in HD conditions
// This is not the same equation as for SST
void mb::BoltzmannSolver::solve_f0_HD() {

	mb::debug_statement("Begin solving governing equation for f0 (HD).");

	omega0 = 0; // initial guess for omega0

	//arma::wall_clock local_timer;

	arma::sword N = p.Nu * p.N_terms;


	arma::Col<double> b(N, arma::fill::zeros);	// solve-against vector (contains normalization condition)
	x_f = arma::colvec(N, arma::fill::zeros);						// solution vector (contains distribution functions)
	arma::SpMat<double> A(N, N);				// coefficient matrix

	// initialize "previous values"
	double omega0_prev = 0; // previous value, tracking for convergence
	double avg_en_prev = 0;

	bool CONVERGED = false;

	for (iter_f0_HD = 1; iter_f0_HD <= p.iter_min || (iter_f0_HD < p.iter_max && !CONVERGED); ++iter_f0_HD) {

		auto t1 = std::chrono::high_resolution_clock::now();

		A.zeros();
		b.zeros();

		for (arma::sword ell = 0; ell < p.N_terms; ++ell) {
			gov_eq_f0_HD(ell, A);
		}
		A = A - A_scattering;

		// rescale to avoid very small numbers?
		double MAX_NUM = arma::max(arma::max(A));
		A = A / MAX_NUM;


		// Normalization condition
		A.row(0).zeros();
		A.head_cols(p.Nu).row(0) += sqrt(g.ue / mb::QE).t() * g.Du / mb::QE;
		b.at(0) = 1.0;


		//// Boundary conditions
		for (arma::sword ell = 1; ell < p.N_terms; ++ell) {
			boundary_conditions(ell, A, b);
		}

		auto t2 = std::chrono::high_resolution_clock::now();
		wtime_Abuild += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		//timer = arma::wall_clock(); // essentially, turn off the timer by blipping it


		avg_en_prev = calculate_avg_en(); // would use previous f0

		t1 = std::chrono::high_resolution_clock::now();
		bool solve_success = arma::spsolve(x_f, A, b, "superlu", multibolt_superlu_opts()); // solve
		t2 = std::chrono::high_resolution_clock::now();
		wtime_spsolve += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		N_spsolve++;
		//timer = arma::wall_clock(); // essentially, turn off the timer by blipping it


		if (!solve_success) {
			mb::error_statement("Solution failing: spsolve unsuccesful.");
			this->solution_succesful = false;
			return;
		}

		double avg_en = calculate_avg_en();

		//x_f(g.ell0_idx) = x_f(g.ell0_idx) / check_f0_normalization();

		omega0_prev = omega0;
		omega0 = calculate_k_iz_eff_N(); // by definition, omega0 = k_iz_eff_N
		omega0 = (1.0 - present_weight_f0) * omega0_prev + present_weight_f0 * omega0;



		// initial check: is the eedf wildly wrong?
		if (iter_f0_HD > p.iter_min && !NO_ITERATION) {
		
			double fnorm = check_f0_normalization();
			if (fnorm > 1.1 || fnorm < 0.9) { // if norm is disobeyed, solution might be unstable

				present_weight_f0 = present_weight_f0 / 2.0;

				mb::normal_statement("EEDF normalization fails, EEDF may be poorly-behaved.");
				
				if (present_weight_f0 > 1e-4) {

					mb::normal_statement("Attempting to save solution by reducing weight_f0.");


					present_weight_f0 = std::max(present_weight_f0, 1e-4);
					mb::normal_statement("weight_f0 is now: " + mb::mb_format(present_weight_f0) + ".");

					omega0_prev = 0;
					omega0 = 0;
					avg_en = 0;
					avg_en_prev = 0;

					iter_f0_HD = 1;

					continue;
				}
			}

		}




		// check if ionization is ocurring at all - converge against avg_en instead if need be
		if (this->NO_ITERATION || omega0 == 0) {
			mb::normal_statement("Single iteration case: no ionization or attachment found.");
			mb::display_iteration_banner("f0", iter_f0_HD, "avg_en [eV]", avg_en, 0);
			CONVERGED = true;
			break;
		}
		else {
			mb::display_iteration_banner("f0", iter_f0_HD, "avg_en [eV]", avg_en, err(avg_en, avg_en_prev));
			CONVERGED = mb::converged(p.conv_err, avg_en, avg_en_prev);

			//mb::display_iteration_banner("f0", iter_f0_HD, "omega0 [m**-3 s^-1]", omega0, err(omega0, omega0_prev));
			//CONVERGED = mb::converged(p.conv_err, omega0, omega0_prev);
		}

		// extra check: assume that if you are past iter_min, you're willing to let SHYness boot the solution out early
		if (iter_f0_HD > p.iter_min) {
			if (this->check_is_solution_failing(0) == true) {
				return;
			}
		}

	}

	mb::check_did_not_converge("f0", iter_f0_HD, p.iter_max);

	mb::debug_statement("Exit solving governing equation for f0 (HD).");
	

}


