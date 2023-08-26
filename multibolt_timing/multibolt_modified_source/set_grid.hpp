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


// Set the grids on every cross section to that of ue and uo
void mb::BoltzmannSolver::set_grid() {

	//arma::wall_clock local_timer;
	auto t1 = std::chrono::high_resolution_clock::now();

	this->g = mb::Grid(); // reset new grid

	g.Du = (present_eV_max / (p.Nu - 1)) * mb::QE;

	arma::vec temp_ue = arma::linspace(0.5 * g.Du, g.Du * (double(p.Nu) - 0.5), p.Nu);
	arma::vec temp_uo = arma::linspace(g.Du, p.Nu * g.Du, p.Nu);
	g.ue = arma::colvec(p.Nu, arma::fill::zeros);
	g.uo = arma::colvec(p.Nu, arma::fill::zeros);

	// had issues setting the vector arma-style, dummy fix
	for (int i = 0; i < p.Nu; ++i) {
		g.ue.at(i) = temp_ue.at(i);
		g.uo.at(i) = temp_uo.at(i);
	}



	g.ell0_idx = arma::regspace<arma::uvec>(0, 1, p.Nu - 1); // 0, 1, 2.....Nu-1 -> corresponds to ell=0 indices
	g.ell1_idx = g.ell0_idx + p.Nu; // Nu, Nu+1, Nu+2.....2*Nu-1 -> corresponds to ell=1 indices

	g.t_idx = arma::regspace<arma::uvec>(1, 1, p.Nu - 1); // truncated once from zero // useful for (k, k-1)
	g.idx_t = arma::regspace<arma::uvec>(0, 1, p.Nu - 2);// truncated once from u_max //  useful for (k+1, k)
	g.t_idx_t = arma::regspace<arma::uvec>(1, 1, p.Nu - 2);


	mb::debug_statement("Applying even and odd grids (linear-spaced, per-energy).");

	for (auto& spec : lib.allspecies) {

		for (auto& x : spec->allcollisions) {
			x->set_grid(g.EVEN, g.ue / mb::QE);
			x->set_grid(g.ODD, g.uo / mb::QE);
		}
	}

	g.uo = g.ue;


	bool any_has_nan = false;
	// check for nan entries in grids
	// if this error flips, it is likely a development error
	for (auto& spec : lib.allspecies) {
		for (auto& x : spec->allcollisions) {

			if (x->gridded_s[g.EVEN].has_nan() || x->gridded_s[g.ODD].has_nan()) {
				mb::normal_statement("The following xsec (.print()) is gridded with nan:");
				x->print();
				any_has_nan = true;
			}
		}
	}
	if (any_has_nan) {
		mb::throw_err_statement("Solution cannot continue: Gridding procedure failed such that nan exists in at least one Xsec.");
	}



	// Based on this grid, figure out if iterations are actually needed
	// in order to not waste time iterating

	double accum_iz = 0;
	double accum_att = 0;

	for (auto& spec : lib.allspecies) {
		for (auto& x : spec->iz) {

			accum_iz += arma::accu(x->gridded_s[g.EVEN]);
		}
		for (auto& x : spec->att) {
			accum_att += arma::accu(x->gridded_s[g.EVEN]);
		}
	}

	if ((doubles_are_same(accum_att, 0) && doubles_are_same(accum_iz, 0))) {
		NO_ITERATION = true;
		mb::debug_statement("Detected that grid is all-conservative processes: all calculations will be single-iteration.");
	}
	else {
		NO_ITERATION = false; // ***
	}

	



	mb::normal_statement("Begin loading collision matrix.");
	
	// pre-load the scattering matrix.

	this->A_scattering = arma::SpMat<double>(p.Nu * p.N_terms, p.Nu * p.N_terms); // set to zeros
	for (arma::sword ell = 0; ell < p.N_terms; ++ell) {
		full_collision_operator(ell, this->A_scattering);
	}

	mb::normal_statement("Finishing loading collision matrix.");

	auto t2 = std::chrono::high_resolution_clock::now();
	wtime_grid += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	//timer = arma::wall_clock(); // essentially, turn off the timer by blipping it
	N_grid++;

};


