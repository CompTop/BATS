#pragma once
/*
utilities for S^n and RP^n
*/

#include "data.hpp"
#include "data_gen.hpp"
#include <cmath>

namespace bats {

// normalizes columns of data to have norm 1
template <typename T>
void normalize_entries(DataSet<T>& data) {
	for (size_t j = 0; j < data.size(); j++) {
		T vnorm = norm(data[j]);
		data[j] /= vnorm;
	}
	return;
}

// force repulsion on vectors on sphere using RP distance
// goal is to get good coverage of projective plane
// v - d x n array of vectors on unit sphere (assume norm 1)
// step_max - limits step length
// repel points using inverse square law
template <typename T>
void force_repel_rp(
	DataSet<T> &v,
	T step_max
) {
	size_t d = v.dim();
	size_t n = v.size();

	// preallocate array for force vectors
	Matrix<T> dv(n,d);
	fill_zeros(dv);
	Matrix<T> f(1, d);
	fill_zeros(f);

	// step 1: compute force
	for (size_t j = 0; j < n; j++) {
		for (size_t i = 0; i < j; i++) {
			T normf = T(0);
			bool same_hemisphere = true;
			for (size_t k = 0; k < d; k++) {
				f(k) = v(j,k) - v(i,k);
				normf += f(k)*f(k);
			}
			// check to see if representative was from far side of sphere
			if (normf > d) {
				normf = T(0);
				same_hemisphere = false;
				for (size_t k = 0; k < d; k++) {
					f(k) = v(j,k) + v(i,k);
					normf += f(k)*f(k);
				}
			}

			normf = std::pow(normf, T(1.5)); // for inverse square law
			// add to total forces on i, j
			if (same_hemisphere) {
				#pragma omp simd
				for (size_t k = 0; k < d; k++) {
					T ck = f(k) / normf;
					dv(j,k) += ck;
					dv(i,k) -= ck;
				}
            } else {
				#pragma omp simd
				for (size_t k =0; k < d; k++) {
					T ck = f(k) / normf;
					dv(j,k) += ck;
					dv(i,k) += ck;
				}
			}
		}
	}

	// project forces to be tangent to sphere
	// enforce maximum step size
	// take step
	// project onto unit ball
	for (size_t j = 0; j < n; j++) {
		T ip = T(0);
		#pragma omp simd
		for (size_t k = 0; k < d; k++) {
			ip += v(j,k) * dv(j,k);
		}

		T dvn = T(0); // dv norm
		#pragma omp simd
		for (size_t k = 0; k < d; k++) {
			dv(j,k) -= ip*v(j,k); // project off normal component
			dvn += dv(j,k) * dv(j,k);
		}

		// enforce maximum step size
		if (dvn > step_max * step_max) {
			T cproj = step_max / std::sqrt(dvn);
			// project dv to be size step_max
			#pragma omp simd
			for (size_t k = 0; k < d; k++) {
				dv(j,k) *= cproj; // project
			}
		}

		// update position
		#pragma omp simd
		for (size_t k = 0; k < d; k++) {
			v(j,k) += dv(j,k);
		}

		// project onto ball
		T vnorm = norm(v[j]);
		v[j] /= vnorm;
	}

	// clean up temporary arrays
	dv.free();
	f.free();

	return;
}

} // namespace bats
