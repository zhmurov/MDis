#pragma once

#include "Cuda.h"
#include "LinearAlgebra.h"
#include "UnitsSystemDoc.h"

static const int QuantityDim = 6;

//! Dimention of a quantity
//! \note additive group is used instead of multiplicative one
//! for example: L^1 * T^(-2)  <==>  _L - 2*_T
typedef Vector<double, QuantityDim> Quantity;

namespace {

Quantity _0 = {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }}; //!< Non-dimentional magnitude

//! Basis of dimentional space
Quantity _L = {{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 }}; //!< Distance
Quantity _M = {{ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 }}; //!< Mass
Quantity _T = {{ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 }}; //!< Time
Quantity _I = {{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }}; //!< Current
Quantity _K = {{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }}; //!< Temperature
Quantity _N = {{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }}; //!< Amount of substance

Quantity _Energy = {{ 2.0, 1.0, -2.0, 0.0, 0.0, 0.0 }};
Quantity _Force  = {{ 1.0, 1.0, -2.0, 0.0, 0.0, 0.0 }};
Quantity _Torque = {{ 2.0, 1.0, -2.0, 0.0, 0.0, 0.0 }};

} // anonymous namespace

//! Logariphmic coefficients of units system transformation
//! Corresponding transformation of basic variables in equations have form:
//!   l = exp(u(0)) * l'  for length
//!   m = exp(u(1)) * m'  for mass
//!   t = exp(u(2)) * t'  for time
//!   i = exp(u(3)) * i'  for electric current 
//!   k = exp(u(4)) * k'  for temperature
//!   n = exp(u(5)) * n'  for amount of substance
//! where {l, m, t, i, k, n} - variables in old system
//! and {l', m', t', i', k', n'} - variables in new system.
//! Transformation of arbitrary dimentional value have form:
//!   Let [x] = L^lo * M^mo * T^to * I^io * K^ko * N^no
//!   Then qnt(x) = (lo, mo, to, io, ko, no)
//!
//!     x = exp(u * qnt(x)) * x'
//!
//! NOTE: In most cases old system is SI units system or other widely known system (e.g. CGS)
typedef Vector<double, QuantityDim> UnitsSystem;

namespace unitssystem {

//! Returns transformation coefficient:
//! v = coef * v', where [v] = L^lo * M^mo * T^to * I^io * K^ki * N^ni
	CUDA_FUNC inline float Coefficient(const UnitsSystem& system, const Quantity& qnt) {
		return (float)exp(-qnt * system);
	}

	CUDA_FUNC inline float InvCoefficient(const UnitsSystem& system, const Quantity& qnt) {
		return (float)exp(qnt * system);
	}

	//! Transforms given dimentional \a value
	CUDA_FUNC inline float Convert(const UnitsSystem& system, float value, const Quantity& qnt) {
		return (float)exp(-qnt * system) * value;
	}

	//! Transforms given dimentional \a value inversely
	CUDA_FUNC inline float InvConvert(const UnitsSystem& system, float value, const Quantity& qnt) {
		return (float)exp(qnt * system) * value;
	}

	//! Transforms given dimentional \a value from one system to another
	CUDA_FUNC inline float Convert(const UnitsSystem& from_system, const UnitsSystem& to_system, float value, const Quantity& qnt) {
		return (float)exp((from_system - to_system) * qnt) * value;
	}

	//! Determine unit system coefficients using unity values (i.e. equal to 1) for basic quantities (length, mass and etc)
	//! Arguments are values in basis units system
	CUDA_FUNC inline UnitsSystem MakeSystemWithUnity(
			UnitsSystem basis,
			double length_unit,
			double mass_unit,
			double time_unit,
			double current_unit,
			double temperature_unit,
			double amount_of_substance_unit)
	{
		UnitsSystem units = {{
				log(length_unit),
				log(mass_unit),
				log(time_unit),
				log(current_unit),
				log(temperature_unit),
				log(amount_of_substance_unit),
			}};
		units -= basis;
		return units;
	}

} // namespace unitssystem
