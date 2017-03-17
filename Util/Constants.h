#pragma once

#include "UnitsSystem.h"

namespace unitssystem {
static UnitsSystem si = {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
}

#define DEFINE_CONSTANT(name, sys, value, qnt) \
  namespace constant { \
    inline float name(const UnitsSystem& units) { \
      return unitssystem::Convert(unitssystem::sys, units, float(value), (qnt)); \
    } \
    inline double Log ## name(const UnitsSystem& units) { \
      return (unitssystem::sys - units) * (qnt) + log(value); \
    } \
  } \
  namespace quantity { \
    inline Quantity name() { return (qnt); } \
  } \
  /**/

DEFINE_CONSTANT(Electric, si, 8.854187817620e-12, -3*_L-_M+4*_T+2*_I)
DEFINE_CONSTANT(Magnetic, si, 1.2566370614e-6, _Force - 2*_I)
DEFINE_CONSTANT(LightSpeed, si, 299792456.2, _L - _T)
DEFINE_CONSTANT(ElectronCharge, si, 1.602176487e-19, _T+_I)
DEFINE_CONSTANT(Dalton, si, 1.6605402e-27, _M)
DEFINE_CONSTANT(Boltzmann, si, 1.38065042424e-23, _Energy - _K)
DEFINE_CONSTANT(Avogadro, si, 6.0221417930e23, -_N)

namespace unitssystem {
static UnitsSystem esu = MakeSystemWithUnity(si, 1e-2, 1e-3, 1.0, 0.1/constant::LightSpeed(si), 1.0, 1.0);
static UnitsSystem gromacs = MakeSystemWithUnity(si, 1e-9, constant::Dalton(si), 1e-12, 1e12*constant::ElectronCharge(si), 1.0f, 1.0f);
}
