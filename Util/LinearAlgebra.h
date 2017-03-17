#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector_types.h>
#include <vector_functions.h>
#include "Cuda.h"

//
// Generic vector in eucledean space with dimentionality equals Dim
//
template <class T, unsigned Dim>
struct Vector {
  T value[Dim];

  CUDA_FUNC T operator()(unsigned i) const { return value[i]; }
  CUDA_FUNC T& operator()(unsigned i) { return value[i]; }
};

// Vector sum
template <class T, unsigned Dim>
CUDA_FUNC inline Vector<T, Dim> operator + (const Vector<T, Dim>& lh, const Vector<T, Dim>& rh) {
  Vector<T, Dim> ret(lh);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) += rh(i);
  return ret;
}

template <class T, unsigned Dim>
CUDA_FUNC inline Vector<T, Dim>& operator += (Vector<T, Dim>& lh, const Vector<T, Dim>& rh) {
  for (unsigned i = 0; i < Dim; i++)
    lh(i) += rh(i);
  return lh;
}

template <class T, unsigned Dim>
CUDA_FUNC inline Vector<T, Dim> operator - (const Vector<T, Dim>& lh, const Vector<T, Dim>& rh) {
  Vector<T, Dim> ret(lh);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) -= rh(i);
  return ret;
}

template <class T, unsigned Dim>
CUDA_FUNC inline Vector<T, Dim>& operator -= (Vector<T, Dim>& lh, const Vector<T, Dim>& rh) {
  for (unsigned i = 0; i < Dim; i++)
    lh(i) -= rh(i);
  return lh;
}

template <class T, unsigned Dim>
CUDA_FUNC inline Vector<T, Dim> operator - (const Vector<T, Dim>& v) {
  Vector<T, Dim> ret(v);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) = -ret(i);
  return ret;
}


// Scalar-Vector product
template <class T, unsigned Dim, class A>
CUDA_FUNC inline Vector<T, Dim> operator * (A lh, const Vector<T, Dim>& rh) {
  Vector<T, Dim> ret(rh);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) *= lh;
  return ret;
}

template <class T, unsigned Dim, class A>
CUDA_FUNC inline Vector<T, Dim> operator * (const Vector<T, Dim>& lh, A rh) {
  Vector<T, Dim> ret(lh);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) *= rh;
  return ret;
}

template <class T, unsigned Dim, class A>
CUDA_FUNC inline Vector<T, Dim>& operator *= (Vector<T, Dim>& lh, A rh) {
  for (unsigned i = 0; i < Dim; i++)
    lh(i) *= rh;
  return lh;
}

template <class T, unsigned Dim, class A>
CUDA_FUNC inline Vector<T, Dim> operator / (const Vector<T, Dim>& lh, A rh) {
  Vector<T, Dim> ret(lh);
  for (unsigned i = 0; i < Dim; i++)
    ret(i) /= rh;
  return ret;
}

template <class T, unsigned Dim, class A>
CUDA_FUNC inline Vector<T, Dim>& operator /= (Vector<T, Dim>& lh, A rh) {
  for (unsigned i = 0; i < Dim; i++)
    lh(i) /= rh;
  return lh;
}

namespace impl {

template <class T, unsigned Dim>
struct DistanceSqrdImpl {
  CUDA_FUNC static T Apply(const T* x, const T* y) {
    return (*x - *y) * (*x - *y) + DistanceSqrdImpl<T, Dim - 1>::Apply(x + 1, y + 1);
  }
};

template <class T>
struct DistanceSqrdImpl<T, 1> {
  CUDA_FUNC static T Apply(const T* x, const T* y) {
    return (*x - *y) * (*x - *y);
  }
};

} // namespace impl

// Distance
template <unsigned Dim, class T>
CUDA_FUNC inline T DistanceSqrd(const Vector<T, Dim>& x, const Vector<T, Dim>& y) {
  return impl::DistanceSqrdImpl<T, Dim>::Apply(x.value, y.value);
}

template <unsigned Dim, class T>
CUDA_FUNC inline T Distance(const Vector<T, Dim>& x, const Vector<T, Dim>& y) {
  return sqrt(DistanceSqrd(x, y));
}

// Vector length
template <class T, unsigned Dim>
CUDA_FUNC inline T Abs(const Vector<T, Dim>& v) {
  T ret(0);
  for (unsigned i = 0; i < Dim; i++)
    ret += v(i) * v(i);
  return sqrt(ret);
}

//
// Matrix
//
template <class T, unsigned Rows, unsigned Cols>
struct Matrix {
  enum { rows = Rows, cols = Cols, size = Rows * Cols };
  T value[size];

  CUDA_FUNC T operator()(unsigned r, unsigned c) const { return value[r * cols + c]; }
  CUDA_FUNC T& operator()(unsigned r, unsigned c) { return value[r * cols + c]; }

  CUDA_FUNC Vector<T, cols>& row(unsigned r) { return reinterpret_cast<Vector<T, cols>*>(value)[r]; }
  CUDA_FUNC Vector<T, rows> col(unsigned c) const {
    Vector<T, rows> ret;
    for (unsigned j = 0, i = c; j < rows; j++, i += rows)
      ret(j) = value[i];
    return ret;
 }
};

// Matrix sum
template <class T, unsigned R, unsigned C>
CUDA_FUNC inline Matrix<T, R, C> operator + (const Matrix<T, R, C>& lh, const Matrix<T, R, C>& rh) {
  Matrix<T, R, C> ret = lh;
  for (unsigned i = 0; i < lh.size; i++) {
    ret.value[i] += rh.value[i];
  }
  return ret;
}

template <class T, unsigned R, unsigned C>
CUDA_FUNC inline Matrix<T, R, C> operator - (const Matrix<T, R, C>& lh, const Matrix<T, R, C>& rh) {
  Matrix<T, R, C> ret = lh;
  for (unsigned i = 0; i < lh.size; i++) {
    ret.value[i] -= rh.value[i];
  }
  return ret;
}

// Matrix norms
template <class T, unsigned R, unsigned C>
CUDA_FUNC T NormL1(const Matrix<T, R, C>& m) {
  T norm = T();
  for (unsigned i = 0; i < m.size; i++)
    norm += m.value[i];
  return norm;
}

template <class T, unsigned R, unsigned C>
CUDA_FUNC T NormLInf(const Matrix<T, R, C>& m) {
  T norm = T();
  for (unsigned i = 0; i < m.size; i++) {
    T v = m.value[i];
    if (v < 0)
      v = -v;
    if (norm < v)
      norm = v;
  }
  return norm;
}

// Matrix-Vector product
template <class T, unsigned Rows, unsigned Cols>
CUDA_FUNC inline Vector<T, Rows> operator * (const Matrix<T, Rows, Cols>& lh, const Vector<T, Cols>& rh) {
  Vector<T, Rows> ret;
  for (unsigned i = 0; i < Rows; i++) {
    ret(i) = 0;
    for (unsigned j = 0; j < Cols; j++)
      ret(i) += lh(i, j) * rh(j);
  }
  return ret;
}

// Matrix transpose
template <class T, unsigned Rows, unsigned Cols>
CUDA_FUNC inline Matrix<T, Cols, Rows> operator * (const Matrix<T, Rows, Cols>& m) {
  Matrix<T, Cols, Rows> ret;
  for (unsigned i = 0; i < Rows; i++)
    for (unsigned j = 0; j < Cols; j++)
      ret(i, j) = m(j, i);
  return ret;
}

// Transpose(Matrix)-Vector product
template <class T, unsigned Rows, unsigned Cols>
CUDA_FUNC inline Vector<T, Rows> TransMul(const Matrix<T, Cols, Rows>& lh, const Vector<T, Cols>& rh) {
  Vector<T, Rows> ret;
  for (unsigned i = 0; i < Rows; i++) {
    ret(i) = 0;
    for (unsigned j = 0; j < Cols; j++)
      ret(i) += lh(j, i) * rh(j);
  }
  return ret;
}

// Matrix-Matrix product
template <class T, unsigned R1, unsigned C1R2, unsigned C2>
CUDA_FUNC inline Matrix<T, R1, C2> operator * (const Matrix<T, R1, C1R2>& lh, const Matrix<T, C1R2, C2>& rh) {
  Matrix<T, R1, C2> ret;
  for (unsigned r = 0; r < R1; r++) {
    for (unsigned c = 0; c < C2; c++) {
      ret(r, c) = 0;
      for (unsigned i = 0; i < C1R2; i++)
        ret(r, c) += lh(r, i) * rh(i, c);
    }
  }
  return ret;
}

// Transpose(Matrix)-Matrix product
template <class T, unsigned R1, unsigned C1R2, unsigned C2>
CUDA_FUNC inline Matrix<T, R1, C2> TransMul(const Matrix<T, C1R2, R1>& lh, const Matrix<T, C1R2, C2>& rh) {
  Matrix<T, R1, C2> ret;
  for (unsigned r = 0; r < R1; r++) {
    for (unsigned c = 0; c < C2; c++) {
      ret(r, c) = 0;
      for (unsigned i = 0; i < C1R2; i++)
        ret(r, c) += lh(i, r) * rh(i, c);
    }
  }
  return ret;
}

// Scalar product
template <class T, unsigned Dim>
CUDA_FUNC inline T operator * (const Vector<T, Dim>& lh, const Vector<T, Dim>& rh) {
  T ret = T();
  for (unsigned i = 0; i < Dim; i++)
    ret += lh(i)*rh(i);
  return ret;
}

// Vector product in 3D
template <class T>
CUDA_FUNC inline Vector<T, 3> VectorProduct(const Vector<T, 3>& lh, const Vector<T, 3>& rh) {
  Vector<T, 3> ret = {{lh(1)*rh(2) - lh(2)*rh(1), lh(2)*rh(0) - lh(0)*rh(2), lh(0)*rh(1) - lh(1)*rh(0)}};
  return ret;
}

// 3x3 Matrix inversion
template <class T>
CUDA_FUNC inline Matrix<T, 3, 3> Inv(const Matrix<T, 3, 3>& m) {
  // FIXME: this invertion algorithm is very slow
  Matrix<T, 3, 3> ret;
  Vector<T, 3> x0(m.col(0)), x1(m.col(1)), x2(m.col(2));
  T det = x0 * VectorProduct(x1, x2);
  //assert(det != 0) // how to handle it on device?
  
  T invdet = 1 / det;
  ret.row(0) = invdet * VectorProduct(x1, x2);
  ret.row(1) = invdet * VectorProduct(x2, x0);
  ret.row(2) = invdet * VectorProduct(x0, x1);
  return ret;
}

//
// Zero filling
//
template <class T, unsigned Dim>
CUDA_FUNC inline void Zero(Vector<T, Dim>& v) {
  for (unsigned i = 0; i < Dim; i++)
    v(i) = T();
}

template <class T, unsigned Rows, unsigned Cols>
CUDA_FUNC inline void Zero(Matrix<T, Rows, Cols>& m) {
  for (unsigned r = 0; r < Rows; r++)
    for (unsigned c = 0; c < Cols; c++)
      m(r, c) = T();
}

//
// Helpers
//
template <class T>
CUDA_FUNC inline Vector<T, 3> MakeVector3(T x, T y, T z) {
  Vector<T, 3> ret = {{x, y, z}};
  return ret;
}

template <class T, unsigned Dim, unsigned Index>
CUDA_FUNC inline Vector<T, Dim> BasisVectorImpl() {
  Vector<T, 3> ret;
  Zero(ret);
  ret(Index) = T(1);
  return ret;
}

template <class T, unsigned Dim, unsigned Index>
inline const Vector<T, Dim>& BasisVector() {
  // NOTE: this is not thread-safe
  static Vector<T, 3> ret = BasisVectorImpl<T, Dim, Index>();
  return ret;
}

template <class T>
CUDA_FUNC inline Matrix<T, 3, 3> MakeMatrix3(T x1, T x2, T x3, T y1, T y2, T y3, T z1, T z2, T z3) {
  Matrix<T, 3, 3> ret = {{x1, x2, x3, y1, y2, y3, z1, z2, z3}};
  return ret;
}

template <unsigned Dim>
CUDA_FUNC inline Vector<float, Dim>& AsVector(float4& x) {
  cuda_assert(Dim <= 4);
  return *reinterpret_cast<Vector<float, Dim>*>(&x);
}

template <unsigned Dim>
CUDA_FUNC inline const Vector<float, Dim>& AsVector(const float4& x) {
  cuda_assert(Dim <= 4);
  return *reinterpret_cast<const Vector<float, Dim>*>(&x);
}

template <unsigned Dim>
CUDA_FUNC inline Vector<float, Dim>& AsVector(float3& x) {
  cuda_assert(Dim <= 3);
  return *reinterpret_cast<Vector<float, Dim>*>(&x);
}

template <unsigned Dim>
CUDA_FUNC inline const Vector<float, Dim>& AsVector(const float3& x) {
  cuda_assert(Dim <= 3);
  return *reinterpret_cast<const Vector<float, Dim>*>(&x);
}

template <class T, unsigned Dim>
CUDA_FUNC inline Matrix<T, Dim, Dim> IdentityMatrix() {
  Matrix<T, Dim, Dim> ret;
  for (unsigned r = 0; r < Dim; r++)
    for (unsigned c = 0; c < Dim; c++)
      ret(r, c) = T(r == c? 1: 0);
  return ret;
}

template <class T>
CUDA_FUNC inline Matrix<T, 3, 3> RotationMatrix(Vector<T, 3> raxis) {
  T angle = Abs(raxis);
  if (fabs(angle) > 1e-6f) {
    Vector<T, 3> naxis = raxis / angle;
    T X = naxis(0), Y = naxis(1), Z = naxis(2);
    T cosA = cos(angle), sinA = sin(angle);
    Matrix<T, 3, 3> ret = {{
      cosA + (1-cosA)*X*X,   (1-cosA)*X*Y - sinA*Z, (1-cosA)*X*Z + sinA*Y,
      (1-cosA)*Y*X + sinA*Z, cosA + (1-cosA)*Y*Y,   (1-cosA)*Y*Z - sinA*X,
      (1-cosA)*Z*X - sinA*Y, (1-cosA)*Z*Y + sinA*X, cosA + (1-cosA)*Z*Z
    }};
    return ret;
  }
  return IdentityMatrix<T, 3>();
}

//
// Output
//
template <class T, unsigned D>
std::ostream& operator << (std::ostream& os, Vector<T, D> v) {
  os << "[";
  for (unsigned i = 0; i < D; i++)
    os << " " << v(i) << " ";
  os << "]";
  return os;
}

template <class T, unsigned R, unsigned C>
std::ostream& operator << (std::ostream& os, Matrix<T, R, C> m) {
  os << "[";
  for (unsigned r = 0; r < R; r++) {
    os << " [";
    for (unsigned c = 0; c < C; c++)
      os << " " << m(r, c) << " ";
    os << "] ";
  }
  os << "]";
  return os;
}
