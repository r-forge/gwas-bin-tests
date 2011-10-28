/*
 * Matrix3D3d.hpp
 *
 *  Created on: Jun 21, 2010
 *      Author: kforner
 */

#ifndef Matrix3D3D_HPP_
#define Matrix3D3D_HPP_

#include "misc.hpp"

/*
 * a Matrix3D implementation of  c-like arrays that uses 2 additional arrays of pointers
 * so the usual access syntax [x][y][z] works without any overhead
 * , also the .c_contiguous_array()[x*ny*nz + y*nz + z] syntax can be used for direct access
 *
 * N.B: there is a memory overhead of O(nx*(ny+1))
 */

template<typename T>
class Matrix3D
{
public:
	// ========== LIFECYCLE ==========
	Matrix3D(int nx, int ny, int nz, T value = 0);

	// Based on the Law Of The Big Three:
	//	~Matrix3D();
	Matrix3D(Matrix3D const& source); // no need, default is fine
	Matrix3D<T>& operator = (Matrix3D<T> const& source);

public:
	// =========== ACCESSORS ==========

	int nx() const;
	int ny() const;
	int nz() const;
	int length() const;

	T* begin() { return c_contiguous_array(); }
	T const* begin() const { return c_contiguous_array(); }
	T* end() { return c_contiguous_array() + length();  }
	T const* end() const { return c_contiguous_array() + length(); }

	T* c_contiguous_array();
	T const* c_contiguous_array() const;

	T*** x_begin() { return &_x[0]; }
	T const*const*const* x_begin() const { return  &_x[0]; }
	T*** x_end() { return x_begin() + _nx; }
	T const*const* x_end() const { return x_begin() + _nx; }

	T** operator[](int x);
	T const* const*const operator[](int x) const;

	T const& operator()(int x, int y, int z) const;
	T& operator()(int x, int y, int z);

public:
	// =========== METHODS ==========

	bool hasSameDimensions(const Matrix3D<T>& m) const;

	bool operator == (const Matrix3D<T>& m) const;

	bool operator != (const Matrix3D<T>& m) const { return ! (*this == m); }

	void fill(T value);

//	template <class Function>
//	void fill(Function f);
//
//	template <class Function>
//	void apply(Function f);

	std::ostream& printOn(std::ostream& o) const;
	friend std::ostream& operator<<(std::ostream& o, Matrix3D<T> const& m) {
		return m.printOn(o);
	}

protected:
	void init();

private:
	// ========== INSTANCE DATA ==========
	vector<T> _actual_array;
	vector<T**> _x;
	vector<T*> _y;

	int _nx, _ny, _nz;
};

template<typename T> inline Matrix3D<T>::Matrix3D(int nx, int ny, int nz, T value) :
	_actual_array(nx * ny * nz, value), _x(nx), _y(nx*ny), _nx(nx),
			_ny(ny), _nz(nz) {
	assert( nx*ny*nz != 0);
	init();
}

template<typename T> inline Matrix3D<T>::Matrix3D(const Matrix3D<T>& s) :
	_actual_array(s._actual_array),
	_x(s._x.size()), _y(s._y.size())
	, _nx(s.nx()), _ny(s.ny()), _nz(s.nz())  {
	init();
}

template<typename T> inline void Matrix3D<T>::init() {
	_x.resize(_nx);
	_y.resize(_nx *_ny);

	T* current_y = &_actual_array[0]; // safe way to get pointer on underlying array
	// compute y pointers
	const int nb2 = _nx *_ny;
	for (int i = 0; i < nb2; ++i) {
		_y[i] = current_y;
		current_y += _nz;
	}

	// compute x pointers
	for (int i = 0; i < _nx; ++i) {
		_x[i] = &_y[i*_ny];
	}
}

template<typename T> inline int Matrix3D<T>::length() const {
	return _actual_array.size();
}


template<typename T> inline int Matrix3D<T>::nx() const {
	return _nx;
}

template<typename T> inline int Matrix3D<T>::ny() const {
	return _ny;
}

template<typename T> inline int Matrix3D<T>::nz() const {
	return _nz;
}



template<typename T> inline const T *Matrix3D<T>::c_contiguous_array() const {
	return _y[0];
}

template<typename T> inline T *Matrix3D<T>::c_contiguous_array() {
	return  _y[0];
}


template<typename T> inline const T & Matrix3D<T>::operator ()(int x, int y, int z) const {
	return _x[x][y][z];
}

template<typename T> inline T & Matrix3D<T>::operator ()(int x, int y, int z) {
	return _x[x][y][z];
}

template<typename T> inline T** Matrix3D<T>::operator [](int x) {
	return _x[x];
}

template<typename T> inline const T * const* const Matrix3D<T>::operator [](int x) const {
	return _x[x];
}

template<typename T>
inline
bool Matrix3D<T>::hasSameDimensions(const Matrix3D<T>& m) const
{
	return _nx == m.nx() && _ny == m.ny() && _nz == m.nz();
}

template<typename T>
inline
bool Matrix3D<T>::operator == (const Matrix3D<T>& m) const {
	return hasSameDimensions(m)	&& _actual_array == m._actual_array;
}

template<typename T> inline std::ostream & Matrix3D<T>::printOn(std::ostream & o) const {
	o << "Matrix3D(" << nx() << ", " << ny() << ", " << nz() << ")\n";

	for (int x = 0; x < _nx; ++x) {
		o << "x=" << x << ":\n";
		for(int y = 0; y < _ny; ++y) {
			for(int z = 0; z < _nz; ++z)
				o << (*this)(x,y,z) << "\t\t";
			o << "\n";
		}
	}

	return o;
}

template<typename T>
inline Matrix3D<T>&  Matrix3D<T>::operator = (Matrix3D<T> const& s)
{
	if ( &s != this) {
		_nx = s._nx;
		_ny = s._ny;
		_actual_array = s._actual_array;
		init();
	}
	return *this;
}


template<typename T>
void Matrix3D<T>::fill(T value)
{
	T* cursor = begin();
	T* the_end = end();
	while ( cursor < the_end)
		*cursor++ = value;

}



#endif /* Matrix3D3D_HPP_ */
