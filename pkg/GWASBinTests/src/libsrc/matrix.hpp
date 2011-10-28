/*
 * matrix.hpp
 *
 *  Created on: Jun 18, 2010
 *      Author: kforner
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "misc.hpp"

/*
 * a Matrix implementation that uses row order major c-like arrays but with
 * pointer on each row to avoid arithmetic operation to access elements
 * so the usual access syntax [i][j] works without any overhead
 *  , also the [row*ncols + col] syntax can be used for direct access
 *
 * N.B: there is a memory overhead of O(nrows)

 */

template<typename T>
class Matrix
{
public:
	// ========== LIFECYCLE ==========
	Matrix(int nrows, int ncols, T value = 0);
	Matrix(const T*const* m, int nrows, int ncols);
	template<int NBCOLS> Matrix(T m[][NBCOLS], int nrows);


	// Based on the Law Of The Big Three:
	//	~Matrix();
	Matrix(Matrix const& source); // no need, default is fine
	Matrix<T>& operator = (Matrix<T> const& source);

public:
	// =========== ACCESSORS ==========

	int nb_rows() const;
	int nb_cols() const;

	T* row_by_row_begin() { return c_contiguous_array(); }
	T const* row_by_row_begin() const { return c_contiguous_array(); }
	T* row_by_row_end() { return c_contiguous_array() + _nb_rows*_nb_cols; }
	T const* row_by_row_end() const { return c_contiguous_array() + _nb_rows*_nb_cols; }

	T** c_arrays();
	T const* const * c_arrays() const;
	T** operator()();
	T const* const * operator()() const;

	T* begin();
	T const* begin() const;
	T* end();
	T const* end() const;
	T* c_contiguous_array();
	T const* c_contiguous_array() const;

	T** row_begin() { return &_rows[0]; }
	T const*const* row_begin() const { return &_rows[0]; }
	T** row_end() { return &_rows[0] + _nb_rows; }
	T const*const* row_end() const { return &_rows[0] + _nb_rows; }

	T* operator[](int row);
	T const* const operator[](int row) const;

	T const& operator()(int i, int j) const;
	T& operator()(int i, int j);

public:
	// =========== METHODS ==========

	bool operator == (const Matrix<T>& m) const;

	bool operator != (const Matrix<T>& m) const { return ! (*this == m); }

	void fill(T value);

	template <class Function>
	void fill(Function f);

	template <class Function>
	void apply(Function f);

	std::ostream& printOn(std::ostream& o) const;
	friend std::ostream& operator<<(std::ostream& o, Matrix<T> const& m) {
		return m.printOn(o);
	}

protected:
	void init();

private:
	// ========== INSTANCE DATA ==========
	vector<T> _actual_array;
	vector<T*> _rows;
	int _nb_rows;
	int _nb_cols;

};

template<typename T> inline Matrix<T>::Matrix(int nrows, int ncols, T value) :
	_actual_array(nrows * ncols, value), _rows(nrows), _nb_rows(nrows),
			_nb_cols(ncols) {
	init();
}

template<typename T> inline Matrix<T>::Matrix(const T*const* m,int nrows, int ncols) :
	_actual_array(nrows * ncols), _rows(nrows), _nb_rows(nrows),
			_nb_cols(ncols) {
	init();
	for (int i = 0; i < nrows; ++i)
		for (int j = 0; j < ncols; ++j)
			_rows[i][j] = m[i][j];
}

template<typename T>
template<int NBCOLS>  inline Matrix<T>::Matrix (T m[][NBCOLS], int nrows)
	: _actual_array(), _rows(),  _nb_rows(nrows), _nb_cols(NBCOLS)
{
	_actual_array.resize(_nb_rows * _nb_cols );
	_rows.resize(_nb_rows);
	init();
	for (int i = 0; i < _nb_rows; ++i)
		for (int j = 0; j < _nb_cols; ++j) {
			_rows[i][j] = m[i][j];
		}

}

template<typename T> inline Matrix<T>::Matrix(const Matrix<T>& source) :
	_actual_array(source._actual_array), _rows(source.nb_rows()), _nb_rows(source.nb_rows()),
			_nb_cols(source.nb_cols()) {
	init();
}

template<typename T> inline void Matrix<T>::init() {
	_rows.resize(_nb_rows);
	T* current_row = &_actual_array[0]; // safe way to get pointer on underlying array
	// compute row pointers
	for (int i = 0; i < _nb_rows; ++i) {
		_rows[i] = current_row;
		current_row += _nb_cols;
	}
}


template<typename T> inline int Matrix<T>::nb_cols() const {
	return _nb_cols;
}

template<typename T> inline int Matrix<T>::nb_rows() const {
	return _nb_rows;
}

template<typename T> inline const T *Matrix<T>::c_contiguous_array() const {
	return _rows[0];
}

template<typename T> inline T *Matrix<T>::c_contiguous_array() {
	return _rows[0];
}

template<typename T> inline const T *Matrix<T>::begin() const {
	return _rows[0];
}

template<typename T> inline T *Matrix<T>::begin() {
	return _rows[0];
}

template<typename T> inline const T *Matrix<T>::end() const {
	return row_by_row_end();
}

template<typename T> inline T *Matrix<T>::end() {
	return row_by_row_end();
}

template<typename T> inline T** Matrix<T>::c_arrays() {
	return &_rows[0];
}

template<typename T> inline const T * const * Matrix<T>::c_arrays() const {
	return &_rows[0];
}

template<typename T> inline T** Matrix<T>::operator ()() {
	return &_rows[0];
}

template<typename T> inline const T * const * Matrix<T>::operator ()() const {
	return &_rows[0];
}

template<typename T> inline const T & Matrix<T>::operator ()(int i, int j) const {
	return _rows[i][j];
}

template<typename T> inline T & Matrix<T>::operator ()(int i, int j) {
	return _rows[i][j];
}

template<typename T> inline T *Matrix<T>::operator [](int row) {
	return _rows[row];
}

template<typename T> inline const T * const Matrix<T>::operator [](int row) const {
	return _rows[row];
}

template<typename T>
inline
bool Matrix<T>::operator == (const Matrix<T>& m) const {
	return _nb_rows == m._nb_rows && _nb_cols == m._nb_cols
			&& _actual_array == m._actual_array;
			//&& std::equal(_actual_array.begin(), _actual_array.end(), m._actual_array.begin() );
}

template<typename T> inline std::ostream & Matrix<T>::printOn(std::ostream & o) const {
	o << "Matrix(" << _nb_rows << ", " << _nb_cols << ")\n";

	for (int i = 0; i < _nb_rows; ++i) {
		for (int j = 0; j < _nb_cols; ++j)
			o << (*this)(i,j) << "\t\t";
		o << "\n";
	}
	return o;
}

template<typename T>
inline Matrix<T>&  Matrix<T>::operator = (Matrix<T> const& source)
{
	if ( &source != this) {
		_nb_rows = source._nb_rows;
		_nb_cols = source._nb_cols;
		_actual_array = source._actual_array;
		init();
	}
	return *this;
}

template<typename T>
template <class Function>
void Matrix<T>::fill(Function f) {
	apply( _3 = f(_1, _2));
}

template<typename T>
void Matrix<T>::fill(T value)
{
	T* cursor = row_by_row_begin();
	T* end = row_by_row_end();
	while ( cursor < end)
		*cursor++ = value;

}


template<typename T>
template <class Function>
void Matrix<T>::apply(Function f) {
	for (int i = 0; i < _nb_rows; ++i)
		for (int j = 0; j < _nb_cols; ++j)
			f(i,j, (*this)(i,j));
}


#endif /* MATRIX_HPP_ */
