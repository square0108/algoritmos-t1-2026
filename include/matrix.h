#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <vector>
#include <iostream>
#include <utility>
#include <memory>

/* Visión:
*	- Matrix es una matriz "común" de toda la vida
* - VirtualPartition es una submatriz que en vez de crear otro container,
* 	utiliza el puntero de su matriz padre para acceder datos
*	- Ambas heredan de MatrixOperand para que las operaciones de suma, resta y multiplicación
* 	acepten objetos MatrixOperand y traten a matrices y particiones de igual forma.
*	Esto significa que escrituras a VirtualPartition en realidad escriben directo a su matriz padre, y
* resulta en menos pasos de copia. Asumí que algo así era a lo que hacía referencia el CLRS con usar
* "index ranges" para reducir un poco el trabajo Theta(n^2)
*/

template <typename T>
class MatrixOperand {
protected:
	size_t rows;
	size_t cols;
public:
	virtual ~MatrixOperand() = default;
	virtual T& operator()(size_t row, size_t col) = 0;
	size_t get_rows() {
		return rows;
	}
	size_t get_cols() {
		return cols;
	}
};

template <typename T>
class Matrix : public MatrixOperand<T> {
	std::vector<T> Data;
public:
	Matrix(size_t rows, size_t cols) {
		this->rows = rows;
		this->cols = cols;
		// Asignar memoria a la matriz
		(this->Data).resize(rows*cols);
	}

	// Acceso a elemento
	T& operator()(size_t row, size_t col) override {
		return Data[row * this->cols + col];
	}

	// Debug
  void print_contents() {
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < this->cols; j++) {
				std::cout << (*this)(i,j) << ",\t";
			}
			std::cout << std::endl;
		}
	}
};

template <typename T>
// "Virtual view" of a matrix block/partition inside of a larger matrix, `original_mat`
class VirtualPartition : public MatrixOperand<T> {
	private:
		Matrix<T>* parent_mat = nullptr; // En vez de container Data, usar Data de matriz particionada
		std::pair<size_t, size_t> partition_pos; // Posicion (0,0) de la submatriz dentro de la matriz padre.
	public:
		VirtualPartition(Matrix<T>& original_mat, size_t sub_rows, size_t sub_cols, std::pair<size_t, size_t> original_position) {
			if (original_position.first + sub_rows > original_mat.get_rows() || original_position.second + sub_cols > original_mat.get_cols()) {
				throw std::runtime_error("Partition failed: Block exceeds original matrix dimensions");
			}
			this->parent_mat = &original_mat;
			this->rows = sub_rows;
			this->cols = sub_cols;
			this->partition_pos = original_position;
		}
		// Acceso a elemento de matriz original, pero con partition_pos siendo el (0,0) de la partición.
		T& operator()(size_t row, size_t col) override {
			return (*parent_mat)(partition_pos.first + row, partition_pos.second + col);
		}
		size_t get_parent_rows() {
			return parent_mat->get_rows();
		}
		size_t get_parent_cols() {
			return parent_mat->get_cols();
		}
		Matrix<T>* get_parent_matrix() {
			return parent_mat;
		}
		std::pair<size_t, size_t> get_pos_within_parent() {
			return partition_pos;
		}
		// Debug
		void print_contents() {
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					std::cout << (*this)(i,j) << ",\t";
				}
				std::cout << std::endl;
			}
		}
};

// Result = Mat_A + Mat_B
template <typename T>
void matrix_sum(MatrixOperand<T>& destination, MatrixOperand<T>& mat_A, MatrixOperand<T>& mat_B) {
	if (mat_A.get_rows() != mat_B.get_rows() || mat_A.get_cols() != mat_B.get_cols() || mat_A.get_rows() != destination.get_rows() || mat_A.get_cols() != destination.get_cols())
		throw std::runtime_error("Sum error: Mismatch in matrix dimensions");
	for (int i = 0; i < mat_A.get_rows(); i++) {
		for (int j = 0; j < mat_A.get_cols(); j++) {
			destination(i,j) = mat_A(i,j) + mat_B(i,j);
		}
	}
}

// Result = Mat_A - Mat_B
template <typename T>
void matrix_sub(MatrixOperand<T>& destination, MatrixOperand<T>& mat_A, MatrixOperand<T>& mat_B) {
	if (mat_A.get_rows() != mat_B.get_rows() || mat_A.get_cols() != mat_B.get_cols() || mat_A.get_rows() != destination.get_rows() || mat_A.get_cols() != destination.get_cols())
		throw std::runtime_error("Subtraction error: Mismatch in matrix dimensions");
	for (unsigned int i = 0; i < mat_A.get_rows(); i++) {
		for (unsigned int j = 0; j < mat_A.get_cols(); j++) {
			destination(i,j) = mat_A(i,j) - mat_B(i,j);
		}
	}
}

template <typename T>
void matrix_multRowCol(MatrixOperand<T>& destination, MatrixOperand<T>& A, MatrixOperand<T>& B) {
	if (A.get_cols() != B.get_rows())
		throw std::runtime_error("RowCol multiplication: Mismatch in matrix dimensions. Ensure A's cols == B's rows");
	else if (destination.get_rows() != A.get_rows() || destination.get_cols() != B.get_cols())
		throw std::runtime_error("RowCol multiplication: Mismatch in expected dimensions of output variable.");
	else {
		for (size_t i = 0; i < destination.get_rows(); i++) {
			for (size_t j = 0; j < destination.get_cols(); j++) {
				T result = (T) 0;
				for (size_t k = 0; k < A.get_cols(); k++) {
					result += A(i,k) * B(k,j);
				}
				destination(i,j) = result;
			}
		}
	}
}

// Copy a nxm block from the source partition into a destination matrix, and also specifying which position (i,j) to start copying on (left to right, top to bottom).
// The (nxm) block must fulfill the following conditions:
// 1. n < dest_part.get_rows()
// 2. m < dest_part.cols
/*
template <typename T>
Matrix<T>& matrix_block_copy(Matrix<T>& dest, ShallowPartition<T>& source_part, std::pair<unsigned int, unsigned int> dest_position) {
	auto dest_row = dest_position.first;
	auto dest_col = dest_position.second;
	if (source_part.sub_rows + dest_row > dest.get_rows() || source_part.sub_cols + dest_col > dest.cols) 
		throw std::runtime_error("Copy error: Destination cannot fit the source block");

	// begin copy
	for (unsigned int i = 0; i < source_part.sub_rows; i++) {
		for (unsigned int j = 0; j < source_part.sub_cols; j++) {
			dest(i+dest_row,j+dest_col) = source_part(i,j);
		}
	}
	return dest;
}
*/

#endif