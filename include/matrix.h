#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <vector>
#include <iostream>
#include <utility>
#include <memory>
#include <type_traits>

template <typename DataType>
struct Matrix {
	std::vector<DataType> Data;
	size_t rows;
	size_t cols;

	Matrix(size_t rows, size_t cols) {
		this->rows = rows;
		this->cols = cols;
		// Asignar memoria a la matriz
		(this->Data).resize(rows*cols);
	}
	// Acceso a elemento
	DataType& operator()(size_t row, size_t col) {
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

template <typename DataType>
/* Se comporta igual que una Matrix<DataType> en aritmética de matrices, pero para acceder a elementos, utiliza
el miemrbo Data de la matriz que está particionando ("matriz padre")*/
struct PartitionView {
		Matrix<DataType>* parent_mat = nullptr;
		std::pair<size_t, size_t> offset; // Posicion (0,0) de la submatriz dentro de la matriz padre.
		size_t rows;
		size_t cols;

		template <typename MatrixType>
		PartitionView(MatrixType& mat, size_t sub_rows, size_t sub_cols, std::pair<size_t, size_t> sub_offset) {
			// Partición de otra partición
			if constexpr (std::is_same_v<MatrixType, PartitionView<DataType>>) {
				this->parent_mat = mat.parent_mat;
				this->offset.first = mat.offset.first + sub_offset.first;
				this->offset.second = mat.offset.second + sub_offset.second;
			}
			// Partición de una matriz con data propia
			else if constexpr (std::is_same_v<MatrixType, Matrix<DataType>>) {
				this->parent_mat = &mat;
				this->offset.first = sub_offset.first;
				this->offset.second = sub_offset.second;
			}
			this->rows = sub_rows;
			this->cols = sub_cols;

			// Error checking
			if (sub_rows < 1 || sub_cols < 1) {
				std::cerr << "Can't create a size 0 partition" << std::endl;
				exit(-1);
			}
			if (sub_offset.first + sub_rows > parent_mat->rows || sub_offset.second + sub_cols > parent_mat->cols) {
				std::cerr << "Partition exceeds matrix dimensions. Input: " 
				<< sub_rows << ", " << sub_cols 
				<< "; Matrix size: " 
				<< parent_mat->rows << ", " << parent_mat->cols 
				<< "; Offset pos: (" 
				<< sub_offset.first << "," << sub_offset.second << ")" << std::endl;
				throw std::runtime_error("Partition failed: Block exceeds original matrix dimensions");
			}
		}

		// Acceso a elemento de matriz original, pero con partition_pos siendo el (0,0) de la partición.
		DataType& operator()(size_t row, size_t col) {
			return (*parent_mat)(offset.first + row, offset.second + col);
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
template <typename DestType, typename MatA, typename MatB>
void matrix_sum(DestType& destination, MatA& mat_A, MatB& mat_B) {
	if (mat_A.rows != mat_B.rows || mat_A.cols != mat_B.cols || mat_A.rows != destination.rows || mat_A.cols != destination.cols) {
		std::cerr << ("Sum error: Mismatch in matrix dimensions");
		exit(-1);
	}
	for (size_t i = 0; i < mat_A.rows; i++) {
		for (size_t j = 0; j < mat_A.cols; j++) {
			destination(i,j) = mat_A(i,j) + mat_B(i,j);
		}
	}
}

// Result = Mat_A - Mat_B
template <typename DestType, typename MatA, typename MatB>
void matrix_sub(DestType& destination, MatA& mat_A, MatB& mat_B) {
	if (mat_A.rows != mat_B.rows || mat_A.cols != mat_B.cols || mat_A.rows != destination.rows || mat_A.cols != destination.cols) {
		std::cerr << ("Subtraction error: Mismatch in matrix dimensions");
		exit(-1);
	}
	for (size_t i = 0; i < mat_A.rows; i++) {
		for (size_t j = 0; j < mat_A.cols; j++) {
			destination(i,j) = mat_A(i,j) - mat_B(i,j);
		}
	}
}

// Moving all this shit to templates cus fuck OOP
template <typename DataType, typename DestType, typename MatA, typename MatB>
void matrix_multRowCol(DestType& destination, MatA& A, MatB& B) {
	if (A.cols != B.rows) {
		std::cerr << ("RowCol multiplication: Mismatch in matrix dimensions. Ensure A's cols == B's rows");
		exit(-1);
	}
	else if (destination.rows != A.rows || destination.cols != B.cols) {
		std::cerr << ("RowCol multiplication: Mismatch in expected dimensions of output variable.");
		exit(-1);
	}
	 else {
		for (size_t i = 0; i < destination.rows; i++) {
			for (size_t j = 0; j < destination.cols; j++) {
				DataType result = (DataType) 0;
				for (size_t k = 0; k < A.cols; k++) {
					result += A(i,k) * B(k,j);
				}
				destination(i,j) = result;
			}
		}
	}
}

// Copy a nxm block from the source partition into a destination matrix, and also specifying which position (i,j) to start copying on (left to right, top to bottom).
// The (nxm) block must fulfill the following conditions:
// 1. n < dest_part.rows
// 2. m < dest_part.cols
/*
template <typename T>
Matrix<DataType>& matrix_block_copy(Matrix<DataType>& dest, ShallowPartition<DataType>& source_part, std::pair<unsigned int, unsigned int> dest_position) {
	auto dest_row = dest_position.first;
	auto dest_col = dest_position.second;
	if (source_part.sub_rows + dest_row > dest.rows || source_part.sub_cols + dest_col > dest.cols) 
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