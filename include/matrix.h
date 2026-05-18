#include <cstdlib>
#include <vector>
#include <iostream>
#include <utility>
#include <memory>

template <typename T>
class Matrix {
	typedef std::vector<T> Row;
	std::vector<Row> Data;

public:
	unsigned int rows = 0;
	unsigned int cols = 0;

	Matrix(unsigned int rows, unsigned int cols) {
		this->rows = rows;
		this->cols = cols;
		// Asignar memoria a la matriz
		(this->Data).resize(rows);
		for (int i = 0; i < rows; i++) ((this->Data)[i]).resize(cols);
	}

	// Acceso a fila
	Row& operator[](unsigned int row) {
		return Data[row];
	}
	// Acceso a elemento
	T& operator()(unsigned int row, unsigned int col) {
		return Data[row][col];
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
class ShallowPartition {
	public:
		Matrix<T>* original_mat = nullptr;
		std::pair<unsigned int, unsigned int> original_position;
		unsigned int sub_rows = 0;
		unsigned int sub_cols = 0;
	
		ShallowPartition(Matrix<T>& original_mat, unsigned int sub_rows, unsigned int sub_cols, std::pair<unsigned int, unsigned int> original_position) {
			if (original_position.first + sub_rows > original_mat.rows || original_position.second + sub_cols > original_mat.cols) {
				throw std::runtime_error("Partition failed: Block exceeds original matrix dimensions");
			}
			this->original_mat = &original_mat;
			this->sub_rows = sub_rows;
			this->sub_cols = sub_cols;
			this->original_position = original_position;
		}

		// Acceso a elemento de matriz original, pero con original_position siendo el (0,0).
		T& operator()(unsigned int row, unsigned int col) {
			return (*original_mat)(original_position.first + row, original_position.second + col);
		}
		// Debug
		void print_contents() {
			for (int i = 0; i < this->sub_rows; i++) {
				for (int j = 0; j < this->sub_cols; j++) {
					std::cout << (*this)(i,j) << ",\t";
				}
				std::cout << std::endl;
			}
		}
};

// Result = Mat_A + Mat_B
template <typename T>
void matrix_sum(Matrix<T>& destination, Matrix<T>& mat_A, Matrix<T>& mat_B) {
	if (mat_A.rows != mat_B.rows || mat_A.cols != mat_B.cols || mat_A.rows != destination.rows || mat_A.cols != destination.cols)
		throw std::runtime_error("Sum error: Mismatch in matrix dimensions");
	for (int i = 0; i < mat_A.rows; i++) {
		for (int j = 0; j < mat_A.cols; j++) {
			destination(i,j) = mat_A(i,j) + mat_B(i,j);
		}
	}
}

// Result = Mat_A - Mat_B
template <typename T>
void matrix_sub(Matrix<T>& destination, Matrix<T>& mat_A, Matrix<T>& mat_B) {
	if (mat_A.rows != mat_B.rows || mat_A.cols != mat_B.cols || mat_A.rows != destination.rows || mat_A.cols != destination.cols)
		throw std::runtime_error("Subtraction error: Mismatch in matrix dimensions");
	for (int i = 0; i < mat_A.rows; i++) {
		for (int j = 0; j < mat_A.cols; j++) {
			destination(i,j) = mat_A(i,j) - mat_B(i,j);
		}
	}
}

// Copy a nxm block from the source partition into a destination matrix, and also specifying which position (i,j) to start copying on (left to right, top to bottom).
// The (nxm) block must fulfill the following conditions:
// 1. n < dest_part.rows
// 2. m < dest_part.cols
template <typename T>
Matrix<T>& matrix_block_copy(Matrix<T>& dest, ShallowPartition<T>& source_part, std::pair<unsigned int, unsigned int> dest_position) {
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
