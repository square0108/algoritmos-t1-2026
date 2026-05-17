#include "matrix.h"
#include <iostream>

int main() {
	Matrix<int> *A = new Matrix<int>(8,8);
	Matrix<int> *B = new Matrix<int>(8,8);
	for (int i = 0; i < A->rows; i++) {
		for (int j = 0; j < A->cols; j++) {
			(*A)[i][j] = i+j+1;
			(*B)[i][j] = i+j+1;
		}
	}
	std::cout << (*A)[0][0] << std::endl;
	std::cout << "Matrix rows: " << A->rows << "; Matrix cols: " << A->cols << std::endl;
	A->print_contents();
	std::cout << "Partition test" << std::endl;
	ShallowPartition<int> part(*A,4,4,std::make_pair(4,4));
	part.print_contents();
	std::cout << "copy_block_matrix test" << std::endl;
	matrix_block_copy(*A,part,std::make_pair(7,7));
	A->print_contents();

	return 0;
}
