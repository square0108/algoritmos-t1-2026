#include "matrix.h"
#include "strassen.h"
#include <string>
#include <cstdlib>
#include <iostream>

int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cerr << "Usage: ./executable [strassen | rowcol | hybrid] [n_matrix_size]" << std::endl;
		exit(-1);
	}
	const std::string mult_type = argv[1];
	const size_t MATRIX_SIZE = atoi(argv[2]);
	/*
	std::cout << "Partition test" << std::endl;
	ShallowPartition<int> part(*A,4,4,std::make_pair(4,4));
	part.print_contents();
	std::cout << "copy_block_matrix test" << std::endl;
	matrix_block_copy(*A,part,std::make_pair(3,3));
	A->print_contents();
	*/
	std::cout << "matrix mult test (no print)..." << std::endl;
	Matrix<int> *D = new Matrix<int>(MATRIX_SIZE,MATRIX_SIZE);
	Matrix<int> *E = new Matrix<int>(MATRIX_SIZE,MATRIX_SIZE);
	Matrix<int> *F = new Matrix<int>(MATRIX_SIZE,MATRIX_SIZE);
	for (int i = 0; i < D->rows; i++) {
		for (int j = 0; j < D->cols; j++) {
			(*D)(i,j) = i+j+1;
			(*E)(i,j) = i+j+1;
		}
	}
	if (mult_type == "strassen") strassen_mult<int>(*F,*D,*E);
	else if (mult_type == "rowcol") matrix_multRowCol<int>(*F,*D,*E);
	std::cout << "...success!" << std::endl;
	delete D; delete E; delete F;
	return 0;
}
