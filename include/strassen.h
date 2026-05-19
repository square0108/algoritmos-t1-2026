#include "matrix.h"
#include <memory>
#include <utility>

#ifndef STRASSEN_H
#define STRASSEN_H

enum Quadrant { q11, q12, q21, q22 };
const std::vector<Quadrant> quadrants = {q11,q12,q21,q22};

std::pair<size_t, size_t> get_quad_coords(Quadrant q, size_t square_mat_size);

template <typename DataType, typename MatA, typename MatB, typename DestType>
// Computes AB = output
void strassen_mult(DestType& output, MatA& A, MatB& B) {
	// Error checking de dimensions matriciales
	if (A.rows != A.cols || (A.rows & ((A.rows)-1)) != 0 /* <==> if not potencia de 2*/
		|| B.rows != B.cols || (B.rows & ((B.rows)-1)) != 0) {
		std::cerr << "Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2." << std::endl;
		throw	std::runtime_error("Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2.");
	}
	size_t n = A.rows;
	if (n == 1) { 
		output(0,0) = A(0,0) * B(0,0);
	}
	else {
		// Matrices de partición no son objetos nuevos, si no que acceden a data de matriz original
		PartitionView<DataType> A11(A,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> A12(A,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> A21(A,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> A22(A,n/2,n/2,get_quad_coords(q22,n));
		PartitionView<DataType> B11(B,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> B12(B,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> B21(B,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> B22(B,n/2,n/2,get_quad_coords(q22,n));
		PartitionView<DataType> C11(output,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> C12(output,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> C21(output,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> C22(output,n/2,n/2,get_quad_coords(q22,n));

		// Arreglo para matrices M_1 ... M_7
		Matrix<DataType>* Mmat[7];
		for (int i = 0; i < 7; i++) Mmat[i] = new Matrix<DataType>(n/2,n/2);
		// Matrices temporales para computar sumas
		Matrix<DataType>* temp1 = new Matrix<DataType>(n/2,n/2);
		Matrix<DataType>* temp2 = new Matrix<DataType>(n/2,n/2);
		// M_1 = (A11+A22)*(B11+B22)
		matrix_sum(*temp1,A11,A22);
		matrix_sum(*temp2,B11,B22);
		strassen_mult<DataType>(*Mmat[0],*temp1,*temp2);

		// M_2
		matrix_sum(*temp1,A21,A22);
		strassen_mult<DataType>(*Mmat[1],*temp1,B11);
		
		// M_3
		matrix_sub(*temp1,B12,B22);
		strassen_mult<DataType>(*Mmat[2],A11,*temp1);
		
		// M_4
		matrix_sub(*temp1,B21,B11);
		strassen_mult<DataType>(*Mmat[3],A22,*temp1);

		// M_5
		matrix_sum(*temp1,A11,A12);
		strassen_mult<DataType>(*Mmat[4],*temp1,B22);

		// M_6
		matrix_sub(*temp1,A21,A11);
		matrix_sum(*temp2,B11,B12);
		strassen_mult<DataType>(*Mmat[5],*temp1,*temp2);

		// M_7
		matrix_sub(*temp1,A12,A22);
		matrix_sum(*temp2,B21,B22);
		strassen_mult<DataType>(*Mmat[6],*temp1,*temp2);

		// Cálculo de los elementos de C. Como output_parts son particiones virtuales, no es necesario escribir a output.
		// C_1,1
		matrix_sum(C11, *Mmat[0], *Mmat[3]);
		matrix_sub(C11, C11, *Mmat[4]);
		matrix_sum(C11, C11, *Mmat[6]);
		// C_1,2
		matrix_sum(C12, *Mmat[2], *Mmat[4]);
		// C_2,1
		matrix_sum(C21, *Mmat[1], *Mmat[3]);
		// C_2,2
		matrix_sub(C22, *Mmat[0], *Mmat[1]);
		matrix_sum(C22, C22, *Mmat[2]);
		matrix_sum(C22, C22, *Mmat[5]);

		// Cleanup de memoria
		for (int i = 0; i < 7; i++) delete Mmat[i];
		delete temp1;
		delete temp2;
	}
	// aqui se "retorna" output
}

std::pair<size_t, size_t> get_quad_coords(Quadrant q, size_t square_mat_size) {
	switch (q) {
		case q11:
			return std::make_pair(0,0);
		case q12:
			return std::make_pair(0,square_mat_size/2);
		case q21:
			return std::make_pair(square_mat_size/2,0);
		case q22:
			return std::make_pair(square_mat_size/2, square_mat_size/2);
		default:
			std::cerr << "Partition coordinate assignment error";
			exit(-1);
	};
}

#endif