#include "matrix.h"
#include <memory>
#include <utility>

#ifndef STRASSEN_H
#define STRASSEN_H

enum Quadrant { q11, q12, q21, q22 };
const std::vector<Quadrant> quadrants = {q11,q12,q21,q22};

std::pair<size_t, size_t> get_quad_coords(Quadrant q, size_t square_mat_size);

/* Utilizamos templates para los types de las matrices A, B y C para que strassen_mult acepte tanto PartitionViews como Matrices
En la práctica, MatA, MatB y MatC deben ser alguna permutación de los tipos Matrix y PartitionView
*/
template <typename DataType, typename MatA, typename MatB, typename MatC>
void strassen_mult(MatC& output, MatA& A, MatB& B, size_t n0_threshold) {
	// Error checking de dimensiones matriciales
	if (A.rows != A.cols || (A.rows & ((A.rows)-1)) != 0 /* <==> if not potencia de 2*/
		|| B.rows != B.cols || (B.rows & ((B.rows)-1)) != 0) {
		std::cerr << "Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2." << std::endl;
		exit(-1);
	}
	size_t n = A.rows;
	if (n <= n0_threshold) { 
		matrix_multRowCol<DataType>(output,A,B);
	}
	else {
		// Matrices de partición no son objetos nuevos, si no que acceden a data de matriz original
		// Cada partición construida recibe: el objeto particionado, las dimensiones, y las coordenadas de offset (la posición 0,0 de la partición respecto a la matriz original)
		PartitionView<DataType> A11(A,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> A12(A,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> A21(A,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> A22(A,n/2,n/2,get_quad_coords(q22,n));
		PartitionView<DataType> B11(B,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> B12(B,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> B21(B,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> B22(B,n/2,n/2,get_quad_coords(q22,n));
		PartitionView<DataType> C11(output,n/2,n/2,get_quad_coords(q11,n)); PartitionView<DataType> C12(output,n/2,n/2,get_quad_coords(q12,n));
		PartitionView<DataType> C21(output,n/2,n/2,get_quad_coords(q21,n)); PartitionView<DataType> C22(output,n/2,n/2,get_quad_coords(q22,n));

		// Matrices temporales para computar sumas
		Matrix<DataType>* temp1 = new Matrix<DataType>(n/2,n/2);
		Matrix<DataType>* temp2 = new Matrix<DataType>(n/2,n/2);
		
		// M_1 = (A11+A22)*(B11+B22)
		matrix_sum(*temp1,A11,A22);
		matrix_sum(*temp2,B11,B22);
		Matrix<DataType>* M1 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M1,*temp1,*temp2,n0_threshold);
		// M_2
		matrix_sum(*temp1,A21,A22);
		Matrix<DataType>* M2 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M2,*temp1,B11,n0_threshold);	
		// M_3
		matrix_sub(*temp1,B12,B22);
		Matrix<DataType>* M3 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M3,A11,*temp1,n0_threshold);
		// M_4
		matrix_sub(*temp1,B21,B11);
		Matrix<DataType>* M4 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M4,A22,*temp1,n0_threshold);
		// M_5
		matrix_sum(*temp1,A11,A12);
		Matrix<DataType>* M5 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M5,*temp1,B22,n0_threshold);
		// M_6
		matrix_sub(*temp1,A21,A11);
		matrix_sum(*temp2,B11,B12);
		Matrix<DataType>* M6 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M6,*temp1,*temp2,n0_threshold);
		// M_7
		matrix_sub(*temp1,A12,A22);
		matrix_sum(*temp2,B21,B22);
		Matrix<DataType>* M7 = new Matrix<DataType>(n/2,n/2);
		strassen_mult<DataType>(*M7,*temp1,*temp2,n0_threshold);

		delete temp1; delete temp2;
		// Cálculo de los elementos de C. Como las particiones son vistas a la data de la matriz original, escribir en C11 escribe directamente a C.
		// C_1,1
		matrix_sum(C11, *M1, *M4);
		matrix_sub(C11, C11, *M5);
		matrix_sum(C11, C11, *M7);
		// C_1,2
		matrix_sum(C12, *M3, *M5);
		// C_2,1
		matrix_sum(C21, *M2, *M4);
		// C_2,2
		matrix_sub(C22, *M1, *M2);
		matrix_sum(C22, C22, *M3);
		matrix_sum(C22, C22, *M6);

		delete M1; delete M2; delete M3; delete M4; delete M5; delete M6; delete M7;
	}
	// aqui se "retorna" output
}

/* Retorna las coordenadas de origen (offset) de la matriz partición cuadrada de su cuadrante respecto a la matriz padre */
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