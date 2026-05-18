#include "matrix.h"
#include <memory>
#include <utility>

#ifndef STRASSEN_H
#define STRASSEN_H

enum QuadPartition { q11, q12, q21, q22 };
const std::vector<QuadPartition> quadrants = {q11,q12,q21,q22};

std::pair<unsigned int, unsigned int> get_quad_coords(QuadPartition q, unsigned int square_mat_size);

template <typename T>
// Computes AB = output
void strassen_mult(Matrix<T>& output, Matrix<T>& A, Matrix<T>& B) {
	if (A.rows != A.cols || (A.rows & ((A.rows)-1)) != 0
		|| B.rows != B.cols || (B.rows & ((B.rows)-1)) != 0) {
		std::cerr << "Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2." << std::endl;
		throw	std::runtime_error("Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2.");
	}
	unsigned int n = A.rows;

	if (n == 1) { 
		output(0,0) = A(0,0) * B(0,0);
	}
	else {
		// Matrices de partición se crean como objetos nuevos
		std::unique_ptr<Matrix<T>> A_parts[4];
		std::unique_ptr<Matrix<T>> B_parts[4];
		std::unique_ptr<Matrix<T>> output_parts[4];
		std::pair<unsigned int, unsigned int> copy_destination = std::make_pair(0,0);
		for (QuadPartition q : quadrants) {
			A_parts[q] = std::make_unique<Matrix<T>>(n/2, n/2);
			B_parts[q] = std::make_unique<Matrix<T>>(n/2, n/2);
			output_parts[q] = std::make_unique<Matrix<T>>(n/2, n/2);
			
			// ShallowPartition selecciona el bloque, y matrix_block_copy lo copia a una matriz nueva n/2 x n/2
			ShallowPartition<T> partA(A,n/2,n/2,get_quad_coords(q,n));  
			matrix_block_copy(*(A_parts[q]),partA,copy_destination);
			ShallowPartition<T> partB(B,n/2,n/2,get_quad_coords(q,n));  
			matrix_block_copy(*(B_parts[q]),partB,copy_destination);
		}
		// Reserva y cálculo de matrices suma
		std::unique_ptr<Matrix<T>> Smat[10];
		for (int i = 0; i < 10; i++) {
			Smat[i] = std::make_unique<Matrix<T>>(n/2,n/2);
		}
		// matrix_sum(destino, operando1, operando2)
		matrix_sum(*(Smat[0]),*A_parts[QuadPartition::q11],*A_parts[QuadPartition::q22]);
		matrix_sum(*(Smat[1]),*B_parts[QuadPartition::q11],*B_parts[QuadPartition::q22]);
		matrix_sum(*(Smat[2]),*A_parts[QuadPartition::q21],*A_parts[QuadPartition::q22]);
		matrix_sub(*(Smat[3]),*B_parts[QuadPartition::q12],*B_parts[QuadPartition::q22]);
		matrix_sub(*(Smat[4]),*B_parts[QuadPartition::q21],*B_parts[QuadPartition::q11]);
		matrix_sum(*(Smat[5]),*A_parts[QuadPartition::q11],*A_parts[QuadPartition::q12]);
		matrix_sub(*(Smat[6]),*A_parts[QuadPartition::q21],*A_parts[QuadPartition::q11]);
		matrix_sum(*(Smat[7]),*B_parts[QuadPartition::q11],*B_parts[QuadPartition::q12]);
		matrix_sub(*(Smat[8]),*A_parts[QuadPartition::q12],*A_parts[QuadPartition::q22]);
		matrix_sum(*(Smat[9]),*B_parts[QuadPartition::q21],*B_parts[QuadPartition::q22]);

		// M_1 ... M_7
		std::unique_ptr<Matrix<T>> Mmat[7];
		for (int i = 0; i < 7; i++) Mmat[i] = std::make_unique<Matrix<T>>(n/2,n/2);
		strassen_mult(*Mmat[0],*(Smat[0]),*(Smat[1]));
		strassen_mult(*Mmat[1],*(Smat[2]),*(B_parts[QuadPartition::q11]));
		strassen_mult(*Mmat[2],*(A_parts[QuadPartition::q11]),*(Smat[3]));
		strassen_mult(*Mmat[3],*(A_parts[QuadPartition::q22]),*(Smat[4]));
		strassen_mult(*Mmat[4],*(Smat[5]),*(B_parts[QuadPartition::q22]));
		strassen_mult(*Mmat[5],*(Smat[6]),*(Smat[7]));
		strassen_mult(*Mmat[6],*(Smat[8]),*(Smat[9]));

		// Cálculo de los elementos de C
		// C_1,1
		matrix_sum(*output_parts[q11], *Mmat[0], *Mmat[3]);
		matrix_sub(*output_parts[q11], *output_parts[q11], *Mmat[4]);
		matrix_sum(*output_parts[q11], *output_parts[q11], *Mmat[6]);
		// C_1,2
		matrix_sum(*output_parts[q12], *Mmat[2], *Mmat[4]);
		// C_2,1
		matrix_sum(*output_parts[q21], *Mmat[1], *Mmat[3]);
		// C_2,2
		matrix_sub(*output_parts[q22], *Mmat[0], *Mmat[1]);
		matrix_sum(*output_parts[q22], *output_parts[q22], *Mmat[2]);
		matrix_sum(*output_parts[q22], *output_parts[q22], *Mmat[5]);
		// Combinación en C
		for (unsigned int i = 0; i < n/2; i++) {
			for (unsigned int j = 0; j < n/2; j++) {
				output(i,j) = (*output_parts[q11])(i,j);
				output(i,j+n/2) = (*output_parts[q12])(i,j);
				output(i+n/2,j) = (*output_parts[q21])(i,j);
				output(i+n/2,j+n/2) = (*output_parts[q22])(i,j);
			}
		}

		// Cleanup de memoria se delega a los smart pointerssSSs
	}
	// aqui se "retorna" output
}

std::pair<unsigned int, unsigned int> get_quad_coords(QuadPartition q, unsigned int square_mat_size) {
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
			throw std::runtime_error("Partition coordinate assignment error");
	};
}

#endif