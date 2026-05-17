#include "matrix.h"
#include <memory>
#include <utility>

enum Quadrants { q11 = 0, q12 = 1, q21 = 2, q22 = 3 };

template <typename T>
// Computes AB = C
std::unique_ptr<Matrix<T>> strassen_mult(Matrix<T>& A, Matrix<T>& B) {
	if (A.rows != A.cols || (A.rows & ((A.rows)-1)) != 0) {
		throw	std::runtime_error("Invalid matrix multiplication parameters: Ensure A and B have dimensions n x n and that n is a power of 2.");
	}
	auto C = std::make_unique<Matrix<T>(A.rows, A.rows)>;
	unsigned int n = A.rows;

	if (n == 1) { 
		(*C)[0][0] = A[0][0] * B[0][0];
	}
	else {
		// Definir coordenadas de origen (0,0) de cada matriz partición respecto a su supermatriz
		std::pair<unsigned int, unsigned int> partition_coords[] = {
			std::make_pair(0,0),		// top left
			std::make_pair(0,n/2),	// top right
			std::make_pair(n/2,0),	// bot left
			std::make_pair(n/2,n/2)	// bot right
		};
		// Matrices de partición se crean como objetos nuevos
		std::unique_ptr<Matrix<T>> A_parts[4];
		std::unique_ptr<Matrix<T>> B_parts[4];
		std::pair<unsigned int, unsigned int> copy_dest_coord = std::make_pair(0,0);
		for (int i = 0; i < 4; i++) {
			A_parts[i] = std::make_unique<Matrix<T>(n/2, n/2)>;
			B_parts[i] = std::make_unique<Matrix<T>(n/2, n/2)>;
			
			// ShallowPartition selecciona el bloque en tiempo O(1), y matrix_block_copy lo copia a una matriz n/2.
			ShallowPartition<T> partA(A,n/2,n/2,partition_coords[i]);  
			matrix_block_copy(*(A_parts[i]),partA,copy_dest_coord);
			ShallowPartition<T> partB(B,n/2,n/2,partition_coords[i]);  
			matrix_block_copy(*(B_parts[i]),partB,copy_dest_coord);
			// C parts no se copia pues se sobreescribirá con las sumas de M_1 ... M_7
		}
		// Reserva y cálculo de matrices suma
		std::unique_ptr<Matrix<T>> S[10];
		for (int i = 0; i < 10; i++) {
			S[i] = std::make_unique<Matrix<T>(n/2,n/2)>;
		}
		// matrix_sum(destino, operando1, operando2)
		matrix_sum(*(S[0]),A_parts[Quadrants::q11],A_parts[Quadrants::q22]);
		matrix_sum(*(S[1]),B_parts[Quadrants::q11],B_parts[Quadrants::q22]);
		matrix_sum(*(S[2]),A_parts[Quadrants::q21],A_parts[Quadrants::q22]);
		matrix_sub(*(S[3]),B_parts[Quadrants::q12],B_parts[Quadrants::q22]);
		matrix_sub(*(S[4]),B_parts[Quadrants::q21],B_parts[Quadrants::q11]);
		matrix_sum(*(S[5]),A_parts[Quadrants::q11],A_parts[Quadrants::q12]);
		matrix_sub(*(S[6]),A_parts[Quadrants::q21],A_parts[Quadrants::q11]);

		// las 7 esferas del dragon

	}
}


