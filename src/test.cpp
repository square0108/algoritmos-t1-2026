#include <string>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <queue>
#include <random>

#include "matrix.h"
#include "strassen.h"

const std::vector<std::string> VALID_OPERATIONS = {"strassen","rowcol"};
const size_t DEFAULT_N0 = 1;
const long long RANDOM_INT64_BOUNDS[2] = {1,5};

void error_msg_exit();

// Para generar valores aleatorios en matriz prueba
std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<std::int64_t> u_distr(RANDOM_INT64_BOUNDS[0],RANDOM_INT64_BOUNDS[1]);

int main(int argc, char* argv[]) {
	// Flags de opciones de input
	bool print_flag = false;

	// Variables de multiplicación recibidas por argumentos
	size_t n0 = DEFAULT_N0;
	size_t matrix_size;
	std::string mult_type;

	// Captación de inputs
	if (argc < 3 || argc > 5) error_msg_exit();
	std::queue<std::string> args;
	for (int i_arg = 1; i_arg < argc; i_arg++) {
		args.push(std::string(argv[i_arg]));
	}

	// Validar tipo de multiplicacion
	mult_type = args.front(); args.pop();
	if (std::find(VALID_OPERATIONS.begin(), VALID_OPERATIONS.end(), mult_type) == VALID_OPERATIONS.end()) error_msg_exit();

	// Validar tamaño de matriz es 2^k
	std::string matrix_size_arg = args.front(); args.pop();
	std::stringstream sstream(matrix_size_arg);
	sstream >> matrix_size;
	if ((matrix_size & ((matrix_size)-1)) != 0) error_msg_exit();

	// Validar argumentos opcionales
	while (!args.empty()) {
		std::string op = args.front(); args.pop();
		if (op != "--print-result") {
			std::stringstream sstream2(matrix_size_arg);
			size_t n0;
			sstream >> n0;
		}
		if (op == "--print-result") 
			print_flag = true;
	}

	/* Matrix creation */
	Matrix<long long> *A = new Matrix<long long>(matrix_size,matrix_size);
	Matrix<long long> *B = new Matrix<long long>(matrix_size,matrix_size);
	Matrix<long long> *C = new Matrix<long long>(matrix_size,matrix_size);
	for (int i = 0; i < A->rows; i++) {
		for (int j = 0; j < A->cols; j++) {
			// placeholder values, add varied tests later
			(*A)(i,j) = u_distr(rng);
			(*B)(i,j) = u_distr(rng);
		}
	}
	if (mult_type == "strassen") {
		strassen_mult<long long>(*C,*A,*B,n0);
	}
	else if (mult_type == "rowcol") matrix_multRowCol<long long>(*C,*A,*B);
	if (print_flag) {
		std::cout << "Matrix A: " << std::endl;
		A->print_contents();
		std::cout << "Matrix B: " << std::endl;
		B->print_contents();
		std::cout << "Matrix AB: " << std::endl;
		C->print_contents();
	}

	delete A; delete B; delete C;
	return 0;
}

void error_msg_exit() {
	std::cerr << "Usage: ./executable <strassen | rowcol> <MATRIX_SIZE> [STRASSEN_THRESHOLD] [--print-result]" << std::endl;
	std::cerr << "Requirements:\n"
						<< "- Ensure MATRIX_SIZE is a power of 2\n"
						<< "- Ensure STRASSEN_THRESHOLD is a power of 2" 
						<< std::endl;
	exit(-1);
}