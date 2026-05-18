#include "matrix.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <utility>

int main(int argc, char* argv[]) {
	Matrix<int> A(4,4);
  Matrix<int> B(8,8);
  int counter = 1;
  for (size_t i = 0; i < A.get_rows(); i++) {
    for (size_t j = 0; j < A.get_cols(); j++) {
      A(i,j) = counter;
      B(i,j) = counter;
      counter++;
    }
  }
  counter = 1;
  for (size_t i = 0; i < B.get_rows(); i++) {
    for (size_t j = 0; j < B.get_cols(); j++) {
      B(i,j) = counter;
      counter++;
    }
  }
  VirtualPartition<int> C(B,4,4,std::make_pair(0,0));
  B.print_contents();
  A.print_contents();
  C.print_contents();
  matrix_sum(A,A,C);
  A.print_contents();
  matrix_sub(C,A,A);
  C.print_contents();
  B.print_contents();
	return 0;
}
