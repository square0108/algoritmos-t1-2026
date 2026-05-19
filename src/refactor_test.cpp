#include "matrix.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <utility>

int main(int argc, char* argv[]) {
  Matrix<int> A(8,8);
  Matrix<int> B(8,8);
  int counter = 1;
  for (size_t i = 0; i < B.rows; i++) {
    for (size_t j = 0; j < B.cols; j++) {
      B(i,j) = counter;
      counter++;
    }
  }
  std::cout << "partitioning a matrix" << std::endl;
  PartitionView<int> asd(B,4,4,std::make_pair(4,4));
  asd.print_contents();
  std::cout << "recursive partitioning test" << std::endl;
  PartitionView<int> top(B,8,8,std::make_pair(0,0));
  B.print_contents();
  std::cout << std::endl;
  for (size_t i = 0; i < B.rows-1; i++) {
    top = PartitionView<int>(top,7-i,7-i,std::make_pair(1,1));
    top.print_contents();
    matrix_sum(B, B, B);
    B.print_contents();
    
    std::cout << std::endl;
  }
	return 0;
}
