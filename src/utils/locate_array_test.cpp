#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "locate_array.h"

int main() {
  std::vector<double> vals{-10, 0, 1, 3.3, 9, 10, 11, 20};

  locate_array l;
  std::cout << "------- 20 zeros -------" << std::endl;
  l.init(20);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- lin: start 1, stop 10, step 1 -------" << std::endl;
  l.init(1., 10., 1.);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- log: start 1, stop 10, step 0.2 -------" << std::endl;
  l.log_init(1, 10, 0.2);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- lin: start 1, stop 10, 10 points -------" << std::endl;
  l.init(1., 10., 10);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- flex, array -2 to 7 -------" << std::endl;
  std::vector<double> x;
  x.resize(7);
  for (int i = 0; i < 7; i++) x[i] = i;
  l.init(x, -2.);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- flex, array -1 to 4 -------" << std::endl;
  double* a = &(x[0]);
  l.init(a, 5., -1.);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- none, start stop -------" << std::endl;
  l.init(0., 0., 0.);
  l.print();
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- lin, two elems ------" << std::endl;
  l.init(10., 20., 2);
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }

  std::cout << "------- flex, array 1e13 to 1e14 -------" << std::endl;
  l.init(0., 9.07e13, 1.234e12);
  l.print();
  std::cerr << "3.3e13 " << l.locate(3.3e13) << std::endl;
  std::cerr << "3.4e13 " << l.locate(3.4e13) << std::endl;
  for (auto i_val = vals.begin(); i_val != vals.end(); i_val++) {
    std::cerr << *i_val << " " << l.locate(*i_val) << std::endl;
  }
}
