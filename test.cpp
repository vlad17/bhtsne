// compile with
// g++ -std=c++11 -o test test.cpp
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double randn() {
  double x, y, radius;
  do {
    x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    radius = (x * x) + (y * y);
  } while((radius >= 1.0) || (radius == 0.0));
  radius = sqrt(-2 * log(radius) / radius);
  x *= radius;
  y *= radius;
  return x;
}


int main() {
  srand(1324);
  std::ofstream out("/tmp/data1.dat");
  for (int i = 0; i < 1000; ++i) {
    for (int j = 0; j < 10; ++j) {
      out << randn() << '\t';
    }
    out << '\n';
  }
}
