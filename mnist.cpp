// compile with
// g++ -std=c++11 -o mnist mnist.cpp
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

template<typename Stream, typename T>
void write(Stream& s, T* t) {
  s.write(reinterpret_cast<char*>(t), sizeof *t);
}

int main() {
  srand(1324);
  system("wget -O /tmp/mnist.idx3.gz "
         "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz");
  system("gunzip /tmp/mnist.idx3.gz");
  FILE* mnist = fopen("/tmp/mnist.idx3", "rb");
  fseek(mnist, 4 + 4 + 4 + 4, 0);
  std::ofstream out("/tmp/mnist.dat", std::ios::binary);
  int n = 10000, d = 28 * 28, no_dims = 2;
  double perplexity = 30, theta = 0.5;
  write(out, &n);
  write(out, &d);
  write(out, &theta);
  write(out, &perplexity);
  write(out, &no_dims);
  for (int i = 0; i < 10000; ++i) {
    for (int j = 0; j < 28 * 28; ++j) {
      unsigned char c = fgetc(mnist);
      write(out, &c);
    }
  }
}
