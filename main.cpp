/* Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of
 *    Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Amended by Vladimir Feinberg on October 2015 to allow for generic types. */

#include <iostream>
#include <string>

#include "tsne.hpp"
#include "tuple_iter.hpp"

typedef std::tuple<double, double, double, double, double,
                   double, double, double, double, double> Double10;

void line_assert(const std::string& expect, bool b) {
  if (!b) {
    throw std::invalid_argument("Bad line parse: " + expect);
  }
}

std::istream& operator>>(std::istream& in, Double10& d10) {
  int last = EOF;
  auto readtok = [&](double& x) {
    in >> x;
    last = in.get();
    line_assert("expected tab or EOF", last == EOF || last == '\t');
  };
  for_each(d10, readtok);
  return in;
}

double L2(const Double10& d1, const Double10& d2) {
  double sq = 0;
  for_each2(d1, d2, [&](double d1, double d2) { sq += pow(d1 - d2, 2); });
  return sq;
}

struct Img28x28 {
  unsigned char pixels_[28 * 28];
};

double L2(const Img28x28& i1, const Img28x28& i2) {
  double ret = 0;
  for (int i = 0; i < 28 * 28; ++i) {
    double x = i1.pixels_[i] - i2.pixels_[i];
    ret += x * x;
  }
  return ret;
}

namespace std {
  template<> class tuple_size<Img28x28> {
  public:
    static const size_t value = 28 * 28;
  };
}

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " infile outfile"
              << std::endl;
    return 1;
  }

  const char* infile = argv[1];
  const char* outfile = argv[2];
  //  auto no_dims = stoi(argv[3]);
  //  auto theta = stod(argv[4]);
  //  auto perplexity = stod(argv[5]);

  //  auto rand_seed = stoi(argv[6]);
  //  if (rand_seed == 0) rand_seed = time(NULL);
  //srand(rand_seed);
  srand(0);
  //  std::cout << "Using random seed " << rand_seed << std::endl;

  TSNE<Img28x28, L2> tsne;
  int no_dims = 0, rand_seed = 0;
  double theta = 0, perplexity = 0;

  try {
    std::cout << "Loading 28x28-tuples from " << infile << "..." << std::endl;
    auto X = tsne.load_data(infile, &theta, &perplexity, &no_dims, &rand_seed);

    // Normalize
    /*
    Img28x28 mean;
    for (auto& x : X)
      for_each2(mean, x, [](double& m, double& d) { m += d; });
    for_each(mean, [&](double& m) { m /= X.size(); });
    for (auto& x : X)
      for_each2(x, mean, [](double& d, double& m) { d -= m; });

    double max_X = 0.;
    for (const auto& x : X)
      for_each(x, [&](double d) { max_X = std::max(max_X, fabs(d)); });
    for (auto& x : X)
      for_each(x, [&](double& d) { d /= max_X; });
    */
    auto Y = tsne.run(std::move(X), no_dims, theta, perplexity);
    tsne.save_data(outfile, no_dims, Y);

  } catch (const std::invalid_argument& e) {
    std::cerr << "Error: invalid argument: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
