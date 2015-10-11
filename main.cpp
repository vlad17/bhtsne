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
#include "gdelt.hpp"

template<std::size_t S>
struct read_cell { void operator()(GDELTMini& m, std::istream& in) {
  in >> std::get<S>(m);
  int last = in.get();
  if (last != EOF && last != '\t')
    throw std::invalid_argument("expected tab or EOF");
}};

std::istream& operator>>(std::istream& in, GDELTMini& gdm) {
  std::string line;
  std::getline(in, line);
  std::stringstream line_stream(line);
  for_eachi<read_cell>(gdm, line_stream);
  return in;
}

template<std::size_t I>
struct get_dist { void operator() (GDELTMini& m1, GDELTMini& m2, double& sum) {
  sum += std::get<I>(m1).distance(std::get<I>(m2));
}
};

double Rho(const GDELTMini& m1, const GDELTMini& m2) {
  double dist = 0;
  for_eachi<get_dist>(const_cast<GDELTMini&>(m1),
                      const_cast<GDELTMini&>(m2), dist);
  return dist;
}

int main(int argc, const char* argv[]) {
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " infile outfile no_dims theta perplexity rand_seed"
              << std::endl;
    return 1;
  }

  const char* infile = argv[1];
  const char* outfile = argv[2];
  auto no_dims = stoi(argv[3]);
  auto theta = stod(argv[4]);
  auto perplexity = stod(argv[5]);

  auto rand_seed = stoi(argv[6]);
  if (rand_seed == 0) rand_seed = time(NULL);
  srand(rand_seed);
  std::cout << "Using random seed " << rand_seed << std::endl;

  TSNE<GDELTMini, Rho> tsne;

  try {
    std::cout << "Loading " << std::tuple_size<GDELTMini>::value
              << "-tuples from " << infile << "..." << std::endl;
    auto X = tsne.load_data(infile);

    // Normalize
    /*
    Double10 mean;
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
