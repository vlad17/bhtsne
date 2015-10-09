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

#ifndef TSNE_H
#define TSNE_H

#include <tuple>
#include <vector>

static inline double sign(double x) {
  return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0));
}

// Methods may throw if data or parameters are ill-formed. The class
// provides strong exception safety. Throws a std::invalid_argument exception.
// Class i
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
class TSNE {
public:
  static const size_t D = std::tuple_size<Tuple>::value;

  // Runs t-SNE on an N-length array of high-dimensional tuples X,
  // returning a flattened two-dimensional (and thus no_dims*N-length) vector of
  // doubles.
  //
  // Note that because this is type-agnostic the user is responsible for
  // making sure distances are normalized.
  //
  // Uses rand(), so can be seeded.
  std::vector<double> run(std::vector<Tuple> X, int no_dims, double theta,
                          double perplexity);

  // Loads a list of tuples according to this type's template from a TSV
  // file.
  std::vector<Tuple> load_data(
    const std::string& filename, double* theta, double* perplexity,
    int* no_dims, int* rand_seed);

  // Saves a flattened double-matrix with no_dims dimensions in TSV format.
  void save_data(const std::string& filename, int no_dims,
                 const std::vector<double>& Y);
private:
  int N;
  void zeroMean(std::vector<double>* Y, int no_dims);
  // Corresponds to an approximation for row and column of the
  // P matrix in the paper.
  struct HiDimSimilarityDist {
    std::vector<unsigned> row_;
    std::vector<unsigned> col_;
    std::vector<double> val_;
  };
  float seconds(clock_t start, clock_t end);
  void bail(const std::string& s);
  HiDimSimilarityDist symmetrizeMatrix(HiDimSimilarityDist P, int N);
  std::vector<double> computeGradient(
      HiDimSimilarityDist* P, std::vector<double>* Y, int no_dims,
      double theta);
  double evaluateError(HiDimSimilarityDist* P, std::vector<double>* Y,
                       int no_dims, double theta);
  HiDimSimilarityDist computeGaussianPerplexity(
      const std::vector<Tuple>& X, double perplexity, int K);
  typedef std::pair<unsigned, Tuple> IndexedTuple;
  double randn();
};

#include "tsne.tpp"

#endif
