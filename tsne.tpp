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

#ifndef TSNE_TPP
#define TSNE_TPP

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "vptree.h"
#include "sptree.h"
#include "tsne.hpp" // header guards prevent infinite include

template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
float TSNE<Tuple, Distance>::seconds(clock_t start, clock_t end) {
  return static_cast<float>((end - start) / CLOCKS_PER_SEC);
}

template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
void TSNE<Tuple, Distance>::bail(const std::string& s) {
  throw std::invalid_argument(s);
}

template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
std::vector<double> TSNE<Tuple, Distance>::run(
    std::vector<Tuple> X, const int no_dims, double theta, double perplexity) {
  N = X.size();

  // Determine whether we are using an exact algorithm
  if (N - 1 < 3 * perplexity) {
    bail("Perplexity too large for the number of data points.");
  }

  std::cout << "Using no_dims = " << no_dims
            << ", perplexity = " << perplexity
            << ", and theta = " << theta << std::endl;

  if (theta <= 0.0) {
    bail("theta <= 0.0, too exact.");
  }

  // Set learning parameters
  clock_t start = clock(), end;
  int max_iter = 1000, stop_lying_iter = 250, mom_switch_iter = 250;
  double momentum = .5, final_momentum = .8;
  double eta = 200.0;

  // Allocate intermediate and result vectors
  std::vector<double> Y(N * no_dims, 0);
  std::vector<double> uY(N * no_dims, 0.0);
  std::vector<double> gains(N * no_dims, 1.0);

  // Compute asymmetric pairwise input similarities

  auto P = computeGaussianPerplexity(X, perplexity, static_cast<int>(3 * perplexity));

  // Symmetrize input similarities
  P = symmetrizeMatrix(std::move(P), N);
  double sum_P = .0;
  for (int i = 0; i < P.row_[N]; i++) sum_P += P.val_[i];
  for (int i = 0; i < P.row_[N]; i++) P.val_[i] /= sum_P;

  end = clock();

  std::cout << "Input similarities computed in "
            << seconds(start, end)
            << " seconds (sparsity = "
            << static_cast<double>(P.row_[N]) / (N * N)
            << ")!" << std::endl;
  std::cout << "Learning embedding..." << std::endl;

  // Lie about the P-values
  for (int i = 0; i < P.row_[N]; i++) P.val_[i] *= 12.0;

  // Initialize solution (randomly)
  for (int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;

  // Perform main training loop
  start = clock();
  auto very_start = start;
  for (int iter = 0; iter < max_iter; iter++) {
    // Compute (approximate) gradient
    auto dY = computeGradient(&P, &Y, no_dims, theta);

    // Update gains
    for (int i = 0; i < N * no_dims; i++)
      gains[i] = (sign(dY[i]) != sign(uY[i])) ?
        (gains[i] + .2) : (gains[i] * .8);
    for (int i = 0; i < N * no_dims; i++)
      if (gains[i] < .01) gains[i] = .01;

    // Perform gradient update (with momentum and gains)
    for (int i = 0; i < N * no_dims; i++)
      uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
    for (int i = 0; i < N * no_dims; i++) Y[i] = Y[i] + uY[i];

    // Make solution zero-mean
    zeroMean(&Y, no_dims);

    // Stop lying about the P-values after a while, and switch momentum
    if (iter == stop_lying_iter) {
      for (int i = 0; i < P.row_[N]; i++) P.val_[i] /= 12.0;
    }
    if (iter == mom_switch_iter) momentum = final_momentum;

    // Print out progress
    if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
      end = clock();
      // Doing approximate computation here!
      double C = evaluateError(&P, &Y, no_dims, theta);
      std::cout << "Iteration " << iter << ": error is " << C;
      if (iter % 50 == 0) {
        std::cout << " (50 iterations in " << seconds(start, end) << "s)";
        start = clock();
      }
      std::cout << std::endl;
    }
  }

  end = clock();
  std::cout << "Fitting performed in " << seconds(very_start, end)
            << " seconds." << std::endl;
  return Y;
}


// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
vector<double> TSNE<Tuple, Distance>::computeGradient(
    HiDimSimilarityDist* P, std::vector<double>* Y, int no_dims, double theta) {
  // Construct space-partitioning tree on current map
  SPTree tree(no_dims, &(*Y)[0], N);

  // Compute all terms required for t-SNE gradient
  double sum_Q = 0.0;
  vector<double> pos_f(N * no_dims, 0);
  vector<double> neg_f(N * no_dims, 0);
  tree.computeEdgeForces(&P->row_[0], &P->col_[0], &P->val_[0], N, &pos_f[0]);

  for (int n = 0; n < N; n++)
    tree.computeNonEdgeForces(n, theta, &neg_f[n * no_dims], &sum_Q);

  // Compute final t-SNE gradient
  vector<double> dC(N * no_dims);
  for (int i = 0; i < N * no_dims; i++) {
    dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
  }
  return dC;
}

// Evaluate t-SNE cost function (approximately)
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
double TSNE<Tuple, Distance>::evaluateError(
    HiDimSimilarityDist* P, std::vector<double>* Yv, int no_dims, double theta) {
  auto& row_P = P->row_;
  auto& col_P = P->col_;
  auto& val_P = P->val_;
  auto* Y = &(*Yv)[0];


    // Get estimate of normalization term
    SPTree* tree = new SPTree(no_dims, Y, N);
    double* buff = (double*) calloc(no_dims, sizeof(double));
    double sum_Q = .0;
    for (int n = 0; n < N; n++) tree->computeNonEdgeForces(n, theta, buff, &sum_Q);

    // Loop over all edges to compute t-SNE error
    int ind1, ind2;
    double C = .0, Q;
    for (int n = 0; n < N; n++) {
        ind1 = n * no_dims;
        for (int i = row_P[n]; i < row_P[n + 1]; i++) {
            Q = .0;
            ind2 = col_P[i] * no_dims;
            for (int d = 0; d < no_dims; d++) buff[d]  = Y[ind1 + d];
            for (int d = 0; d < no_dims; d++) buff[d] -= Y[ind2 + d];
            for (int d = 0; d < no_dims; d++) Q += buff[d] * buff[d];
            Q = (1.0 / (1.0 + Q)) / sum_Q;
            C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
        }
    }

    // Clean up memory
    free(buff);
    delete tree;
    return C;
}

namespace internal {
  // TODO figure this out
template<typename A, typename B, double (*Distance)(const B&, const B&)>
struct SndDistance {
  typedef std::pair<A, B> P;
  double operator()(const P& p1, const P& p2) const {
      return Distance(p1.second, p2.second);
  };
};
}  // namespace internal

// Compute input similarities with a fixed perplexity using ball trees (this
// function allocates memory another function should free)
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
typename TSNE<Tuple, Distance>::HiDimSimilarityDist
TSNE<Tuple, Distance>::computeGaussianPerplexity(
    const vector<Tuple>& X, double perplexity, int K) {
  if (perplexity > K) bail("Perplexity should be lower than K!\n");

  std::vector<unsigned> row_P(N + 1, 0);
  std::vector<unsigned> col_P(N * K, 0);
  std::vector<double> val_P(N * K, 0);
  std::vector<double> cur_P((N - 1) * K, 0);
  std::vector<IndexedTuple> indexed_X(N);
  for (int i = 0; i < N; ++i) indexed_X[i] = std::make_pair(i, X[i]);

  row_P[0] = 0;
  for (int n = 0; n < N; n++) row_P[n + 1] = row_P[n] + static_cast<unsigned>(K);

  // Build ball tree on data set
  VpTree<IndexedTuple, internal::SndDistance<unsigned, Tuple, Distance>> tree;
  tree.create(std::move(indexed_X));

  // Loop over all points to find nearest neighbors
  std::cout << "  Building tree..." << std::endl;
  vector<IndexedTuple> indices;
  vector<double> distances;
  for (int n = 0; n < N; n++) {
    if (n % 10000 == 0) {
      std::cout << "    - point " << n << " of " << N << std::endl;
    }

    // Find nearest neighbors
    indices.clear();
    distances.clear();
    tree.search(std::make_pair(n, X[n]), K + 1, &indices, &distances);

    // Initialize some variables for binary search
    bool found = false;
    double beta = 1.0;
    double min_beta = -DBL_MAX;
    double max_beta =  DBL_MAX;
    double tol = 1e-5;

    // Iterate until we found a good perplexity
    int iter = 0; double sum_P;
    while (!found && iter < 200) {
      // Compute Gaussian kernel row
      for (int m = 0; m < K; m++) cur_P[m] = exp(-beta * distances[m + 1]);

      // Compute entropy of current row
      sum_P = DBL_MIN;
      for (int m = 0; m < K; m++) sum_P += cur_P[m];
      double H = .0;
      for (int m = 0; m < K; m++) H += beta * (distances[m + 1] * cur_P[m]);
      H = (H / sum_P) + log(sum_P);

      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - log(perplexity);
      if (Hdiff < tol && -Hdiff < tol) {
        found = true;
      }
      else {
        if (Hdiff > 0) {
          min_beta = beta;
          if (max_beta == DBL_MAX || max_beta == -DBL_MAX)
            beta *= 2.0;
          else
            beta = (beta + max_beta) / 2.0;
        }
        else {
          max_beta = beta;
          if (min_beta == -DBL_MAX || min_beta == DBL_MAX)
            beta /= 2.0;
          else
            beta = (beta + min_beta) / 2.0;
        }
      }

      // Update iteration counter
      iter++;
    }

    // Row-normalize current row of P and store in matrix
    for (int m = 0; m < K; m++) cur_P[m] /= sum_P;
    for (int m = 0; m < K; m++) {
      col_P[row_P[n] + m] = static_cast<unsigned>(indices[m + 1].first);
      val_P[row_P[n] + m] = cur_P[m];
    }
  }

  return {std::move(row_P), std::move(col_P), std::move(val_P)};
}


// Symmetrizes a sparse matrix
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
typename TSNE<Tuple, Distance>::HiDimSimilarityDist
TSNE<Tuple, Distance>::symmetrizeMatrix(HiDimSimilarityDist P, int N) {
  // Count number of elements and row counts of symmetric matrix
  std::vector<int> row_counts(N, 0);
  for (int n = 0; n < N; n++) {
    for (int i = P.row_[n]; i < P.row_[n + 1]; i++) {

      // Check whether element (P.col_[i], n) is present
      bool present = false;
      for (int m = P.row_[P.col_[i]]; m < P.row_[P.col_[i] + 1]; m++) {
        if (P.col_[m] == n) present = true;
      }
      if (present) row_counts[n]++;
      else {
        row_counts[n]++;
        row_counts[P.col_[i]]++;
      }
    }
  }

  int no_elem = 0;
  for (int n = 0; n < N; n++) no_elem += row_counts[n];

  // Allocate memory for symmetrized matrix
  std::vector<unsigned> sym_row_P(N + 1, 0);
  std::vector<unsigned> sym_col_P(no_elem, 0);
  std::vector<double> sym_val_P(no_elem, 0);

  // Construct new row indices for symmetric matrix
  sym_row_P[0] = 0;
  for (int n = 0; n < N; n++)
    sym_row_P[n + 1] = sym_row_P[n] + static_cast<unsigned>(row_counts[n]);

  // Fill the result matrix
  std::vector<int> offset(N);
  for (int n = 0; n < N; n++) {
    for (unsigned int i = P.row_[n]; i < P.row_[n + 1]; i++) {
      // considering element(n, P.col_[i])

      // Check whether element (P.col_[i], n) is present
      bool present = false;
      for (unsigned int m = P.row_[P.col_[i]]; m < P.row_[P.col_[i] + 1]; m++) {
        if (P.col_[m] == n) {
          present = true;
          if (n <= P.col_[i]) {
            // make sure we do not add elements twice
            sym_col_P[sym_row_P[n]        + offset[n]]        = P.col_[i];
            sym_col_P[sym_row_P[P.col_[i]] + offset[P.col_[i]]] = n;
            sym_val_P[sym_row_P[n]        + offset[n]]        = P.val_[i] + P.val_[m];
            sym_val_P[sym_row_P[P.col_[i]] + offset[P.col_[i]]] = P.val_[i] + P.val_[m];
          }
        }
      }

      // If (P.col_[i], n) is not present, there is no addition involved
      if (!present) {
        sym_col_P[sym_row_P[n]         + offset[n]]         = P.col_[i];
        sym_col_P[sym_row_P[P.col_[i]] + offset[P.col_[i]]] = n;
        sym_val_P[sym_row_P[n]         + offset[n]]         = P.val_[i];
        sym_val_P[sym_row_P[P.col_[i]] + offset[P.col_[i]]] = P.val_[i];
      }

      // Update offsets
      if (!present || (present && n <= P.col_[i])) {
        offset[n]++;
        if (P.col_[i] != n) offset[P.col_[i]]++;
      }
    }
  }

  // Divide the result by two
  for (int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

  return {std::move(sym_row_P), std::move(sym_col_P), std::move(sym_val_P)};
}

// Makes low-dimensional data zero-mean
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
void TSNE<Tuple, Distance>::zeroMean(std::vector<double>* Y, int no_dims) {
  // Compute data mean
  const int N = Y->size() / no_dims;
  std::vector<double> mean(no_dims, 0);
  int nD = 0;
  for (int n = 0; n < N; n++) {
    for (int d = 0; d < no_dims; d++) {
      mean[d] += (*Y)[nD + d];
    }
    nD += no_dims;
  }
  for (int d = 0; d < no_dims; d++) {
    mean[d] /= static_cast<double>(N);
  }

  // Subtract data mean
  nD = 0;
  for (int n = 0; n < N; n++) {
    for (int d = 0; d < no_dims; d++) {
      (*Y)[nD + d] -= mean[d];
    }
    nD += no_dims;
  }
}

// Generates a Gaussian random number
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
double TSNE<Tuple, Distance>::randn() {
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

// Function that loads data from a t-SNE file
// Note: this function does a malloc that should be freed elsewhere
template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
std::vector<Tuple> TSNE<Tuple, Distance>::load_data(const std::string& filename) {
  std::ifstream in(filename);
  std::vector<Tuple> data;
  std::string line;
  while (std::getline(in, line)) {
    std::stringstream sstr(line);
    data.emplace_back();
    in >> data.back();
  }
  std::cout << "Read the " << data.size() << " x " << D
            << " data matrix successfully!" << std::endl;
  return data;
}

template<typename Tuple, double (*Distance)(const Tuple&, const Tuple&)>
void TSNE<Tuple, Distance>::save_data(const std::string& filename,
                                      const int no_dims,
                                      const std::vector<double>& Y) {
  if (Y.size() % no_dims != 0) {
    bail("Flattened matrix length " + std::to_string(Y.size())
         + " incompatible with row length " + std::to_string(no_dims));
  }

  int nrows = Y.size() / no_dims;
  std::cout << "Saving " << nrows << " x " << no_dims
            << " data matrix to " << filename << ", flat binary fmt"
            << std::endl;

  std::ofstream out(filename, std::ios::binary);
  int ctr = 0;
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < no_dims; ++j) {
      out.write(reinterpret_cast<const char*>(&Y[ctr++]), sizeof (double));
    }
  }
}

#endif /* TSNE_TPP */
