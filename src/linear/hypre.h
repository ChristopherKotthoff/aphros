#pragma once

#include <vector>
#include <mpi.h>
#include <memory>
#include <string>

// FIXME: Convention *_ for private variables ignored

class Hypre {
 public:
  using MIdx = std::vector<int>;
  using Scal = double;
  struct Block { // linear system ax=r
    MIdx l; // lower corner
    MIdx u; // upper corner
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * st.size()
    std::vector<Scal>* r; // rhs of size n
    std::vector<Scal>* x; // solution and initial guess of size n
  };

  // bb: blocks
  // gs: global size
  // per: periodic in each direction
  // tol: tolerance
  // print: print level
  // solver: solver name
  // maxiter: maximum number of iterations
  Hypre(MPI_Comm comm, const std::vector<Block>& bb,
        MIdx gs, std::vector<bool> per);
  Hypre() = delete;
  Hypre(const Hypre&) = delete;

  // Assembles matrix and vectors from bb
  void Update();
  // Solves system and puts result to x
  void Solve(Scal tol, int print, std::string solver, int maxiter);
  // Returns relative residual norm from last Solve()
  Scal GetResidual();
  // Returns the number of iterations from last Solve()
  int GetIter();
  ~Hypre();

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
