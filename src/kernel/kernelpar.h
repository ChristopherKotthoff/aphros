#pragma once

#include <memory>

#include "geom/vect.h"
#include "linear/linear.h"
#include "kernel.h"
#include "parse/vars.h"
#include "kernelmesh.h"

// Kernel aware of mesh with structure Par_ passed to constructor.
template <class M_, class Par_>
class KernelMeshPar : public KernelMesh<M> {
 public:
  using P = KernelMesh<M>; // parent
  using M = M_;
  using Par = Par_;
  static constexpr size_t dim = M::dim;

  KernelMeshPar(Vars& var, const MyBlockInfo& bi, Par& par)
      : KernelMesh<M>(var, bi)
      , par(par) {}
  void Run() override;

 protected:
  using P::var;
  using P::m;
  using P::IsRoot;
  using P::IsLead;
  Par& par_;
};

// Factory for KernelMeshPar.
// M_: mesh
// K_: kernel derived from KernelMeshPar with defined Par
template <class M_, class K_>
class KernelMeshParFactory : public KernelMeshFactory<M_> {
 public:
  using M = M_;
  using K = K_;
  using Par = typename K::Par;
  GKernelMeshFactory(Par& par) : par_(par) {}
  K* Make(Vars& var, const MyBlockInfo& bi) override {
    return new K(var, bi, par_);
  }

 protected:
  Par& par_;
};

