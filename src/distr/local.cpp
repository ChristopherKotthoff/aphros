#include "local.h"

#include "geom/vect.h"
#include "geom/mesh3d.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

template <class KF>
Distr* TryLocal(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Local<KF>(comm, *kfd, par);
  }
  return nullptr;
}

Distr* CreateLocal(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {
  Distr* r = nullptr;
  if (!r) r = TryLocal<KernelMeshFactory<MeshStructured<double, 3>>>(
      comm, kf, par);
  //if (!r) r = Try<KernelMeshFactory<geom3d::MeshStructured<float>>(
  //    comm, kf, par);
  assert(r && "CreateLocal(): KernelFactory not found");
  return r;
}
