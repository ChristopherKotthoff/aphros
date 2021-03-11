// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#include <unistd.h>
#include <fstream>
#include <sstream>

#include "sysinfo.h"
#include "system.h"
#include "logger.h"
#include "util/mpi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if USEFLAG(AMGX)
#include <cuda_runtime.h>
#endif

namespace sysinfo {

size_t GetMem() {
  return SystemGetMem();
}

bool HasHyperthreads() {
  return SystemHasHyperthreads();
}

std::string GetHostname() {
  const size_t kMaxLength = 4096;
  char buf[kMaxLength];
  int err = gethostname(buf, kMaxLength);
  fassert_equal(err, 0);
  return buf;
}

Info GetInfo(InfoSelect select) {
  Info info;
#if USEFLAG(MPI)
  if (select.mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &info.comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &info.comm_rank);
  }
#endif

#ifdef _OPENMP
  if (select.openmp) {
    info.omp_num_threads = omp_get_num_threads();
    info.omp_max_threads = omp_get_max_threads();
  }
#endif

#if USEFLAG(AMGX)
  if (select.cuda) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    info.cuda_enabled = true;
    info.cuda_uuid = *(uint16_t*)(&prop.uuid);
    info.cuda_mem = prop.totalGlobalMem;
  }
#endif

  info.hostname = GetHostname();
  return info;
}

} // namespace sysinfo
