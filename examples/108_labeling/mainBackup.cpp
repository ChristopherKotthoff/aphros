// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include <dump/hdf.h>
#include <func/init.h>
#include <parse/argparse.h>
#include <parse/vars.h>
#include <util/distr.h>
#include <util/filesystem.h>
#include <util/format.h>
#include <util/vof.h>

#include "distri_CCL.hpp"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

static int called_init_ = 0;
static int called_first_dump_ = 0;
static int called_second_dump_ = 0;
static int called_run_ = 0;
static int called_main_ = 0;


void Init(
    GRange<size_t> layers, Multi<FieldCell<Scal>>& fcu, std::string prefix,
    const Vars& var, M& m) {
  //std::cout << "in init" << std::endl;
  called_init_++;
  auto sem = m.GetSem();
  struct {
    Vars var;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    fcu.Reinit(layers, m, 0);
  }
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    if (sem()) {
      t.var.String.Set("init_vf", "list");
      t.var.String.Set("list_path", var.String["init_" + prefix + sl]);
      t.var.Int.Set("dim", m.GetGlobalSize()[2] == 1 ? 2 : 3);
      t.var.Int.Set("list_ls", 1);
    }
    if (sem.Nested()) {
      InitVf(fcu[l], t.var, m, false);
    }
  }
}

void Dump(
    const FieldCell<Scal>& fcu, std::string fieldname, std::string filename,
    M& m) {
  //std::cout << "in first dump" << std::endl;
  called_first_dump_++;
  auto sem = m.GetSem();
  const auto path = filename + ".h5";
  if (sem.Nested()) {
    Hdf<M>::Write(fcu, path, m);
  }
  if (sem() && m.IsRoot()) {
    Hdf<M>::WriteXmf(util::SplitExt(path)[0] + ".xmf", fieldname, path, m);
  }
}

void Dump(
    GRange<size_t> layers, const Multi<FieldCell<Scal>>& fcu,
    std::string fieldname, std::string filename, M& m) {
  //std::cout << "in second dump" << std::endl;
  called_second_dump_++;
  auto sem = m.GetSem();
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    if (sem.Nested()) {
      Dump(fcu[l], fieldname, filename + sl, m);
    }
  }
}

void Run(M& m, Vars& var) {
  //std::cout << "in run" << std::endl;
  called_run_++;
  auto sem = m.GetSem();
  struct {
    GRange<size_t> layers;
    Multi<FieldCell<Scal>> fcvf; // volume fraction
    Multi<FieldCell<Scal>> fccl; // color
    MapEmbed<BCond<Scal>> mebc; // boundary conditions for volume fraction
    FieldCell<Scal> fcvf_sum; // sum of volume fractions from all layers
    FieldCell<Scal> fccl_sum; // colors from all layers
  } * ctx(sem);
  constexpr Scal kClNone = -1;
  auto& t = *ctx;
  if (sem()) {
    t.layers = GRange<size_t>(var.Int["layers"]);
    t.fccl.Reinit(t.layers, m, kClNone);
  }
  if (sem.Nested()) {
    Init(t.layers, t.fcvf, "vf", var, m);
  }
  if (sem()) {
    // Clear color in cells with zero volume fraction
    for (auto l : t.layers) {
      for (auto c : m.AllCells()) {
        t.fccl[l][c] = (t.fcvf[l][c] > 0 ? l : kClNone);
      }
    }
  }
  if (sem.Nested()) {
    const bool verbose = false;
    const bool reduce = true;
    const bool unionfind = false;
    const bool grid = false;
    UVof<M>().Recolor(
        t.layers, t.fcvf, t.fccl, t.fccl, 0, Vect(0), 1e10, t.mebc, verbose,
        unionfind, reduce, grid, m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fcvf, "vf", "vf", m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fccl, "cl", "cl", m);
  }
  if (sem()) {
    // Collect volume fractions and colors from all layers in one field.
    t.fcvf_sum.Reinit(m, 0);
    t.fccl_sum.Reinit(m, kClNone);
    for (auto c : m.Cells()) {
      for (auto l : t.layers) {
        t.fcvf_sum[c] += t.fcvf[l][c];
        if (t.fcvf[l][c] > 0) {
          t.fccl_sum[c] = t.fccl[l][c];
        }
      }
    }
  }
  if (sem.Nested()) {
    Dump(t.fcvf_sum, "vf", "vf_sum", m);
  }
  if (sem.Nested()) {
    Dump(t.fccl_sum, "cl", "cl_sum", m);
  }
}

int main(int argc, const char** argv) {
  //std::cout << "in main" << std::endl;
  called_main_++;

  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser(
      "Example for connected component labeling", mpi.IsRoot());
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<std::string>("config", "a.conf")
      .Help("Path to configuration file");
  parser.AddVariable<int>("--nx", 128).Help("Mesh size in the x-direction");
  parser.AddVariable<int>("--ny", 128).Help("Mesh size in the y-direction");
  parser.AddVariable<int>("--nz", 128).Help("Mesh size in the z-direction");
  parser.AddVariable<int>("--bs", 32).Help("Block size");
  parser.AddVariable<int>("--layers", 4).Help("Number of layers");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::stringstream conf;

  MIdx meshsize(args.Int["nx"], args.Int["ny"], args.Int["nz"]);
  const int bs = args.Int["bs"];
  MIdx blocksize(bs, bs, args.Int["nz"] == 1 ? 1 : bs);
  Subdomains<MIdx> sub(meshsize, blocksize, mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';

  const auto configpath = args.String["config"];
  if (!configpath.empty()) {
    conf << "include " + configpath + '\n';
  }
  conf << "set string backend native\n";
  conf << "set double extent 1\n";
  conf << util::Format("set int layers {}\n", args.Int["layers"]);
  conf << args.String["extra"] << '\n';

  auto temp = RunMpiBasicString<M>(mpi, Run, conf.str());

  std::cout << "amount called_init_ = " << called_init_ << std::endl; 
  std::cout << "amount called_first_dump_ = " << called_first_dump_ << std::endl; 
  std::cout << "amount called_second_dump_ = " << called_second_dump_ << std::endl; 
  std::cout << "amount called_run_ = " << called_run_ << std::endl; 
  std::cout << "amount called_main_ = " << called_main_ << std::endl; 

  return temp;
}
