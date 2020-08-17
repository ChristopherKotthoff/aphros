// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

#include "debug/isnan.h"
#include "dump/dumper.h"
#include "dump/hdf.h"
#include "func/init.h"
#include "func/init_bc.h"
#include "func/init_contang.h"
#include "geom/mesh.h"
#include "kernelmeshpar.h"
#include "parse/curv.h"
#include "parse/parser.h"
#include "parse/proj.h"
#include "parse/simple.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/advection.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/curv.h"
#include "solver/embed.h"
#include "solver/fluid_dummy.h"
#include "solver/multi.h"
#include "solver/normal.h"
#include "solver/particles.h"
#include "solver/proj.h"
#include "solver/reconst.h"
#include "solver/simple.h"
#include "solver/solver.h"
#include "solver/tracer.h"
#include "solver/vof.h"
#include "solver/vofm.h"
#include "util/convdiff.h"
#include "util/events.h"
#include "util/filesystem.h"
#include "util/hydro.h"
#include "util/metrics.h"
#include "util/posthook.h"
#include "util/stat.h"

class GPar {};

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    const FieldFace<typename M::Scal>& ffv, const M& m, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  FieldCell<Scal> fcdiv(m, 0);
  for (auto c : eb.Cells()) {
    Scal div = 0;
    for (auto q : eb.Nci(c)) {
      div += ffv[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcdiv[c] = div / eb.GetVolume(c);
  }
  return fcdiv;
}

template <class M>
FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcvel,
    const MapEmbed<BCond<typename M::Vect>>& me_vel, Embed<M>& eb) {
  auto& m = eb.GetMesh();
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using UEB = UEmbed<M>;

  std::array<FieldCell<Vect>, 3> grad;
  for (size_t d = 0; d < 3; ++d) {
    grad[d].Reinit(eb, Vect(0));
    const auto mebc = GetScalarCond(me_vel, d, m);
    const FieldCell<Scal> fcu = GetComponent(fcvel, d);
    const FieldEmbed<Scal> ffg = UEB::Gradient(fcu, mebc, eb);
    grad[d] = UEB::AverageGradient(ffg, eb);
  }

  FieldCell<Vect> r(m, Vect(0));
  for (auto c : eb.Cells()) {
    r[c][0] = grad[2][c][1] - grad[1][c][2];
    r[c][1] = grad[0][c][2] - grad[2][c][0];
    r[c][2] = grad[1][c][0] - grad[0][c][1];
  }

  return r;
}

template <class M_>
class Hydro : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>; // parent
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Sem = typename M::Sem;
  using Par = GPar;
  template <class T>
  using Multi = Multi<T>;
  using UEB = UEmbed<M>;
  using ParticlesView = typename ParticlesInterface<M>::ParticlesView;
  static constexpr size_t dim = M::dim;
  friend void StepHook<>(Hydro*);

  // TODO: issue warning if variable in Vars was not used
  // but differs from default (like in CMake)
  Hydro(Vars&, const MyBlockInfo&, Par&);
  void Run() override;
  M& GetMesh() {
    return m;
  }

 protected:
  using P::bi_;
  using P::m;
  using P::var;

 private:
  void Init();
  void InitEmbed();
  void InitTracer(Multi<FieldCell<Scal>>& vfcu);
  void InitTracerFields(Multi<FieldCell<Scal>>& vfcu);
  void SpawnTracer();
  void InitParticles();
  void SpawnParticles(ParticlesView& view);
  void OverwriteBc();
  void InitFluid(const FieldCell<Vect>& fc_vel);
  void InitAdvection(const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl);
  void InitStat();
  Vect CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb);
  Vect CalcViscousDrag(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
      const Embed<M>& eb);
  void DumpFields();
  void Dump(bool force);
  // Calc rho, mu and force based on volume fraction
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();
  void ReportStep();
  void ReportStepAdv();
  void ReportStepTracer();
  void ReportStepParticles();
  void ReportIter();
  // Issue sem.LoopBreak if abort conditions met
  void CheckAbort(Sem& sem, Scal& nabort);
  void StepFluid();
  void StepAdvection();
  void StepTracer();
  void StepParticles();
  void StepBubgen();
  void StepEraseVolumeFraction(std::string prefix, Scal& last_t);
  void StepEraseColor(std::string prefix);

  using ASV = Vof<M>; // advection VOF
  using ASVM = Vofm<M>; // advection multi VOF
  using ASVEB = Vof<Embed<M>>; // advection VOF embed
  using ASVMEB = Vofm<Embed<M>>; // advection multi VOF embed
  using EB = Embed<M>;
  static constexpr Scal kClNone = ASVM::kClNone;

  void UpdateAdvectionPar() {
    if (auto as = dynamic_cast<ASV*>(as_.get())) {
      as->SetPar(ParsePar<ASV>()(var));
    } else if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      as->SetPar(ParsePar<ASVM>()(var));
    } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
      as->SetPar(ParsePar<ASVEB>()(var));
    } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      as->SetPar(ParsePar<ASVMEB>()(var));
    }
  }
  // Surface tension time step
  Scal GetStDt() {
    const Scal sig = var.Double["sigma"];
    const Scal* cflst = var.Double.Find("cflst");
    if (cflst && sig != 0.) {
      Scal pi = M_PI;
      Scal h3 = m.GetVolume(IdxCell(0));
      Scal r1 = var.Double["rho1"];
      Scal r2 = var.Double["rho2"];
      return (*cflst) * std::sqrt(h3 * (r1 + r2) / (4. * pi * std::abs(sig)));
    }
    return std::numeric_limits<Scal>::max();
  }
  // Viscosity time step
  Scal GetVisDt() {
    const Scal rho1 = var.Double["rho1"];
    const Scal rho2 = var.Double["rho2"];
    const Scal mu1 = var.Double["mu1"];
    const Scal mu2 = var.Double["mu2"];
    const Scal nu1 = mu1 / rho1;
    const Scal nu2 = mu2 / rho2;
    const Scal num = std::max(nu1, nu2);
    const Scal* cflvis = var.Double.Find("cflvis");
    if (cflvis && num != 0.) {
      const Scal h2 = sqr(m.GetCellSize()[0]); // XXX adhoc cubic cell
      return (*cflvis) * h2 / num;
    }
    return std::numeric_limits<Scal>::max();
  }
  void CalcVort() {
    auto& fcv = fs_->GetVelocity();
    if (eb_) {
      fcom_ = GetVort(fcv, fs_->GetVelocityCond(), *eb_);
    } else {
      fcom_ = GetVort(fcv, fs_->GetVelocityCond(), m);
    }
    fcomm_.Reinit(m);
    for (auto c : m.Cells()) {
      fcomm_[c] = fcom_[c].norm();
    }
  }
  FieldCell<Scal> CalcStrain(const FieldCell<Vect> fcvel) const {
    auto& fcv = fcvel;
    auto ffv = UEB::Interpolate(fcv, fs_->GetVelocityCond(), m);

    std::array<FieldCell<Vect>, dim> g; // g[i][c][j] is derivative du_i/dx_j
    for (size_t i = 0; i < dim; ++i) {
      g[i] = Gradient(GetComponent(ffv, i), m);
    }

    FieldCell<Scal> fcs(m, 0);
    int edim = var.Int["dim"];
    for (auto c : m.Cells()) {
      for (int i = 0; i < edim; ++i) {
        for (int j = 0; j < edim; ++j) {
          fcs[c] += sqr(g[i][c][j]) + g[i][c][j] * g[j][c][i];
        }
      }
      fcs[c] *= 0.5;
    }
    return fcs;
  }
  FieldCell<Scal> GetDiv() {
    auto& ffv = fs_->GetVolumeFlux().GetFieldFace();
    if (eb_) {
      return GetDivergence(ffv, m, *eb_);
    }
    FieldCell<Scal> fc(m, 0); // result
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        fc[c] += ffv[f] * m.GetOutwardFactor(c, q);
      }
      fc[c] /= m.GetVolume(c);
    }
    return fc;
  }

  GRange<size_t> layers;
  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source of mixture volume
  FieldCell<Scal> fc_src2_; // source of second phase volume
  FieldCell<Scal> fc_srcm_; // mass source
  FieldCell<Vect> fc_force_; // force
  FieldCell<Scal> fc_dist_; // distance from eb
  FieldCell<Scal> fc_phi_; // distance from eb
  FieldEmbed<Scal> febp_; // balanced force projections

  MapEmbed<BCondAdvection<Scal>> mf_adv_;
  MapEmbed<BCond<Scal>> mebc_vfsm_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_orig_;
  MapEmbed<size_t> me_group_;
  MapEmbed<Scal> me_contang_;
  std::vector<std::string> vdesc_;
  MapCell<std::shared_ptr<CondCell>> mc_cond_;
  MapCell<std::shared_ptr<CondCellFluid>> mc_velcond_;

  std::unique_ptr<Embed<M>> eb_;
  std::unique_ptr<AdvectionSolver<M>> as_;
  std::unique_ptr<FluidSolver<M>> fs_;

  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_

  FieldCell<Scal> fc_sig_; // surface tension sigma
  FieldCell<Scal> fc_contang_; // contact angle

  FieldCell<Vect> fcvm_; // velocity field time_prev // TODO: revise

  FieldCell<Vect> fcom_; // vorticity
  FieldCell<Scal> fcomm_; // vorticity magnitude
  FieldCell<Scal> fc_strain_; // double inner product of strain rate tensor
  Multi<FieldCell<Scal>> fck_; // curvature
  typename PartStrMeshM<M>::Par psm_par_;
  std::unique_ptr<PartStrMeshM<M>> psm_;

  Scal bgt_ = -1.; // bubgen last time
  Scal erasevf_last_t_ = -std::numeric_limits<Scal>::max();
  Scal erasevf2_last_t_ = -std::numeric_limits<Scal>::max();

  struct StatHydro {
    Scal dt = 0; // dt fluid
    Scal dta = 0; // dt advection
    Scal t = 0;
    size_t step = 0;
    size_t iter = 0;
  };
  StatHydro st_;
  std::ofstream fstat_;
  std::unique_ptr<Stat<M>> stat_;
  Dumper dumper_;
  Dumper dmptraj_; // dumper for traj
  Dumper dmptrep_; // dumper for timer report
  std::unique_ptr<Events> events_; // events from var
  SingleTimer timer_;
  std::unique_ptr<TracerInterface<M>> tracer_;
  std::unique_ptr<ParticlesInterface<M>> particles_;
  std::mt19937 randgen_;
  Scal tracer_dt_;
  Scal particles_dt_;
  std::string vf_save_state_path_;
};

template <class M>
void Hydro<M>::InitEmbed() {
  if (var.Int["enable_embed"]) {
    auto sem = m.GetSem("embed");
    struct {
      FieldNode<Scal> fnl;
    } * ctx(sem);
    if (sem("ctor")) {
      eb_.reset(new Embed<M>(m, var.Double["embed_gradlim"]));
      ctx->fnl = UEB::InitEmbed(m, var, m.IsRoot());
      InitEmbedHook(ctx->fnl, var, m);
    }
    if (sem.Nested("smoothen")) {
      SmoothenNode(ctx->fnl, m, var.Int["embed_smoothen_iters"]);
    }
    if (sem.Nested("init")) {
      eb_->Init(ctx->fnl);
    }
  }
}
template <class M>
void Hydro<M>::SpawnTracer() {
  auto& conf = tracer_->GetConf();
  const auto trl = GRange<size_t>(conf.layers);
  const Vect sphere_c(var.Vect["tracer_spawn_sphere_c"]);
  const Scal sphere_r = var.Double["tracer_spawn_sphere_r"];
  const size_t dim = m.GetEdim();
  auto vfcu = tracer_->GetVolumeFraction();
  for (auto l : trl) {
    const std::string prefix = "tracer" + std::to_string(l);
    auto k = var.Double[prefix + "_factor"];
    for (auto c : m.AllCells()) {
      const auto xc = m.GetCenter(c);
      Vect dx = xc - sphere_c;
      if (dim == 2) {
        dx[2] = 0;
      }
      if (dx.sqrnorm() < sqr(sphere_r)) {
        vfcu[l][c] = k;
      }
    }
  }
  tracer_->SetVolumeFraction(vfcu);
}

template <class M>
void Hydro<M>::OverwriteBc() {
  // piecewise-linear function
  auto piecewise = [&](Scal t, const std::vector<Scal>& times,
                       const std::vector<Scal>& values) {
    fassert_equal(values.size(), times.size());
    if (times.size() == 0) {
      return GetNan<Scal>();
    }
    size_t i = 0;
    while (i < times.size() && times[i] <= t) {
      ++i;
    }
    if (i == 0) { // t < times[0]
      return GetNan<Scal>();
    }
    if (i < times.size()) { // times[i - 1] <= t < times[i]
      const Scal t0 = times[i - 1];
      const Scal t1 = times[i];
      const Scal v0 = values[i - 1];
      const Scal v1 = values[i];
      return t0 < t1 ? v0 + (v1 - v0) * (t - t0) / (t1 - t0) : v0;
    } else {
      return values.back();
    }
  };
  const auto factor = piecewise(
      fs_->GetTime(), var.Vect["overwrite_inlet_times"],
      var.Vect["overwrite_inlet_factors"]);
  if (!IsNan(factor)) {
    mebc_fluid_.LoopPairs([&](auto cf_bc) {
      auto& curr = mebc_fluid_[cf_bc.first];
      const auto& orig = mebc_fluid_orig_[cf_bc.first];
      if (curr.type == BCondFluidType::inlet) {
        fassert(orig.type == BCondFluidType::inlet);
        curr.velocity = orig.velocity * factor;
      }
    });
  }

  const auto p = piecewise(
      fs_->GetTime(), var.Vect["overwrite_inletpressure_times"],
      var.Vect["overwrite_inletpressure_pressure"]);
  if (!IsNan(p)) {
    mebc_fluid_.LoopPairs([&](auto cf_bc) {
      auto& curr = mebc_fluid_[cf_bc.first];
      if (curr.type == BCondFluidType::inletpressure) {
        curr.pressure = p;
      }
    });
  }
}

template <class M>
void Hydro<M>::SpawnParticles(ParticlesView& view) {
  const Vect sphere_c(var.Vect["particles_spawn_sphere_c"]);
  const Scal sphere_r = var.Double["particles_spawn_sphere_r"];
  // particles per unit time
  const size_t dim = m.GetEdim();
  const Scal sphere_vol =
      (dim == 3 ? 4. / 3. * M_PI * std::pow(sphere_r, 3)
                : M_PI * sqr(sphere_r));
  const Vect h = m.GetCellSize();
  const Scal cell_vol = (dim == 3 ? h.prod() : h[0] * h[1]);
  const Vect velocity(var.Vect["particles_spawn_velocity"]);
  const Scal density = var.Double["particles_density"];
  const auto spawn_rate = var.Vect["particles_spawn_rate"];
  const auto diameter = var.Vect["particles_diameter"];
  const auto termvel = var.Vect["particles_termvel"];
  const size_t num_rates = spawn_rate.size();
  fassert_equal(diameter.size(), num_rates);
  fassert_equal(termvel.size(), num_rates);
  std::uniform_real_distribution<Scal> u(0, 1);
  std::uniform_real_distribution<Scal> um(-0.5, 0.5);
  auto& g = randgen_;

  for (auto c : m.Cells()) {
    const auto xc = m.GetCenter(c);
    Vect dx = xc - sphere_c;
    if (dim == 2) {
      dx[2] = 0;
    }
    if (dx.sqrnorm() < sqr(sphere_r)) {
      for (size_t i = 0; i < num_rates; ++i) {
        const Scal prob = particles_dt_ * cell_vol * spawn_rate[i] / sphere_vol;
        if (u(g) < prob) {
          view.x.push_back(m.GetCenter(c) + Vect(um(g), um(g), um(g)) * h);
          view.v.push_back(velocity);
          view.r.push_back(diameter[i] * 0.5);
          view.rho.push_back(density);
          view.termvel.push_back(termvel[i]);
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::InitParticles() {
  if (var.Int["enable_particles"]) {
    typename ParticlesInterface<M>::Conf conf;
    conf.mixture_density = var.Double["rho1"];
    conf.mixture_viscosity = var.Double["mu1"];
    conf.gravity = Vect(var.Vect["gravity"]);

    std::vector<Vect> p_x;
    std::vector<Vect> p_v;
    std::vector<Scal> p_r;
    std::vector<Scal> p_rho;
    std::vector<Scal> p_termvel;
    ParticlesView init{p_x, p_v, p_r, p_rho, p_termvel};
    SpawnParticles(init);
    if (eb_) {
      particles_.reset(new Particles<EB>(m, *eb_, init, fs_->GetTime(), conf));
    } else {
      particles_.reset(new Particles<M>(m, m, init, fs_->GetTime(), conf));
    }
  }
}

template <class M>
void Hydro<M>::InitTracer(Multi<FieldCell<Scal>>& vfcu) {
  if (var.Int["enable_tracer"]) {
    auto multi = [](const std::vector<Scal>& v) {
      Multi<Scal> w(v.size());
      w.data() = v;
      return w;
    };
    typename TracerInterface<M>::Conf conf;
    conf.layers = var.Int["tracer_layers"];
    const auto trl = GRange<size_t>(conf.layers);

    conf.density = multi(var.Vect["tracer_density"]);
    fassert(conf.density.size() >= conf.layers);

    conf.viscosity = multi(var.Vect["tracer_viscosity"]);
    fassert(conf.viscosity.size() >= conf.layers);
    conf.gravity = Vect(var.Vect["gravity"]);

    auto termvel = multi(var.Vect["tracer_termvel"]);
    conf.diameter.resize(conf.layers);
    if (var.Int["tracer_use_termvel"]) {
      const Scal mixture_density = var.Double["rho1"];
      for (auto l : trl) {
        conf.diameter[l] = std::sqrt(std::abs(
            18 * conf.viscosity[l] * termvel[l] /
            (conf.gravity.norm() * (conf.density[l] - mixture_density))));
      }
    } else {
      conf.diameter = multi(var.Vect["tracer_diameter"]);
    }
    fassert(conf.diameter.size() >= conf.layers);

    conf.scheme = GetConvSc(var.String["tracer_scheme"]);

    using SlipType = typename TracerInterface<M>::SlipType;
    conf.slip.resize(conf.layers);
    for (auto l : trl) {
      auto& slip = conf.slip[l];
      std::stringstream arg(var.String["tracer" + std::to_string(l) + "_slip"]);
      std::string type;
      arg >> type;
      if (type == "none") {
        slip.type = SlipType::none;
      } else if (type == "stokes") {
        slip.type = SlipType::stokes;
      } else if (type == "termvel") {
        slip.type = SlipType::constant;
        slip.velocity = conf.gravity * (termvel[l] / conf.gravity.norm());
      } else if (type == "constant") {
        slip.type = SlipType::constant;
        arg >> slip.velocity;
      } else {
        throw std::runtime_error(FILELINE + "Unknown slip='" + type + "'");
      }
    }
    Multi<MapEmbed<BCond<Scal>>> vmebc(conf.layers); // boundary conditions
    mebc_fluid_.LoopPairs([&](auto cf_bc) {
      for (auto l : trl) {
        if (l == 0) {
          vmebc[l][cf_bc.first] =
              BCond<Scal>(BCondType::dirichlet, cf_bc.second.nci, 1.);
        } else {
          vmebc[l][cf_bc.first] =
              BCond<Scal>(BCondType::dirichlet, cf_bc.second.nci, 0.);
        }
      }
    });

    if (eb_) {
      tracer_.reset(new Tracer<EB>(m, *eb_, vfcu, vmebc, fs_->GetTime(), conf));
    } else {
      tracer_.reset(new Tracer<M>(m, m, vfcu, vmebc, fs_->GetTime(), conf));
    }
  }
}

template <class M>
void Hydro<M>::InitFluid(const FieldCell<Vect>& fc_vel) {
  fcvm_ = fc_vel;

  std::string fs = var.String["fluid_solver"];
  if (fs == "simple") {
    auto p = ParsePar<Simple<M>>()(var);
    if (eb_) {
      fs_.reset(new Simple<Embed<M>>(
          m, *eb_, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_,
          &fc_force_, &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    } else {
      fs_.reset(new Simple<M>(
          m, m, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
          &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    }
  } else if (fs == "proj") {
    auto p = ParsePar<Proj<M>>()(var);
    if (eb_) {
      fs_.reset(new Proj<Embed<M>>(
          m, *eb_, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_,
          &fc_force_, &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    } else {
      fs_.reset(new Proj<M>(
          m, m, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
          &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    }
  } else if (fs == "dummy") {
    if (eb_) {
      fs_.reset(new FluidDummy<Embed<M>>(
          m, *eb_, fc_vel, &fc_rho_, &fc_mu_, &fc_force_, &febp_, &fc_src_,
          &fc_srcm_, 0., st_.dt, var));
    } else {
      fs_.reset(new FluidDummy<M>(
          m, m, fc_vel, &fc_rho_, &fc_mu_, &fc_force_, &febp_, &fc_src_,
          &fc_srcm_, 0., st_.dt, var));
    }
  } else {
    throw std::runtime_error("Unknown fluid_solver=" + fs);
  }
}

template <class M>
void Hydro<M>::InitAdvection(
    const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl) {
  std::string as = var.String["advection_solver"];
  if (as == "vof") {
    if (eb_) {
      auto p = ParsePar<ASVEB>()(var);
      as_.reset(new ASVEB(
          m, *eb_, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p));
    } else {
      auto p = ParsePar<ASV>()(var);
      as_.reset(new ASV(
          m, m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p));
    }
    layers = GRange<size_t>(1);
  } else if (as == "vofm") {
    if (eb_) {
      auto p = ParsePar<ASVMEB>()(var);
      auto as = new ASVMEB(
          m, *eb_, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    } else {
      auto p = ParsePar<ASVM>()(var);
      auto as = new ASVM(
          m, m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    }
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }

  fck_.resize(layers);
  fck_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
  auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize().norminf(), var);
  psm_par_ = ParsePar<PartStrMeshM<M>>()(ps, var);
}

template <class M>
void Hydro<M>::InitStat() {
  stat_.reset(new Stat<M>(m, eb_.get()));
  auto& stat = *stat_;
  auto& vf = as_->GetField();
  auto& vel = fs_->GetVelocity();
  auto& rho = fc_rho_;
  auto& p = fs_->GetPressure();
  stat.AddSum(
      "vol", "volume of domain", //
      [](IdxCell c, const M& m) { return m.GetVolume(c); });
  stat.AddSum(
      "vol1", "volume of phase 1", //
      [&vf](IdxCell c, const M& m) { return (1 - vf[c]) * m.GetVolume(c); });
  if (eb_) {
    stat.AddSum(
        "vol_eb", "volume of domain with embed", //
        [](IdxCell c, const EB& eb) { return eb.GetVolume(c); });
    stat.AddSum(
        "vol1_eb", "volume of phase 1 with embed", //
        [&vf](IdxCell c, const EB& eb) {
          return eb.GetVolume(c) - vf[c] * eb.GetMesh().GetVolume(c);
        });
  }
  stat.AddSum(
      "vol2", "volume of phase 2", //
      [&vf](IdxCell c, const M& m) { return vf[c] * m.GetVolume(c); });
  stat.AddSum(
      "ekin", "kinetic energy of mixture", //
      [&vel, &rho](IdxCell c, const M& m) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * m.GetVolume(c);
      });
  stat.AddSum(
      "ekin1", "kinetic energy of phase 1", //
      [&vf, &vel, &rho](IdxCell c, const M& m) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * m.GetVolume(c) * (1 - vf[c]);
      });
  stat.AddSum(
      "ekin2", "kinetic energy of phase 2", //
      [&vf, &vel, &rho](IdxCell c, const M& m) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * m.GetVolume(c) * vf[c];
      });
  stat.AddMax(
      "pmax", "maximum pressure", //
      [&p](IdxCell c, const M&) { return p[c]; });
  stat.AddMin(
      "pmin", "minimum pressure", //
      [&p](IdxCell c, const M&) { return p[c]; });
  stat.AddNone("t", "time", [&fs_ = fs_]() { return fs_->GetTime(); });
  stat.AddNone(
      "dt", "fluid time step", //
      [&fs_ = fs_]() { return fs_->GetTimeStep(); });
  stat.AddNone(
      "dta", "advection time step", //
      [&as_ = as_]() { return as_->GetTimeStep(); });
  stat.AddNone(
      "wt", "wall-clock time", //
      [&timer_ = timer_]() { return timer_.GetSeconds(); });
  stat.AddNone("iter", "iteration", [&st_ = st_]() { return st_.iter; });
  stat.AddNone("step", "step", [&st_ = st_]() { return st_.step; });
  stat.AddNone(
      "diff", "velocity difference between last two iterations", //
      [&fs_ = fs_]() { return fs_->GetError(); });

  stat.AddSumHidden("x*vol1", "", [&vf](IdxCell c, const M& m) { //
    return m.GetCenter(c) * m.GetVolume(c) * (1 - vf[c]);
  });
  stat.AddSumHidden("x*vol2", "", [&vf](IdxCell c, const M& m) { //
    return m.GetCenter(c) * m.GetVolume(c) * vf[c];
  });
  stat.AddSumHidden("v*vol1", "", [&vf, &vel](IdxCell c, const M& m) { //
    return vel[c] * m.GetVolume(c) * (1 - vf[c]);
  });
  stat.AddSumHidden("v*vol2", "", [&vf, &vel](IdxCell c, const M& m) { //
    return vel[c] * m.GetVolume(c) * vf[c];
  });
  stat.AddSumHidden("p*vol1", "", [&p, &vf](IdxCell c, const M& m) { //
    return p[c] * (1 - vf[c]) * m.GetVolume(c);
  });
  stat.AddSumHidden("p*vol2", "", [&p, &vf](IdxCell c, const M& m) { //
    return p[c] * vf[c] * m.GetVolume(c);
  });
  auto& ffv = fs_->GetVolumeFlux();
  if (eb_) {
    stat.AddSum(
        "q_inlet", "inlet volume rate", //
        [&ffv, this, &m=m]() {
          Scal sum = 0;
          mebc_fluid_.LoopBCond(*eb_, [&](auto cf, IdxCell c, auto bc) { //
            if (m.IsInner(c)) {
              if (bc.type == BCondFluidType::inlet ||
                  bc.type == BCondFluidType::inletflux) {
                sum += ffv[cf] * (bc.nci == 0 ? -1 : 1);
              }
            }
          });
          return sum;
        });
    stat.AddSum(
        "q_inletpressure", "inletpressure volume rate", //
        [&ffv, this, &m=m]() {
          Scal sum = 0;
          mebc_fluid_.LoopBCond(*eb_, [&](auto cf, IdxCell c, auto bc) { //
            if (m.IsInner(c)) {
              if (bc.type == BCondFluidType::inletpressure) {
                sum += ffv[cf] * (bc.nci == 0 ? -1 : 1);
              }
            }
          });
          return sum;
        });
    stat.AddSum(
        "q_outlet", "outlet volume rate", //
        [&ffv, this, &m=m]() {
          Scal sum = 0;
          mebc_fluid_.LoopBCond(*eb_, [&](auto cf, IdxCell c, auto bc) { //
            if (m.IsInner(c)) {
              if (bc.type == BCondFluidType::outlet ||
                  bc.type == BCondFluidType::outletpressure) {
                sum += ffv[cf] * (bc.nci == 0 ? 1 : -1);
              }
            }
          });
          return sum;
        });
    stat.AddSum(
        "area_inletpressure", "inletpressure area", //
        [this, &m=m]() {
          Scal sum = 0;
          auto& eb = *eb_;
          mebc_fluid_.LoopBCond(*eb_, [&](auto cf, IdxCell c, auto bc) { //
            if (m.IsInner(c)) {
              if (bc.type == BCondFluidType::inletpressure) {
                sum += eb.GetArea(cf);
              }
            }
          });
          return sum;
        });
    stat.AddSumHidden(
        "p*area_inletpressure", "inletpressure pressure * area", //
        [&p, this, &m=m]() {
          Scal sum = 0;
          auto& eb = *eb_;
          mebc_fluid_.LoopBCond(*eb_, [&](auto cf, IdxCell c, auto bc) { //
            if (m.IsInner(c)) {
              if (bc.type == BCondFluidType::inletpressure) {
                sum += p[c] * eb.GetArea(cf);
              }
            }
          });
          return sum;
        });
  }

  auto div = [](auto v, Scal d) {
    if (d == 0) {
      return v * 0;
    }
    return v / d;
  };

  stat.AddDerived(
      "c1", "centeroid of phase 1", //
      [div](const Stat<M>& stat) {
        return div(stat.vect["x*vol1"], stat["vol1"]);
      });
  stat.AddDerived(
      "c2", "centeroid of phase 2", //
      [div](const Stat<M>& stat) {
        return div(stat.vect["x*vol2"], stat["vol2"]);
      });
  stat.AddDerived(
      "v1", "velocity of phase 1", //
      [div](const Stat<M>& stat) {
        return div(stat.vect["v*vol1"], stat["vol1"]);
      });
  stat.AddDerived(
      "v2", "velocity of phase 2", //
      [div](const Stat<M>& stat) {
        return div(stat.vect["v*vol2"], stat["vol2"]);
      });
  stat.AddDerived(
      "pd", "pressure max-min", //
      [div](const Stat<M>& stat) { return stat["pmax"] - stat["pmin"]; });
  stat.AddDerived(
      "pavg1", "average pressure of phase 1", //
      [div](const Stat<M>& stat) { return div(stat["p*vol1"], stat["vol1"]); });
  stat.AddDerived(
      "pavg2", "average pressure of phase 2", //
      [div](const Stat<M>& stat) { return div(stat["p*vol2"], stat["vol2"]); });
  // XXX: relies on alphabetical order, "vol2_0" < "vol2_diff"
  //      to have "vol2_0" defined before "vol2_diff"
  stat.AddDerived(
      "vol2_0", "initial volume of phase 2", //
      [](const Stat<M>& stat) {
        return stat["vol2_0"] == 0 ? stat["vol2"] : stat["vol2_0"];
      });
  stat.AddDerived(
      "vol2_diff", "relative difference between vol2 and vol2_0", //
      [div](const Stat<M>& stat) {
        return div(stat["vol2"] - stat["vol2_0"], stat["vol2_0"]);
      });
  stat.AddDerived(
      "p_inletpressure", "inletpressure average pressure", //
      [div](const Stat<M>& stat) {
        return div(
            stat["p*area_inletpressure"], stat["area_inletpressure"]);
      });

  if (var.Int["stat_dissip"]) {
    stat.AddSum(
        "dissip", "dissipation rate of mixture", //
        [&str = fc_strain_, &mu = fc_mu_](IdxCell c, const M& m) {
          return 2 * mu[c] * str[c] * m.GetVolume(c);
        });
    stat.AddSum(
        "dissip1", "dissipation rate of phase 1", //
        [&str = fc_strain_, &mu = fc_mu_, &vf](IdxCell c, const M& m) {
          return 2 * mu[c] * str[c] * (1 - vf[c]) * m.GetVolume(c);
        });
    stat.AddSum(
        "dissip2", "dissipation rate of phase 2", //
        [&str = fc_strain_, &mu = fc_mu_, &vf](IdxCell c, const M& m) {
          return 2 * mu[c] * str[c] * vf[c] * m.GetVolume(c);
        });
    stat.AddDerived(
        "edis", "dissipated energy of mixture", //
        [&fs_ = fs_](const Stat<M>& stat) {
          return stat["edis"] + fs_->GetTimeStep() * stat["dissip"];
        });
    stat.AddDerived(
        "edis1", "dissipated energy of phase 1", //
        [&fs_ = fs_](const Stat<M>& stat) {
          return stat["edis1"] + fs_->GetTimeStep() * stat["dissip1"];
        });
    stat.AddDerived(
        "edis2", "dissipated energy of phase 2", //
        [&fs_ = fs_](const Stat<M>& stat) {
          return stat["edis2"] + fs_->GetTimeStep() * stat["dissip2"];
        });
  }
  if (var.Int["enstrophy"]) {
    stat.AddSum(
        "enstr", "enstrophy", //
        [&omm = fcomm_, &rho = fc_rho_](IdxCell c, const M& m) {
          return 0.5 * sqr(omm[c]) * rho[c] * m.GetVolume(c);
        });
  }
  if (var.Int["statvel"]) {
    const Vect vel0(var.Vect["vel"]);
    stat.AddSumHidden(
        "dv*dv*vol", "", [&vel0, &vel](IdxCell c, const M& m) {
          auto dv = vel[c] - vel0;
          return dv * dv * m.GetVolume(c);
        });
    stat.AddSumHidden(
        "dv*dv*vol2", "", [&vel0, &vel, &vf](IdxCell c, const M& m) {
          auto dv = vel[c] - vel0;
          return dv * dv * vf[c] * m.GetVolume(c);
        });
    stat.AddMax(
        "dvel_max", "maximum velocity difference with `vel`",
        [&vel0, &vel](IdxCell c, const M&) {
          auto dv = vel[c] - vel0;
          return dv.abs();
        });
    stat.AddDerived(
        "dvel_l2", "L2 norm of velocity difference with `vel`", //
        [](const Stat<M>& stat) {
          const Vect v = stat.vect["dv*dv*vol"] / stat["vol"];
          return Vect(std::sqrt(v[0]), std::sqrt(v[1]), std::sqrt(v[2]));
        });
    stat.AddMax(
        "dvel2_max", "maximum phase 2 velocity difference with `vel`",
        [&vel0, &vel, &vf](IdxCell c, const M&) {
          auto dv = (vel[c] - vel0) * vf[c];
          return dv.abs();
        });
    stat.AddDerived(
        "dvel2_l2", "L2 norm of phase 2 velocity difference with `vel`", //
        [](const Stat<M>& stat) {
          const Vect dv = stat.vect["dv*dv*vol2"] / stat["vol2"];
          return Vect(std::sqrt(dv[0]), std::sqrt(dv[1]), std::sqrt(dv[2]));
        });
  }
  if (var.Int["stat_vofm"]) {
    auto add_vofm = [this, &stat](auto as) {
      for (auto l : layers) {
        auto sl = std::to_string(l);
        const auto& vf = *as->GetFieldM()[l];
        stat.AddSum(
            "vofm_cells_vf" + sl, "cells with positive volume fraction", //
            [&vf](IdxCell c, const M&) -> Scal {
              return vf[c] > 0 ? 1 : 0;
            });
        const auto& cl  = *as->GetColor()[l];
        stat.AddSum(
            "vofm_cells_cl" + sl, "cells with defined color", //
            [&cl](IdxCell c, const M&) -> Scal {
              return cl[c] > 0 ? 1 : 0;
            });
        stat.AddSum(
            "vofm_vol" + sl, "integral of volume fraction", //
            [&vf](IdxCell c, const M& m) {
              return vf[c] * m.GetVolume(c);
            });
        auto slp = std::to_string(l + 1);
        // TODO: revise with `++hist[cnt]`, now takes L^2 operations.
        stat.AddSum(
            "vofm_hist" + slp,
            "number of cells with " + slp + " non-empty layers", //
            [vfm = as->GetFieldM(), &layers = layers, cnt_target = l + 1](
                IdxCell c, const M&) -> Scal {
              size_t cnt = 0;
              for (auto l : layers) {
                if ((*vfm[l])[c] > 0) {
                  ++cnt;
                }
              }
              return cnt == cnt_target ? 1 : 0;
            });
      }
    };
    if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      add_vofm(as);
    }
    if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      add_vofm(as);
    }
  }
  if (var.Int["stat_area"]) {
    if (auto as = dynamic_cast<ASV*>(as_.get())) {
      auto& n = as->GetNormal();
      auto& a = as->GetAlpha();
      auto& vf = as->GetField();
      const Vect h = m.GetCellSize();
      stat.AddSum(
          "area", "area of the interface", //
          [&n, &a, &vf, &h](IdxCell c, const M&) -> Scal {
            if (vf[c] > 0 && vf[c] < 1 && !IsNan(a[c])) {
              using R = Reconst<Scal>;
              auto s = std::abs(R::GetArea(R::GetCutPoly2(n[c], a[c], h), n[c]));
              if (!IsNan(s)) {
                return s;
              }
            }
            return 0;
          });
    }
  }
  if (eb_) {
    stat.AddSum("pdrag", "pressure drag", [this]() {
      auto d = CalcPressureDrag(fs_->GetPressure(), *eb_);
      if (m.GetEdim() == 2) {
        d /= m.GetCellSize()[2];
      }
      return d;
    });
    stat.AddSum("vdrag", "viscous drag", [this]() {
      auto d = CalcViscousDrag(fs_->GetVelocity(), fc_mu_, *eb_);
      if (m.GetEdim() == 2) {
        d /= m.GetCellSize()[2];
      }
      return d;
    });
    stat.AddDerived("drag", "total drag", [](const Stat<M>& stat) {
      return stat.vect["pdrag"] + stat.vect["vdrag"];
    });
  }
  if (particles_) {
    stat.AddSum(
        "particles_n", "number of particles",
        [&p = particles_]() -> Scal { return p->GetView().x.size(); });
    stat.AddSum(
        "particles_nrecv", "number of particles transferred at last step",
        [&p = particles_]() -> Scal { return p->GetNumRecv(); });
  }

  stat_->SortNames();

  auto blacklist = GetWords(var.String("stat_blacklist", ""));
  for (auto name : stat_->GetNames()) {
    stat_->SetEnabled(name, !blacklist.count(name));
  }

  if (m.IsRoot()) {
    fstat_.open("stat.dat");
    fstat_.precision(16);
    stat_->WriteHeader(fstat_);

    std::ofstream fsum("stat_summary");
    stat_->WriteSummary(fsum, true);
  }
}

template <class M>
void Hydro<M>::Init() {
  using namespace fluid_condition;
  auto sem = m.GetSem("init");
  struct {
    FieldCell<Vect> fcvel; // initial velocity
    FieldCell<Scal> fcvf; // initial volume fraction
    FieldCell<Scal> fccl; // initial color
    FieldCell<Scal> fcbody;
    FieldCell<bool> fcbodymask;
    Multi<FieldCell<Scal>> tracer_vfcu;
    Vars varbody;
  } * ctx(sem);
  auto& fcvel = ctx->fcvel;
  auto& fcvf = ctx->fcvf;
  auto& fccl = ctx->fccl;
  if (sem("flags")) {
    m.flags.linreport = var.Int["linreport"];
    m.flags.check_symmetry = var.Int["check_symmetry"];
    m.flags.check_symmetry_dump_threshold =
        var.Double["check_symmetry_dump_threshold"];
    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.is_periodic[2] = var.Int["hypre_periodic_z"];
    randgen_.seed(m.GetId() + 1);
  }
  if (sem.Nested("embed")) {
    InitEmbed();
  }
  if (sem.Nested()) {
    InitVf(fcvf, var, m);
  }
  if (sem.Nested()) {
    if (var.Int["enable_tracer"]) {
      InitTracerFields(ctx->tracer_vfcu);
    }
  }
  if (sem("fields")) {
    if (eb_) {
      auto& eb = *eb_;
      for (auto c : m.AllCells()) {
        if (eb.GetType(c) == M::Type::excluded) {
          fcvf[c] = 0;
        }
      }
    }

    fc_src_.Reinit(m, 0.);
    fc_src2_.Reinit(m, 0.);
    fc_srcm_.Reinit(m, 0.);

    // initial surface tension sigma
    fc_sig_.Reinit(m, 0);
    auto isig = CreateInitSig<M>(var);
    isig(fc_sig_, m);
    m.Comm(&fc_sig_);

    // initial contact angle
    {
      fc_contang_.Reinit(m, -1);
      const std::string name = var.String["init_contang"];
      if (auto ptr = ModuleInitContang<M>::GetInstance(name)) {
        (*ptr)(fc_contang_, var, m);
      } else {
        if (m.IsRoot()) {
          std::cerr << "Known values of 'init_contang': ";
          for (auto& p : ModuleInitContang<M>::GetInstances()) {
            std::cerr << p.first << " ";
          }
          std::cerr << std::endl;
        }
        throw std::runtime_error(FILELINE + ": Unknown init_contang=" + name);
      }
      m.Comm(&fc_contang_);
    }

    // initial velocity
    fcvel.Reinit(m, Vect(0));
    InitVel(fcvel, var, m);
    if (eb_) {
      InitVelHook(fcvel, var, m, *eb_);
    } else {
      InitVelHook(fcvel, var, m);
    }
    m.Comm(&fcvel);
    fcvel.SetHalo(2);
    fcvel.SetName("fcvel");

    // global mesh size
    MIdx gs = m.GetGlobalSize();

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
      std::cout << "viscosity dt=" << GetVisDt() << std::endl;
    }

    // boundary conditions
    if (eb_) {
      auto p = InitBc(var, *eb_);
      mebc_fluid_ = std::get<0>(p);
      mf_adv_ = std::get<1>(p);
      me_group_ = std::get<2>(p);
      vdesc_ = std::get<3>(p);
      mf_adv_.LoopBCond(*eb_, [&](auto cf, auto c, auto) { //
        me_contang_[cf] = fc_contang_[c];
        mf_adv_[cf].contang = fc_contang_[c];
      });
    } else {
      auto p = InitBc(var, m);
      mebc_fluid_ = std::get<0>(p);
      mf_adv_ = std::get<1>(p);
      me_group_ = std::get<2>(p);
      vdesc_ = std::get<3>(p);
      for (auto& p : mf_adv_.GetMapFace()) {
        const auto& f = p.first;
        auto& bc = p.second;
        const auto c = m.GetCell(f, bc.nci);
        me_contang_[f] = fc_contang_[c];
        bc.contang = fc_contang_[c];
      }
    }
    mebc_fluid_orig_ = mebc_fluid_;

    // boundary conditions for smoothing of volume fraction
    for (auto& p : mebc_fluid_.GetMapFace()) {
      const IdxFace f = p.first;
      auto& bc = p.second;
      const auto nci = bc.nci;
      auto& bcvf = mebc_vfsm_[f];
      bcvf.nci = nci;
      bcvf.type = BCondType::neumann;
    }
  }

  if (var.Int["bc_wall_init_vel"] && sem("bc_wall_init_vel")) {
    // velocity on walls from initial conditions in neighbor cells
    for (auto& p : mebc_fluid_.GetMapFace()) {
      const IdxFace f = p.first;
      auto& bc = p.second;
      if (bc.type == BCondFluidType::wall) {
        const IdxCell c = m.GetCell(f, bc.nci);
        bc.velocity = fcvel[c];
      }
    }
  }

  if (var.Int["initvort"] && sem.Nested("initvort")) {
    InitVort(fcvel, fcvel, mebc_fluid_, m);
  }

  if (var.Int["vel_init_random"] && sem("random")) {
    Scal amp = var.Double["random_amp"];
    Vect vel(var.Vect["random_vel"]);
    std::default_random_engine g(m.GetId());
    std::uniform_real_distribution<Scal> u(-amp, amp);
    for (auto c : m.Cells()) {
      Vect v = vel * u(g) / 7;
      for (auto q : m.Nci(c)) {
        IdxCell cn = m.GetCell(c, q);
        fcvel[cn] += v;
      }
      fcvel[c] += v;
    }
    m.Comm(&fcvel);
  }

  if (sem.Nested("smooth")) {
    Smoothen(fcvf, mebc_vfsm_, m, var.Int["vf_init_sm"]);
  }

  if (sem.Nested("mixture")) {
    CalcMixture(fcvf);
  }

  if (sem("color-ini")) {
    if (var.Int["enable_color"]) {
      // initial color
      // TODO revise with bcast
      auto icl = CreateInitCl<M>(var, m.IsRoot());
      icl(fccl, fcvf, m);
      m.Comm(&fccl);
    } else {
      fccl.Reinit(m, 0.);
    }
  }

  if (sem.Nested("cellcond")) {
    GetFluidCellCond(var, m, mc_velcond_);
  }

  if (sem("body-mask")) {
    auto& varbody = ctx->varbody;
    varbody.String.Set("init_vf", var.String["body_init"]);
    varbody.String.Set("list_path", var.String["body_list_path"]);
    varbody.Int.Set("dim", var.Int["dim"]);
    varbody.Int.Set("list_ls", 3);
  }
  if (sem.Nested("body-mask")) {
    InitVf(ctx->fcbody, ctx->varbody, m);
  }
  /*
  // TODO: implement
  if (sem("body-bc")) {
    // Step-wise approximation of bodies
    const Scal clear0 = var.Double["bcc_clear0"];
    const Scal clear1 = var.Double["bcc_clear1"];
    const Scal inletcl = var.Double["inletcl"];
    const Scal fill_vf = var.Double["bcc_fill"];
    auto& fc = ctx->fcbodymask;
    fc.Reinit(m, false);
    for (auto c : m.AllCells()) {
      fc[c] = (ctx->fcbody[c] > 0.5);
    }
    if (var.Int["body_init_inverse"]) {
      for (auto c : m.AllCells()) {
        fc[c] = !fc[c];
      }
    }
    AppendBodyCond<M>(
        fc, var.String["body_bc"], m, clear0, clear1, inletcl, fill_vf, nullptr,
        mebc_fluid_, mf_adv_);
  }
  */

  if (sem("calcdt0")) {
    const Scal dt = var.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    tracer_dt_ = dt;
    particles_dt_ = dt;
  }
  if (sem("solv")) {
    InitFluid(fcvel);

    InitAdvection(fcvf, fccl);

    InitTracer(ctx->tracer_vfcu);

    InitParticles();

    st_.iter = 0;
    st_.step = 0;

    if (m.IsLead()) {
      this->var_mutable.Int.Set("iter", st_.iter);
    }

    InitStat();

    if (var.Int["fill_halo_nan"]) {
      std::vector<std::pair<IdxFace, size_t>> vf;
      for (auto& p : mf_adv_.GetMapFace()) {
        vf.emplace_back(p.first, p.second.GetNci());
      }
      m.SetNanFaces(vf);
      m.flags.nan_faces_value = var.Double["fill_halo_nan_value"];
    }

    events_ = std::unique_ptr<Events>(
        new Events(this->var_mutable, m.IsRoot(), m.IsLead()));
    events_->AddHandler(
        "vf_save_state",
        [&path = vf_save_state_path_](std::string arg) { //
          path = arg;
        });
    events_->Parse();
  }
  if (sem.Nested("vofm-load")) {
    if (var.Int["init_vf_load_state"]) {
      if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
        as->LoadState(var.String["init_vf_state_dir"]);
      }
    }
  }
  if (var.Int["dumpbc"]) {
    if (vdesc_.size()) {
      if (sem("dump-bcgroups")) {
        if (m.IsRoot()) {
          std::ofstream fdesc("bc_groups.dat");
          for (size_t i = 0; i < vdesc_.size(); ++i) {
            fdesc << i << " " << vdesc_[i] << std::endl;
          }
        }
      }
      if (sem.Nested("bcdump")) {
        if (eb_) {
          UInitEmbedBc<M>::DumpPoly("bc.vtk", me_group_, me_contang_, *eb_, m);
        } else {
          UInitEmbedBc<M>::DumpPoly("bc.vtk", me_group_, me_contang_, m, m);
        }
      }
    }
  }
  if (eb_ && sem.Nested()) {
    eb_->DumpPoly(var.Int["vtkbin"], var.Int["vtkmerge"]);
  }
  if (var.Int["dumpinit"]) {
    if (sem.Nested()) {
      Dump(true);
    }
  }
}

template <class M>
Hydro<M>::Hydro(Vars& var0, const MyBlockInfo& bi, Par& par)
    : KernelMeshPar<M, Par>(var0, bi, par)
    , dumper_(var, "dump_field_")
    , dmptraj_(var, "dump_traj_")
    , dmptrep_(var, "dump_trep_") {}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");
  if (sem("local")) {
    if (var.Int["stat_dissip"]) {
      fc_strain_ = CalcStrain(fs_->GetVelocity());
    }
    if (var.Int["enstrophy"]) {
      CalcVort();
    }
  }
}

template <class M>
void Hydro<M>::CalcDt() {
  auto sem = m.GetSem("calcdt");
  struct {
    Scal dtmin;
  } * ctx(sem);

  if (sem("local")) {
    st_.t = fs_->GetTime();
    ctx->dtmin = fs_->GetAutoTimeStep();
    m.Reduce(&ctx->dtmin, "min");
  }
  if (sem("reduce")) {
    // set from cfl if defined
    if (auto* cfl = var.Double.Find("cfl")) {
      st_.dt = ctx->dtmin * (*cfl);
      st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
    }

    // constraint from surface tension
    st_.dt = std::min<Scal>(st_.dt, GetStDt());

    // constraint from viscosity
    st_.dt = std::min<Scal>(st_.dt, GetVisDt());

    fs_->SetTimeStep(st_.dt);

    // set from cfla if defined
    if (auto* cfla = var.Double.Find("cfla")) {
      st_.dta = ctx->dtmin * (*cfla);
      st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
    }
    // round up dta to such that dt / dta is integer
    const Scal dt = fs_->GetTime() + fs_->GetTimeStep() - as_->GetTime();
    st_.dta = dt / std::max(1, int(dt / st_.dta + 0.5));
    as_->SetTimeStep(st_.dta);

    if (tracer_) {
      if (auto* cflt = var.Double.Find("cflt")) {
        tracer_dt_ = ctx->dtmin * (*cflt);
        // round up dta to such that dt / dta is integer
        const Scal dtwhole =
            fs_->GetTime() + fs_->GetTimeStep() - tracer_->GetTime();
        tracer_dt_ = dtwhole / std::max(1, int(dtwhole / tracer_dt_ + 0.5));
      } else {
        tracer_dt_ = fs_->GetTimeStep();
      }
    }

    if (particles_) {
      if (auto* cflp = var.Double.Find("cflp")) {
        particles_dt_ = ctx->dtmin * (*cflp);
        // round up dta to such that dt / dta is integer
        const Scal dtwhole =
            fs_->GetTime() + fs_->GetTimeStep() - particles_->GetTime();
        particles_dt_ =
            dtwhole / std::max(1, int(dtwhole / particles_dt_ + 0.5));
      } else {
        particles_dt_ = fs_->GetTimeStep();
      }
    }
  }
  if (sem()) {
    // FIXME: empty stage
  }
}

template <class M>
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf0) {
  auto sem = m.GetSem("mixture");

  if (sem("init")) {
    fc_mu_.Reinit(m);
    fc_rho_.Reinit(m);
    fc_force_.Reinit(m, Vect(0));
    febp_.Reinit(m, 0);
    fc_smvf_ = fc_vf0;
  }

  if (sem.Nested("smooth")) {
    Smoothen(fc_smvf_, mebc_vfsm_, m, var.Int["vfsmooth"]);
  }

  if (sem("calc")) {
    FieldCell<Scal>& a = fc_smvf_;
    FieldFace<Scal>& af = ff_smvf_;
    if (eb_) {
      auto& eb = *eb_;
      af = UEB::Interpolate(a, {}, eb).GetFieldFace();
    } else {
      af = UEB::Interpolate(a, mebc_vfsm_, m);
    }

    const Vect force(var.Vect["force"]);
    const Vect grav(var.Vect["gravity"]);
    const Scal rho1(var.Double["rho1"]);
    const Scal rho2(var.Double["rho2"]);
    const Scal mu1(var.Double["mu1"]);
    const Scal mu2(var.Double["mu2"]);

    // Init density and viscosity
    for (auto c : m.AllCells()) {
      const Scal a2 = a[c];
      const Scal a1 = 1 - a2;
      fc_rho_[c] = rho1 * a1 + rho2 * a2;
      fc_mu_[c] = mu1 * a1 + mu2 * a2;
    }
    FieldFace<Scal> ff_rho(m);
    for (auto f : m.AllFaces()) {
      const Scal a2 = af[f];
      const Scal a1 = 1 - a2;
      ff_rho[f] = rho1 * a1 + rho2 * a2;
    }

    if (tracer_ && var.Int["tracer_override_mixture"]) {
      const auto& fc_rho_mix = tracer_->GetMixtureDensity();
      const auto& fc_mu_mix = tracer_->GetMixtureViscosity();
      for (auto c : m.AllCells()) {
        const Scal a2 = a[c];
        const Scal a1 = 1. - a2;
        fc_rho_[c] = fc_rho_mix[c] * a1 + rho2 * a2;
        fc_mu_[c] = fc_mu_mix[c] * a1 + mu2 * a2;
      }
      FieldFace<Scal> ff_rho_mix(m);
      if (eb_) {
        ff_rho_mix = UEB::Interpolate(fc_rho_mix, {}, *eb_).GetFieldFace();
      } else {
        ff_rho_mix = UEB::Interpolate(fc_rho_mix, mebc_vfsm_, m);
      }
      for (auto f : m.AllFaces()) {
        const Scal a2 = af[f];
        const Scal a1 = 1. - a2;
        ff_rho[f] = ff_rho_mix[f] * a1 + rho2 * a2;
      }
    }

    if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      FieldEmbed<Scal> ffvisc(m, 0);
      auto& eb = *eb_;
      auto fccl = as->GetColor();
      auto fcu = as->GetFieldM();
      const Scal musurf(var.Double["musurf"]);
      if (musurf) {
        for (auto f : eb.SuFaces()) {
          const IdxCell cm = eb.GetCell(f, 0);
          const IdxCell cp = eb.GetCell(f, 1);
          std::set<Scal> s;
          for (auto i : layers) {
            const Scal clm = (*fccl[i])[cm];
            const Scal clp = (*fccl[i])[cp];
            if (clm != kClNone) s.insert(clm);
            if (clp != kClNone) s.insert(clp);
          }
          for (auto cl : s) {
            Scal um = 0;
            Scal up = 0;
            for (auto i : layers) {
              if ((*fccl[i])[cm] == cl) {
                um = (*fcu[i])[cm];
              }
              if ((*fccl[i])[cp] == cl) {
                up = (*fcu[i])[cp];
              }
            }
            ffvisc[f] = std::abs(up - um) * musurf;
          }
        }
      }
      auto fcadd = UEB::Interpolate(ffvisc, eb);
      for (auto c : eb.Cells()) {
        fc_mu_[c] += fcadd[c];
      }
    }


    // Append gravity to force
    for (auto f : m.AllFaces()) {
      const Vect n = m.GetNormal(f);
      febp_[f] += force.dot(n);
      febp_[f] += grav.dot(n) * ff_rho[f];
    }

    // Surface tension
    if (var.Int["enable_surftens"] && as_) {
      CalcSurfaceTension(
          m, layers, var, fc_force_, febp_.GetFieldFace(), fc_sig_,
          GetBCondZeroGrad<Scal>(mebc_fluid_), fck_, fc_vf0, af, as_.get());
    }

    // zero force in z if 2D
    if (var.Int["dim"] <= 2) {
      for (auto f : m.Faces()) {
        using Dir = typename M::Dir;
        if (m.GetIndexFaces().GetDir(f) == Dir::k) {
          febp_[f] = 0.; // XXX: zero in z
        }
      }
    }
  }
  // FIXME move, but keep inside nested call
  if (!vf_save_state_path_.empty() && sem.Nested("vf_save_state")) {
    if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      as->SaveState(vf_save_state_path_);
    }
  }
}

template <class M>
void Hydro<M>::DumpFields() {
  auto sem = m.GetSem("dumpfields");
  struct {
    std::array<Multi<FieldCell<Scal>>, dim> im;
    FieldCell<Scal> fc_cellcond;
    FieldCell<Scal> fcdiv; // divergence of velocity
    FieldCell<Scal> fcdis; // energy dissipation
    FieldCell<Scal> fc_ebvf; // embedded boundaries volume fraction
    FieldCell<Scal> fc_tracer_sum; // sum of tracer_ fields starting from 1
  } * ctx(sem);
  if (sem("dump")) {
    if (m.IsRoot()) {
      dumper_.Report(std::cout);
    }

    auto dl = GetWords(var.String["dumplist"]);

    auto dump = [&dl, this](const FieldCell<Scal>& fc, std::string name) {
      if (dl.count(name)) {
        m.Dump(&fc, name);
      }
    };
    auto dumpv = [&dl, this](
                     const FieldCell<Vect>& fc, size_t i, std::string name) {
      if (dl.count(name)) {
        m.Dump(&fc, i, name);
      }
    };

    auto& fcv = fs_->GetVelocity();
    dumpv(fcv, 0, "vx");
    dumpv(fcv, 1, "vy");
    dumpv(fcv, 2, "vz");
    dump(fs_->GetPressure(), "p");
    dump(as_->GetField(), "vf");
    dump(fc_rho_, "rho");
    dump(fc_mu_, "mu");
    dump(fc_sig_, "sig");
    dump(fc_contang_, "contang");
    if (dl.count("cellcond")) {
      auto& fc = ctx->fc_cellcond;
      fc.Reinit(m, 0);
      for (auto& it : mc_velcond_) {
        fc[it.first] = 1;
      }
      m.Dump(&fc, "cellcond");
    }
    if (dl.count("omx") || dl.count("omy") || dl.count("omz") ||
        dl.count("omm") || dl.count("omcalc")) {
      CalcVort();
      dumpv(fcom_, 0, "omx");
      dumpv(fcom_, 1, "omy");
      dumpv(fcom_, 2, "omz");
      dump(fcomm_, "omm");
    }
    if (dl.count("dis") || dl.count("strain")) {
      fc_strain_ = CalcStrain(fs_->GetVelocity());
      if (dl.count("strain")) m.Dump(&fc_strain_, "strain");
      if (dl.count("dis")) {
        ctx->fcdis = fc_strain_;
        for (auto c : m.Cells()) {
          ctx->fcdis[c] *= 2. * fc_mu_[c];
        }
        m.Dump(&ctx->fcdis, "dis");
      }
    }
    if (dl.count("div")) {
      ctx->fcdiv = GetDiv();
      m.Dump(&ctx->fcdiv, "div");
    }
    if (auto as = dynamic_cast<ASV*>(as_.get())) {
      dumpv(as->GetNormal(), 0, "nx");
      dumpv(as->GetNormal(), 1, "ny");
      dumpv(as->GetNormal(), 2, "nz");
      dump(as->GetColor(), "cls");
      dump(fck_[0], "k");
    }
    // TODO reuse ASV code
    if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
      dumpv(as->GetNormal(), 0, "nx");
      dumpv(as->GetNormal(), 1, "ny");
      dumpv(as->GetNormal(), 2, "nz");
      dump(as->GetColor(), "cls");
      dump(fck_[0], "k");
    }
    if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      for (auto l : layers) {
        auto sl = std::to_string(l);
        dump(*as->GetFieldM()[l], "vf" + sl);
        dump(*as->GetColor()[l], "cl" + sl);
        dumpv(*as->GetNormal()[l], 0, "nx" + sl);
        dumpv(*as->GetNormal()[l], 1, "ny" + sl);
        dumpv(*as->GetNormal()[l], 2, "nz" + sl);
        dump(fck_[l], "k" + sl);
      }

      // combined colors
      dump(as->GetColorSum(), "cls");

      // image
      auto conv = [&](size_t d, size_t l,
                      Multi<FieldCell<Scal>>& fc) -> const FieldCell<Scal>& {
        fc.resize(layers);
        fc[l].Reinit(m);
        for (auto c : m.Cells()) {
          fc[l][c] = as->GetImage(l, c)[d];
        }
        return fc[l];
      };
      for (auto d : {0, 1, 2}) {
        for (auto l : layers) {
          std::stringstream st;
          st << "im"
             << "xyz"[d] << l;
          std::string s = st.str();
          dump(conv(d, l, ctx->im[d]), s);
        }
      }
    }
    // TODO add ASVMEB

    if (eb_) {
      auto& eb = *eb_;
      if (dl.count("ebvf")) {
        auto& fc = ctx->fc_ebvf;
        fc.Reinit(m, 0);
        for (auto c : eb.Cells()) {
          fc[c] = eb.GetVolumeFraction(c);
        }
        m.Dump(&fc, "ebvf");
      }
      if (fc_dist_.size()) {
        dump(fc_dist_, "ebdist");
      }
      if (fc_phi_.size()) {
        dump(fc_phi_, "ebphi");
      }
    }

    if (tracer_) {
      for (auto l : tracer_->GetView().layers) {
        dump(tracer_->GetVolumeFraction()[l], "tu" + std::to_string(l));
      }
      if (dl.count("tusum")) {
        ctx->fc_tracer_sum.Reinit(m, 0);
        for (auto l : tracer_->GetView().layers) {
          if (l > 0) {
            const auto& fc = tracer_->GetVolumeFraction()[l];
            for (auto c : m.Cells()) {
              ctx->fc_tracer_sum[c] += fc[c];
            }
          }
        }
        dump(ctx->fc_tracer_sum, "tusum");
      }
    }
  }
  if (var.Int["enable_advection"]) {
    if (var.Int["dumppoly"] && sem.Nested()) {
      as_->DumpInterface(GetDumpName("s", ".vtk", dumper_.GetN()));
    }
    if (var.Int["dumppolymarch"] && sem.Nested()) {
      as_->DumpInterfaceMarch(GetDumpName("sm", ".vtk", dumper_.GetN()));
    }
  }
  if (particles_ && var.Int["dump_particles"]) {
    const std::string path = GetDumpName("part", ".csv", dumper_.GetN());
    if (sem()) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << particles_->GetTime() << " to " << path
                  << std::endl;
      }
    }
    if (sem.Nested()) {
      particles_->DumpCsv(path);
    }
  }
  if (sem()) {
    // XXX: empty stage, otherwise ctx is destroyed before dump
  }
}

template <class M>
void Hydro<M>::Dump(bool force) {
  auto sem = m.GetSem("dump");
  struct {
    Multi<FieldCell<MIdx>> fcim;
  } * ctx(sem);
  if (sem.Nested("fields")) {
    if (dumper_.Try(st_.t, st_.dt) || force) {
      DumpFields();
    }
  }
  if (dmptraj_.Try(st_.t, st_.dt) || force) {
    if (sem("copyimage")) {
      ctx->fcim.resize(layers);
      ctx->fcim.InitAll(FieldCell<MIdx>(m));
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        for (auto c : m.AllCells()) {
          for (auto l : layers) {
            ctx->fcim[l][c] = as->GetImage(l, c);
          }
        }
      } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
        for (auto c : m.AllCells()) {
          for (auto l : layers) {
            ctx->fcim[l][c] = as->GetImage(l, c);
          }
        }
      } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
        for (auto c : m.AllCells()) {
          ctx->fcim[0][c] = as->GetImage(c);
        }
      } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
        for (auto c : m.AllCells()) {
          ctx->fcim[0][c] = as->GetImage(c);
        }
      }
    }
    if (sem.Nested("trajdump")) {
      if (var.Int["enable_color"]) {
        Multi<const FieldCell<Scal>*> fcu(layers);
        Multi<const FieldCell<Scal>*> fccl(layers);
        if (auto as = dynamic_cast<ASVM*>(as_.get())) {
          fcu = as->GetFieldM();
          fccl = as->GetColor();
        } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
          // TODO reuse ASVM code
          fcu = as->GetFieldM();
          fccl = as->GetColor();
        } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
          // TODO reuse ASV code
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        }
        if (eb_) {
          DumpTraj<EB>(
              *eb_, true, var, dmptraj_.GetN(), st_.t, layers, fcu, fccl,
              ctx->fcim, fs_->GetPressure(), fs_->GetVelocity(), fcvm_, st_.dt);
        } else {
          DumpTraj<M>(
              m, true, var, dmptraj_.GetN(), st_.t, layers, fcu, fccl,
              ctx->fcim, fs_->GetPressure(), fs_->GetVelocity(), fcvm_, st_.dt);
        }
      }
    }
  }
  if (sem("dmptrep")) {
    if (m.IsRoot() && dmptrep_.Try(st_.t, st_.dt)) {
      std::string s = GetDumpName("trep", ".log", dmptrep_.GetN());
      m.TimerReport(s);
      std::cout << std::fixed << std::setprecision(8) << "timer report"
                << " t=" << st_.t << " to " << s << std::endl;
    }
  }
  if (sem("dumpstat")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("stat_step_every", 1) == 0 || force) {
        stat_->WriteValues(fstat_);
      }
    }
  }
  if (auto as = dynamic_cast<ASV*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  // TODO reuse ASV code
  if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  if (auto as = dynamic_cast<ASVM*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  // TODO reuse ASVM code
  if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
}

template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");
  struct {
    Scal nabort;
  } * ctx(sem);

  if (sem.Nested("init")) {
    Init();
  }

  sem.LoopBegin();

  if (sem("events")) {
    if (events_) {
      vf_save_state_path_ = "";
      events_->Execute(st_.t);
    }
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"]) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, t > tmax=" << var.Double["tmax"]
                  << std::endl;
      }
      sem.LoopBreak();
    } else if (int(st_.step + 0.5) >= var.Int["max_step"]) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, step > max_step="
                  << var.Int["max_step"] << std::endl;
      }
      sem.LoopBreak();
    } else if (st_.step > 1 && fs_->GetError() < var.Double("stop_diff", 0)) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, diff < stop_diff="
                  << var.Double["stop_diff"] << std::endl;
      }
      sem.LoopBreak();
    } else {
      if (m.IsRoot()) {
        if (st_.step % var.Int("report_step_every", 1) == 0) {
          ReportStep();
        }
      }
      m.SeedSample();
    }
  }

  CheckAbort(sem, ctx->nabort);

  if (sem("updatepar")) {
    if (auto fs = dynamic_cast<Simple<M>*>(fs_.get())) {
      fs->SetPar(ParsePar<Simple<M>>()(var));
    } else if (auto fs = dynamic_cast<Proj<M>*>(fs_.get())) {
      fs->SetPar(ParsePar<Proj<M>>()(var));
    }
    UpdateAdvectionPar();
    fcvm_ = fs_->GetVelocity();
  }
  if (sem.Nested("mixture")) {
    CalcMixture(as_->GetField());
  }
  if (sem.Nested("fs-start")) {
    fs_->StartStep();
  }
  if (sem.Nested("fs-iters")) {
    if (var.Int["enable_fluid"]) {
      StepFluid();
    }
  }
  if (sem.Nested("fs-finish")) {
    fs_->FinishStep();
  }
  if (sem.Nested("as-steps")) {
    if (var.Int["enable_advection"]) {
      StepAdvection();
    }
  }
  if (sem.Nested("tracer-step")) {
    if (tracer_) {
      StepTracer();
    }
  }
  if (sem.Nested("particles-step")) {
    if (particles_) {
      StepParticles();
    }
  }
  if (sem.Nested("stat")) {
    CalcStat();
  }
  if (sem.Nested("stat")) {
    stat_->Update();
  }
  if (sem.Nested()) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity
  }
  if (sem.Nested()) {
    Dump(false);
  }
  if (sem.Nested("stephook")) {
    StepHook(this);
  }
  if (sem("inc")) {
    ++st_.step;
    m.CollectSample("Hydro::Step");
  }
  sem.LoopEnd();

  if (sem.Nested("dumplast")) {
    if (var.Int["dumplast"]) {
      Dump(true);
    }
  }

  if (sem.Nested("posthook")) {
    if (eb_) {
      PostHook(var, fs_->GetVelocity(), m, *eb_);
    } else {
      PostHook(var, fs_->GetVelocity(), m);
    }
  }
}

template <class M>
void Hydro<M>::ReportStep() {
  std::cout << std::fixed << std::setprecision(8) << "STEP=" << st_.step
            << " t=" << st_.t << " dt=" << st_.dt << " ta=" << as_->GetTime()
            << " dta=" << as_->GetTimeStep() << " wt=" << timer_.GetSeconds()
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepAdv() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....adv: t=" << as_->GetTime() << " dt=" << as_->GetTimeStep()
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepTracer() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....tracer: t=" << tracer_->GetTime() << " dt=" << tracer_dt_
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepParticles() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....particles: t=" << particles_->GetTime()
            << " dt=" << particles_dt_ << std::endl;
}

template <class M>
auto Hydro<M>::CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb)
    -> Vect {
  MapEmbed<BCond<Scal>> me_pressure;
  auto& m = eb.GetMesh();
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    const auto nci = bc.nci;
    if (bc.type == BCondFluidType::slipwall ||
        bc.type == BCondFluidType::symm) {
      me_pressure[cf] = BCond<Scal>(BCondType::neumann, nci);
    } else {
      me_pressure[cf] = BCond<Scal>(BCondType::extrap, nci);
    }
  });
  auto fep = UEB::Interpolate(fcp, me_pressure, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
    if (m.IsInner(c)) {
      if (bc.type == BCondFluidType::wall) {
        sum += eb.GetSurface(cf) * fep[cf];
      }
    }
  });
  return sum;
}

template <class M>
auto Hydro<M>::CalcViscousDrag(
    const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
    const Embed<M>& eb) -> Vect {
  auto& m = eb.GetMesh();
  MapEmbed<BCond<Scal>> me_neumann;
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    me_neumann[cf] = BCond<Scal>(BCondType::neumann, bc.nci);
  });
  auto feg = UEB::Gradient(fcvel, fs_->GetVelocityCond(), eb);
  auto femu = UEB::Interpolate(fcmu, me_neumann, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
    if (m.IsInner(c)) {
      if (bc.type == BCondFluidType::wall) {
        sum += feg[cf] * (-eb.GetArea(cf) * femu[cf]);
      }
    }
  });
  return sum;
}

template <class M>
void Hydro<M>::ReportIter() {
  std::cout << std::scientific << std::setprecision(16)
            << ".....iter=" << fs_->GetIter() << ", diff=" << fs_->GetError()
            << std::endl;
}

template <class M>
void Hydro<M>::CheckAbort(Sem& sem, Scal& nabort) {
  if (sem("abort-local")) {
    nabort = 0.;
    try {
      CHECKNAN(as_->GetField(), true)
      CHECKNAN(fs_->GetVelocity(), true)
      CHECKNAN(fs_->GetPressure(), true)
      // check abort TODO: revise,move
      for (auto c : m.Cells()) {
        if (fs_->GetVelocity()[c].sqrnorm() > sqr(var.Double["abortvel"])) {
          std::stringstream g;
          g << "abortvel exceeded at x=" << m.GetCenter(c);
          throw std::runtime_error(g.str());
        }
      }
    } catch (const std::runtime_error& e) {
      std::cout << e.what() << std::endl;
      nabort += 1.;
    }
    m.Reduce(&nabort, "sum");
  }

  if (sem("abort-reduce")) {
    if (nabort != 0.) {
      if (m.IsRoot()) {
        std::cout << "nabort = " << nabort << std::endl;
      }
      sem.LoopBreak();
    }
  }
}

template <class M>
void Hydro<M>::StepFluid() {
  auto sem = m.GetSem("iter"); // sem nested
  if (sem("iter")) {
    OverwriteBc();
  }
  sem.LoopBegin();
  if (sem.Nested("iter")) {
    fs_->MakeIteration();
  }
  if (sem("report")) {
    ++st_.iter;
    if (m.IsLead()) {
      this->var_mutable.Int["iter"] = st_.iter;
    }
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0) {
        ReportIter();
      }
    }
  }
  if (sem("convcheck")) {
    auto it = fs_->GetIter();
    if ((fs_->GetError() < var.Double["tol"] &&
         (int)it >= var.Int["min_iter"]) ||
        (int)it >= var.Int["max_iter"]) {
      sem.LoopBreak();
    }
  }
  // TODO: Suspender loop hangs if (probably) Nested is last
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::StepTracer() {
  auto sem = m.GetSem("tracer-steps"); // sem nested
  sem.LoopBegin();
  if (sem("spawn")) {
    SpawnTracer();
  }
  if (sem.Nested("start")) {
    tracer_->Step(tracer_dt_, fs_->GetVolumeFlux());
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0) {
        ReportStepTracer();
      }
    }
  }
  if (sem("convcheck")) {
    if (tracer_->GetTime() >= fs_->GetTime() - 0.5 * tracer_dt_) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::StepParticles() {
  auto sem = m.GetSem("particles-steps"); // sem nested
  sem.LoopBegin();
  if (sem.Nested("start")) {
    particles_->Step(particles_dt_, fs_->GetVolumeFlux());
  }
  if (sem("spawn")) {
    std::vector<Vect> p_x;
    std::vector<Vect> p_v;
    std::vector<Scal> p_r;
    std::vector<Scal> p_rho;
    std::vector<Scal> p_termvel;
    ParticlesView view{p_x, p_v, p_r, p_rho, p_termvel};
    SpawnParticles(view);
    particles_->Append(view);
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0) {
        ReportStepParticles();
      }
    }
  }
  if (sem("convcheck")) {
    if (particles_->GetTime() >= fs_->GetTime() - 0.5 * particles_dt_) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::StepAdvection() {
  auto sem = m.GetSem("steps"); // sem nested
  sem.LoopBegin();
  if (auto as = dynamic_cast<ASVM*>(as_.get())) {
    const Scal* const voidpenal = var.Double.Find("voidpenal");
    if (voidpenal && sem("void-penal")) {
      auto fccl = as->GetColor();
      auto fcu = as->GetFieldM();
      for (auto f : m.Faces()) {
        const IdxCell cm = m.GetCell(f, 0);
        const IdxCell cp = m.GetCell(f, 1);
        Scal um = 0;
        Scal up = 0;
        for (auto l : layers) {
          if ((*fccl[l])[cm] != kClNone) {
            um += (*fcu[l])[cm];
          }
          if ((*fccl[l])[cp] != kClNone) {
            up += (*fcu[l])[cp];
          }
        }
        um = std::min(1., um);
        up = std::min(1., up);
        FieldFace<Scal>& ffv =
            const_cast<FieldFace<Scal>&>(fs_->GetVolumeFlux().GetFieldFace());
        ffv[f] += -(up - um) * (*voidpenal) * m.GetArea(f);
      }
    }
  }
  if (sem.Nested("start")) {
    as_->StartStep();
  }
  if (sem.Nested("iter")) {
    as_->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as_->FinishStep();
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0) {
        ReportStepAdv();
      }
    }
  }
  if (sem("convcheck")) {
    if (as_->GetTime() >= fs_->GetTime() - 0.5 * as_->GetTimeStep()) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
  if (sem.Nested("as-post")) {
    as_->PostStep();
  }
  if (sem.Nested("curv")) {
    psm_ = UCurv<M>::CalcCurvPart(as_.get(), psm_par_, fck_, m);
  }
  if (var.Int["enable_bubgen"]) {
    if (sem.Nested("bubgen")) {
      StepBubgen();
    }
  }
  if (var.Int["enable_erasevf"]) {
    if (sem("erasevf")) {
      StepEraseVolumeFraction("erasevf", erasevf_last_t_);
      StepEraseVolumeFraction("erasevf2", erasevf2_last_t_);
    }
  }
  if (sem.Nested("erasecl")) {
    StepEraseColor("erasecl");
  }
}

template <class M>
void Hydro<M>::InitTracerFields(Multi<FieldCell<Scal>>& vfcu) {
  auto sem = m.GetSem("tracerfields");
  struct {
    Vars vart;
  } * ctx(sem);
  if (sem("init")) {
    vfcu.resize(var.Int["tracer_layers"]);
    vfcu.InitAll(FieldCell<Scal>(m, 0));
  }
  for (int l = 0; l < var.Int["tracer_layers"]; ++l) {
    const std::string prefix = "tracer" + std::to_string(l);
    if (sem("var" + std::to_string(l))) {
      ctx->vart.String.Set("init_vf", var.String[prefix + "_init"]);
      ctx->vart.String.Set("list_path", var.String[prefix + "_list_path"]);
      ctx->vart.Int.Set("dim", var.Int["dim"]);
      ctx->vart.Int.Set("list_ls", 3);
    }
    if (sem.Nested("field" + std::to_string(l))) {
      InitVf(vfcu[l], ctx->vart, m);
    }
    if (sem("factor" + std::to_string(l))) {
      auto k = var.Double[prefix + "_factor"];
      for (auto c : m.AllCells()) {
        vfcu[l][c] *= k;
      }
    }
  }
  if (sem()) {
  }
}

template <class M>
void Hydro<M>::StepBubgen() {
  auto sem = m.GetSem("bubgen");
  struct {
    FieldCell<Scal> fcvf; // volume fraction
    Vars var;
  } * ctx(sem);
  auto& fcvf = ctx->fcvf;
  const Scal t0 = var.Double["bubgen_t0"];
  const Scal tper = var.Double["bubgen_per"];
  bool bg = (st_.t > t0 && st_.t - bgt_ >= tper);
  if (bg) {
    if (sem("as-bubgen-var")) {
      ctx->var.String.Set("init_vf", "list");
      ctx->var.String.Set("list_path", var.String["bubgen_path"]);
      ctx->var.Int.Set("dim", var.Int["dim"]);
      ctx->var.Int.Set("list_ls", var.Int["list_ls"]);
    }
    if (sem.Nested("as-bubgen-initvf")) {
      InitVf(fcvf, ctx->var, m);
    }
    if (sem("as-bubgen-apply")) {
      bgt_ = st_.t;
      auto apply_vof = [&](auto* as, const auto& eb) {
        if (as) {
          auto& u = const_cast<FieldCell<Scal>&>(as->GetField());
          for (auto c : eb.AllCells()) {
            if (fcvf[c] > 0) {
              u[c] = std::max(u[c], fcvf[c]);
            }
          }
        }
      };
      auto apply_vofm = [&](auto* as, const auto& eb) {
        if (as) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[0]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[0]);
          for (auto c : eb.AllCells()) {
            if (fcvf[c] > 0) {
              u[c] = std::max(u[c], fcvf[c]);
              cl[c] = 1.;
            }
          }
        }
      };
      if (eb_) {
        auto& eb = *eb_;
        for (auto c : m.AllCells()) {
          fcvf[c] = std::min(fcvf[c], eb.GetVolumeFraction(c));
        }
        apply_vofm(dynamic_cast<ASVMEB*>(as_.get()), *eb_);
        apply_vof(dynamic_cast<ASVEB*>(as_.get()), *eb_);
      } else {
        apply_vofm(dynamic_cast<ASVM*>(as_.get()), m);
        apply_vof(dynamic_cast<ASV*>(as_.get()), m);
      }
    }
    if (sem()) {
      // FIXME: empty stage to finish communication to keep ctx
    }
  }
}

template <class M>
void Hydro<M>::StepEraseVolumeFraction(std::string prefix, Scal& last_t) {
  if (!var.Double.Contains(prefix + "_t0")) {
    return;
  }
  const Vect rect_x0(var.Vect[prefix + "_rect_x0"]);
  const Vect rect_x1(var.Vect[prefix + "_rect_x1"]);
  const Rect<Vect> rect(rect_x0, rect_x1);
  const Scal t0 = var.Double[prefix + "_t0"];
  const Scal tper = var.Double[prefix + "_per"];
  if (st_.t > t0 && st_.t - last_t >= tper) {
    if (m.IsRoot() && tper > fs_->GetTimeStep()) {
      std::cout << prefix + " t=" << st_.t << std::endl;
    }
    last_t = st_.t;
    auto apply_vof = [this,&rect](auto* as, const auto& eb) {
      if (as) {
        auto& u = const_cast<FieldCell<Scal>&>(as->GetField());
        for (auto c : eb.AllCells()) {
          const auto x = m.GetCenter(c);
          if (rect.IsInside(x) && u[c] > 0) {
            u[c] = 0;
          }
        }
      }
    };
    auto apply_vofm = [this,&rect](auto* as, const auto& eb) {
      if (as) {
        for (auto l : layers) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[l]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[l]);
          for (auto c : eb.AllCells()) {
            const auto x = m.GetCenter(c);
            if (rect.IsInside(x) && u[c] > 0) {
              u[c] = 0;
              cl[c] = kClNone;
            }
          }
        }
      }
    };
    if (eb_) {
      apply_vofm(dynamic_cast<ASVMEB*>(as_.get()), *eb_);
      apply_vof(dynamic_cast<ASVEB*>(as_.get()), *eb_);
    } else {
      apply_vofm(dynamic_cast<ASVM*>(as_.get()), m);
      apply_vof(dynamic_cast<ASV*>(as_.get()), m);
    }
  }
}

template <class M>
void Hydro<M>::StepEraseColor(std::string prefix) {
  if (!var.Vect.Contains(prefix + "_rect_x0")) {
    return;
  }
  const Vect rect_x0(var.Vect[prefix + "_rect_x0"]);
  const Vect rect_x1(var.Vect[prefix + "_rect_x1"]);
  const Rect<Vect> rect(rect_x0, rect_x1);

  auto sem = m.GetSem();
  struct {
    Scal cl = std::numeric_limits<Scal>::max();
  } * ctx(sem);
  auto& t = *ctx;
  auto apply_vofm = [this, &rect, &t](
                        auto* as, const auto& eb, M& m, Sem& sem) {
    if (as) {
      if (sem()) {
        for (auto l : layers) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[l]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[l]);
          for (auto c : eb.AllCells()) {
            const auto x = m.GetCenter(c);
            if (rect.IsInside(x) && u[c] > 0 && cl[c] != kClNone) {
              t.cl = cl[c];
            }
          }
        }
        m.Reduce(&t.cl, "min");
      }
      if (sem()) {
        for (auto l : layers) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[l]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[l]);
          for (auto c : eb.AllCells()) {
            if (cl[c] == t.cl) {
              u[c] = 0;
              cl[c] = kClNone;
            }
          }
        }
      }
    }
  };
  if (eb_) {
    apply_vofm(dynamic_cast<ASVMEB*>(as_.get()), *eb_, m, sem);
  } else {
    apply_vofm(dynamic_cast<ASVM*>(as_.get()), m, m, sem);
  }
}
