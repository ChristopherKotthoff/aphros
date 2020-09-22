#!/usr/bin/env python3

import re
import os

text = '''<!--- Generated by [[GEN]] -->
# Aphrós [<img src="https://circleci.com/gh/cselab/aphros.svg?style=svg">](https://github.com/cselab/aphros/commits/master)

<img src="doc/images/foam.png" width=300 align="right">

Finite volume solver for incompressible multiphase flows with surface tension.

Key features:

- implementation in C++14
- scalability to thousands of compute nodes with the
  [Cubism](https://gitlab.ethz.ch/mavt-cse/Cubism)
  library for distributed computing on structured grids
- coroutines to enable encapsulation in the block-wise processing framework
- fluid solver based on SIMPLE or Bell-Colella-Glaz methods
- conservative split PLIC advection solver
- novel particle method for curvature estimation improving the accuracy at low resolutions
[[demo]](https://cselab.github.io/aphros/curv.html)
[[ref:partstr]]

### Clone

    git clone https://github.com/cselab/aphros.git

### Documentation

<https://cselab.github.io/aphros/doc>

Generated in [doc/sphinx](doc/sphinx).

### Requirements

C++14, cmake, MPI, hdf5 (parallel), python3, python3-numpy.

Bundled dependencies:

* [hypre](https://github.com/hypre-space/hypre)
* [eigen](https://gitlab.com/libeigen/eigen) (optional)
* [overlap](https://github.com/severinstrobl/overlap) (optional)
* [vofi](https://github.com/VOFTracking/Vofi) (optional)
* [fpzip](https://github.com/LLNL/fpzip) (optional)

### Build and install

Follow instructions from [deploy/README.md](deploy/README.md) to
prepare environment and install dependencies.

Configure, build, install and run tests:

    cd src
    make -j4
    make test

### Docker

Instead of building the code in your system, you can build and run a Docker
container

    docker build github.com/cselab/aphros --tag aphros
    docker run -i aphros

## Videos

Examples of simulations visualized using
[ParaView](https://www.paraview.org/) and [OSPRay](https://www.ospray.org/)
in collaboration with Jean M. Favre at [CSCS](https://www.cscs.ch).

|    |    |
:---:|:---:
[<img src="doc/images/coalescence.jpg" width=384>]([[VIDEOS]]/coalescence.mp4) | [<img src="doc/images/taylor_green.jpg" width=250>]([[VIDEOS]]/taylor_green.mp4)
Coalescence of bubbles [[ref:partstr]] | Taylor-Green vortex with bubbles [[ref:pasc19]] [[ref:datadriven]]
[<img src="doc/images/vortex_bubble.jpg" width=200>]([[VIDEOS]]/vortex_bubble.mp4) | [<img src="doc/images/plunging_jet.jpg" width=200>]([[VIDEOS]]/plunging_jet.mp4)
Bubble trapped by vortex ring [[ref:datadriven]] | Plunging jet [[ref:pasc19]]
[<img src="doc/images/reactor.jpg" width=384>]([[VIDEOS]]/reactor.mp4) | [<img src="doc/images/mesh_bubbles.jpg" width=384>]([[VIDEOS]]/mesh_bubbles.mp4)
Electrochemical reactor [[ref:ees]] | Bubbles through mesh
[<img src="doc/images/rising_bubbles.jpg" width=384>]([[VIDEOS]]/rising_bubbles.mp4) | [<img src="doc/images/foaming_waterfall.jpg" width=384>]([[VIDEOS]]/foaming_waterfall.mp4)
 Rising bubbles clustering on the surface [[ref:aps]] [[ref:cscs]] | Foaming waterfall [[ref:pasc20]]

|     |
|:---:|
|[<img src="doc/images/breaking_waves.jpg" width=795>](https://www.youtube.com/watch?v=iGdphpztCJQ)|
|APS Gallery of Fluid Motion 2019 award winner: Breaking waves: to foam or not to foam? [[ref:aps]]|

## Developers

Aphros is developed and maintained by researchers at ETH Zurich

* [Petr Karnakov](https://www.cse-lab.ethz.ch/member/petr-karnakov/)
* [Dr. Sergey Litvinov](https://www.cse-lab.ethz.ch/member/sergey-litvinov/)
* [Fabian Wermelinger](https://www.cse-lab.ethz.ch/member/fabian-wermelinger/)

under the supervision of

* [Prof. Petros Koumoutsakos](https://www.cse-lab.ethz.ch/member/petros-koumoutsakos/)

## Publications

[[item:ees]] S. M. H. Hashemi, P. Karnakov, P. Hadikhani, E. Chinello, S.
  Litvinov, C.  Moser, P. Koumoutsakos, and D. Psaltis, "A versatile and
  membrane-less electrochemical reactor for the electrolysis of water and
  brine", _Energy & environmental science_, 2019
  [10.1039/C9EE00219G](https://doi.org/10.1039/C9EE00219G)
[[item:pasc19]] P. Karnakov, F. Wermelinger, M. Chatzimanolakis, S. Litvinov,
  and P.  Koumoutsakos, "A high performance computing framework for multiphase,
  turbulent flows on structured grids" in _Proceedings of the platform for
  advanced scientific computing conference on – PASC ’19_, 2019
  [10.1145/3324989.3325727](https://doi.org/10.1145/3324989.3325727)
[[item:icmf]] P. Karnakov, S. Litvinov, P. Koumoutsakos
  "Coalescence and transport of bubbles and drops"
  _10th International Conference on Multiphase Flow (ICMF)_, 2019
[[item:partstr]] P. Karnakov, S. Litvinov, and P. Koumoutsakos, "A hybrid
  particle volume-of-fluid method for curvature estimation in multiphase
  flows”, _International journal of multiphase flow_, 2020
  [10.1016/j.ijmultiphaseflow.2020.103209](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103209)
[[item:datadriven]] Z. Wan, P. Karnakov, P. Koumoutsakos, T. Sapsis, "Bubbles in
  Turbulent Flows: Data-driven, kinematic models with history terms”,
  _International journal of multiphase flow_, 2020
  [10.1016/j.ijmultiphaseflow.2020.103286](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103286)
[[item:aps]] P. Karnakov, S. Litvinov, J. M. Favre, P. Koumoutsakos
  "V0018: Breaking waves: to foam or not to foam?"
  _Gallery of Fluid Motion Award_
  [10.1103/APS.DFD.2019.GFM.V0018](https://doi.org/10.1103/APS.DFD.2019.GFM.V0018)
[[item:cscs]] Annual report 2019 of the Swiss National Supercomputing Centre (cover page)
  [[link]](https://www.cscs.ch/publications/annual-reports/cscs-annual-report-2019)
[[item:pasc20]] P. Karnakov, F. Wermelinger, S. Litvinov,
  and P.  Koumoutsakos, "Aphros: High Performance Software for Multiphase Flows with Large Scale
  Bubble and Drop Clusters" in _Proceedings of the platform for
  advanced scientific computing conference on – PASC ’20_, 2020
  [10.1145/3394277.3401856](https://doi.org/10.1145/3394277.3401856)
'''

m_refs = list(re.finditer("\[\[ref:[^]]*\]\]", text))
m_items = list(re.finditer("\[\[item:[^]]*\]\]", text))
refs = [m_ref.group(0) for m_ref in m_refs]

gen = text
gen = gen.replace('[[GEN]]', os.path.basename(__file__))
gen = gen.replace('[[VIDEOS]]', "https://cselab.github.io/aphros/videos")
for i, m_item in enumerate(m_items):
    item = m_item.group(0)
    name = re.match("\[\[item:([^]]*)\]\]", item).group(1)
    ref = "[[ref:{}]]".format(name)
    start = m_item.start(0)
    end = m_items[i + 1].start(0) if i + 1 < len(m_items) else len(text)
    m_url = re.search("\((http[^)]*)\)", text[start:end])

    gen = gen.replace(item, "{:}.".format(i + 1))

    if m_url:
        url = m_url.group(1)
        gen = gen.replace(ref, "[[{:}]]({})".format(i + 1, url))
    else:
        if ref in refs:
            print("Warning: no URL found for '{}'".format(item))
        gen = gen.replace(ref, "[{:}]".format(i + 1))

with open("README.md", 'w') as f:
    f.write(gen)
