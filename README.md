# Open source codes employing lattice Boltzmann methods
A curated list of some open source frameworks, libraries and softwares employing lattice Boltzmann methods. This list is by no means complete. So if you would like me to add something, please send me a link to sthavishthabr@gmail.com or perform a pull request. 

---------------------------
## Numerical frameworks

* [ch4-project](https://github.com/ecalzavarini/ch4-project) [[1](#ch4-project)] - a parallelized numerical framework written in C++ which can simulate <ins>multiphase flows, turbulence, lagrangian particles (point-sized, tracers etc.)</ins> and <ins>thermal flows</ins>. Developed by the [Unité de Mécanique de Lille](http://uml.univ-lille.fr/fr/presentation-de-lunite), [University of Lille](https://www.univ-lille.fr/home/). 
* [HemeLB](https://github.com/UCL/hemelb/tree/master/Code) [[2](#hemelb)] - a parallelized numerical framework written in C which can simulate flows in <ins>complex geometries including cerebral blood flow scenarios, ranging from normal to neuropathological conditions, aneurysms and arterio-venous malformations</ins>, as well as the <ins>entire intra-cranial vasculature</ins>. Developed at the [Centre for Computational Science](http://ccs.chem.ucl.ac.uk/), [University College London](https://www.ucl.ac.uk/). 
* [HemoCell](https://github.com/UvaCsl/HemoCell) [[3](#hemocell)] - a parallelized numerical framework written in C++ for simulating <ins>deformable blood suspensions</ins> using a combined Immersed boundary-lattice Boltzmann method. This is built on top of the [Palabos](https://palabos.unige.ch/) library. Developed at the [Computational Science lab](https://uva.computationalscience.nl/), [University of Amsterdam](https://www.uva.nl/). 
* [LBDEMcoupling-public](https://github.com/ParticulateFlow/LBDEMcoupling-public) [[4](#LBDEMcoupling)] - a parallelized coupling package between the [Palabos](https://palabos.unige.ch/) library and the discrete element code [LIGGGHTS](https://www.cfdem.com/liggghts-open-source-discrete-element-method-particle-simulation-code) using C++. Can perform <ins>resolved simulations of fluid-particle systems</ins>. Developed by the [Johannes Kepler University Linz](https://www.jku.at/). 
* [LBSim](https://github.com/fabioskomori/lbsim) - a C++ fluid dynamics simulator which can simulate <ins>complex flows, multiphase flows, thermal flows</ins> and <ins>microfluidics</ins>. Developed by the [Software Development Center](http://www.usp.br/nds/) at the [University of Sao Paulo](https://www5.usp.br/). 
* [LBsoft](https://github.com/copmat/LBsoft) [[5](#LBsoft)] - an MPI parallelized numerical solver written in FORTRAN90 for simulating <ins>colloidal emulsions</ins> by coupling lattice Boltzmann methods for the fluid and discrete particle dynamics for the colloids. 
* [LB3D](http://ccs.chem.ucl.ac.uk/lb3d) [[6](#LB3D)] - an MPI parallelized numerical framework written in Fortran90 to simulate <ins>3D simple, binary oil/water</ins> and <ins>ternary oil/water/amphiphile fluids</ins> using the Shan-Chen model for binary fluid interactions. Also offers the possibilities of studying flows involving <ins>micro-mixing, porous media, fluid surface interactions</ins> and other <ins>multicomponent flows</ins>. Developed at the [University College London](https://www.ucl.ac.uk/), [University of Stuttgart](https://www.uni-stuttgart.de/en/) and [Eindhoven University of Technology](https://www.tue.nl/en/). 
* [LIFE](https://github.com/joconnor22/LIFE) - a parallelized numerical framework written in C++ employing lattice Boltzmann-immersed boundary-finite element solver for simulating <ins>fluid-structure interaction problems involving slender structures</ins>. 
* [LUDWIG](https://github.com/ludwig-cf/ludwig) [[7](#LUDWIG_1), [8](#LUDWIG_2)] - a parallelized numerical framework written in C which can simulate <ins>complex fluids including mixtures, symmetric binary fluids, colloidal suspensions, polar gels, charged fluids</ins> and <ins>liquid crystals</ins>. Developed at the [University of Edinburgh](https://www.ed.ac.uk/). 
* [LUMA](https://github.com/aharwood2/LUMA) [[9](#LUMA)] - an MPI parallelized numerical solver written in C/C++ for simulating <ins>3D complex fluid structure interaction</ins> and <ins>turbulent flow problems</ins>. Developed at the [University of Manchester](https://www.manchester.ac.uk/).
* [Microflow 3D](http://www.microflow.pwr.edu.pl/index.php/download2/official-releases) [[10](#Microflow3D_1), [11](#Microflow3D_2)] - an open-source parallel numerical framework (CUDA and GPU) for simulating <ins>transient flows in comlex geometries, inlcuding microflows in channels of microsystems or bio-chips</ins> and <ins>flows in narrow blood vessels</ins>. Developed at the [Wroclaw University of Science and Technology](http://pwr.edu.pl/en/). 
* [MultiphasePorousMediaPalabos](https://github.com/je-santos/MultiphasePorousMediaPalabos#author-s-publications) - an open source parallel numerical framework using the [Palabos](https://palabos.unige.ch/) library for modeling single and multiphase flows in complex porous geometries. Developed at the [University of Texas, Austin](https://www.utexas.edu/). 
* [Musubi](https://hg.osdn.net/view/apes/musubi) [[12](#Musubi)] - a parallel numerical solver written in Fortran 2003 for simulating multicomponent flows, blood flows, porous media and liquid mixtures. Developed at the [Simulation Techniques and Scientific Computing](https://www.mb.uni-siegen.de/sts/) group at [Universität Siegen](https://www.uni-siegen.de/start/). 
* [NekLBM](https://neklbm.mcs.anl.gov/index.php/Main_Page) - a high-order parallelized numerical solver written in Fortran 90 and C, based on spectral element discontinuous Galerkin methods. Developed at the [Mathematics and Computer Science Division](https://www.anl.gov/mcs) of [Argonne National Laboratory](https://www.anl.gov/).
* [OpenLB](https://www.openlb.net/) [[13](#openlb_1), [14](#openlb_2)] - an MPI parallelized numerical framework written in C++ which can simulate <ins>multiphase flows, multicomponent flows, turbulence, particulate flows, complex flows, thermal flows, non-Newtonian flows</ins> and <ins>porous media</ins> to name a few. Developed by the [Lattice Boltzmann research group](http://www.lbrg.kit.edu/) at [Kalsruhe Institute of Technology](https://www.kit.edu/english/). 
* [Palabos](https://palabos.unige.ch/) [[15](#Palabos)] - an MPI parallelized numerical framework written in C++ which can simulate <ins>multiphase flows, multicomponent flows, turbulence, particulate flows, complex flows, thermal flows, non-Newtonian flows</ins> and <ins>grid refinement</ins>. Developed by the [Scientific and Parallel Computing group](https://spc.unige.ch/en/) at [University of Geneva](https://www.unige.ch/). 
* [SailFish](http://sailfish.us.edu.pl/) [[16](#Sailfish)] - a numerical framework written in Python and CUDA C/OpenCL for GPUs (CUDA, OpenCL) which can simulate <ins>complex flows, multicomponent flows</ins> and <ins>turbulence</ins>. 
* [Taxila LBM](https://github.com/ecoon/Taxila-LBM) [[17](#TaxilaLBM_1), [18](#TaxilaLBM_2)] - a parallel numerical framework written in Fortran90 which can simulate flows in </ins>porous and geometrically complex media</ins> (including single and multiphase flows). Developed by the [Los Alamos National laboratory](https://www.lanl.gov/). 
* [TCLB](https://github.com/CFD-GO/TCLB) - a parallel numerical solver (MPI+CUDA or MPI+CPU) for simulating <ins>electrokinetic flows, thermal flows, multiphase flows, non-Newtonian flows</ins> and <ins>shallow water flows</ins>. Developed at the [C-CFD group](https://c-cfd.meil.pw.edu.pl/), [Warsaw University of Technology](https://pw.edu.pl/). 
* [waLBerla](https://www.walberla.net/) [[19](#waLBerla)] - a parallelized (both MPI and hybrid MPI/OpenMP are supported) numerical framework written in C++ supporting <ins>block-structured adaptive mesh refinement, fully resolved particulate flows</ins> and <ins>free surface flows</ins> to name a few. Developed by the [Chair for System Simulation](https://www.cs10.tf.fau.de/), [Friedrich–Alexander–Universität Erlangen–Nürnberg](https://www.fau.eu/).
* [wlb](https://github.com/weierstrass/wlb) [[20](#wlb)] - a github repository hosting C++ codes to simulate <ins>2D electrokinetic flows</ins>, by coupling Navier-Stokes, Nernst-Planck and Poisson’s equation of electrostatics. 
* [2d-lbm-dem](https://github.com/cb-geo/2d-lbm-dem) - a 2D coupled lattice Boltzmann and discrete element method written in C for simulating <ins>granular flows</ins>. Developed by the [Computational Geomechanics group](https://www.cb-geo.com/research/lbm/#main-features), a co-operation between two research groups at the [University of Cambridge](https://www.cam.ac.uk/) and [University of California, Berkeley](https://www.berkeley.edu/). 
 
---------------------------
## Miscellaneous codes

* [EduLB](http://sourceforge.net/projects/edulb/) - an educational C++ code to show the implementation of lattice Boltzmann method by simulating flow over an obstacle in a channel. 
* [LatBo.jl](https://github.com/UCL/LatBo.jl) - an code developed in Julia programming language. 
* [Lattice-Boltzmann-fluid-flow-in-Tensorflow](https://github.com/loliverhennigh/Lattice-Boltzmann-fluid-flow-in-Tensorflow) - a github repository hosting lattice Boltzmann simulation results written in Tensorflow. 
* [LatticeBoltzmannMethod](https://github.com/shurikkuzmin/LatticeBoltzmannMethod) - a github repository hosting some excellent codes (C++) showcasing multiphase flows, microflows and immersed boundary-lattice Boltzmann methods to name a few. 
* [LBMCode](https://www.researchgate.net/publication/299367818_Lattice-Boltzmann_code_for_flow_simulation_in_a_simple_straight_channel) - a FORTRAN90 code solving the shallow water equations to simulate flows in a straight channel. 
* [LBM-1D](https://github.com/sthavishtha/LBM-1D) - a github repository hosting some simple MATLAB codes to simulate 1D advection-diffusion and Navier-Stokes equations. 
* [lbmles](https://github.com/jabhiji/lbmles) - a github repository hosting a 2D lattice Boltzmann code to solve fluid flow in lid-driven cavity at very large Reynolds numbers. Both CPU and GPU (C++ and CUDA) versions are available. 
* [lbm_matlab](https://github.com/rlee32/lbm_matlab) - a github repository hosting some MATLAB codes showcasing grid refinement, viscosity counteraction approach (for achieving numerical stability), RANS Spalart-Allmaras turbulence model and a few more.
* [lbm-principles-practice](https://github.com/lbm-principles-practice/code) - a github repository hosting the codes (C++, MATLAB) used in the lattice Boltzmann book [[21](#lbmbook)]. 
* [3D-LBM-AMR](https://github.com/blackoffee/3D-LBM-AMR) - a github repository hosting C++ codes showcasing adaptive mesh refinement in 3D problems. 

---------------
## References

1. <a name="ch4-project"></a> Calzavarini, E., _Eulerian–Lagrangian fluid dynamics platform: The ch4-project_, Software Impacts, Vol. 1, 2019. [Link](https://doi.org/10.1016/j.simpa.2019.100002)

1. <a name="hemelb"></a> Mazzeo, M.D., Coveney, P.V., _HemeLB: A high performance parallel lattice-Boltzmann code for large scale fluid flow in complex geometries_, Computer Physics Communications, Vol. 178 (12), pp. 894-914, 2008. [Link](https://doi.org/10.1016/j.cpc.2008.02.013)

1. <a name="hemocell"></a> Závodszky, G. et al., _Cellular Level In-silico Modeling of Blood Rheology with An Improved Material Model for Red Blood Cells_, Frontiers in Physiology, Vol. 8, pp. 563, 2017. [Link](https://dx.doi.org/10.3389%2Ffphys.2017.00563)

1. <a name="LBDEMcoupling"></a> Seil, P. and Pirker, S., _LBDEMcoupling: Open-Source Power for Fluid-Particle Systems_, Proceedings of the 7th International Conference on Discrete Element Methods (DEM) 2016, Springer Proceedings in Physics, Vol. 188, 2017. [Link](https://www.researchgate.net/deref/http%3A%2F%2Fdx.doi.org%2F10.1007%2F978-981-10-1926-5_70)

1. <a name="LBsoft"></a> Bonaccorso, F. et al., _LBsoft: a parallel open-source software for simulation of colloidal systems_, arXiv, 2020. [Link](https://arxiv.org/abs/2004.04080)

1. <a name="LB3D"></a> Schmieschek, S. et al., _LB3D: A parallel implementation of the Lattice-Boltzmann method for simulation of interacting amphiphilic fluids_, Computer Physics Communications, Vol. 217, pp. 149-161, 2017. [Link](https://doi.org/10.1016/j.cpc.2017.03.013)

1. <a name="LUDWIG_1"></a> Desplata, J.-C., Pagonabarraga, I. and Bladon, P., _LUDWIG: A parallel Lattice-Boltzmann code for complex fluids_, Computer Physics Communications, Vol. 134 (3), pp. 273-290, 2001. [Link](https://www.sciencedirect.com/science/article/pii/S0010465500002058)

1. <a name="LUDWIG_2"></a> Gray, A. and Stratford, K., _Ludwig: multiple GPUs for a complex fluid lattice Boltzmann application_, Designing Scientific Applications on GPUs, Chapman and Hall/CRC, 2013.

1. <a name="LUMA"></a> Harwood, A.R.G. et al., _LUMA: A many-core, Fluid–Structure Interaction solver based on the Lattice-Boltzmann Method_, SoftwareX, Vol. 7, pp. 88-94, 2018. [Link](https://doi.org/10.1016/j.softx.2018.02.004)

1. <a name="Microflow3D_1"></a> Tomczak, T., Szafran, R., _A new GPU implementation for lattice-Boltzmann simulations on sparse geometries_, Computer Physics Communications, Vol. 235, pp. 258-278, 2019. [Link](https://doi.org/10.1016/j.cpc.2018.04.031)

1. <a name="Microflow3D_2"></a> Tomczak, T., Szafran, R., _Sparse geometries handling in lattice Boltzmann method implementation for graphic processors_, IEEE Transactions on Parallel and Distributed Systems, Vol. 29(8), pp. 1865 - 1878, 2018. [Link](https://doi.org/10.1109/TPDS.2018.2810237)

1. <a name="Musubi"></a> Hasert, M. et al., _Complex fluid simulations with the parallel tree-based Lattice Boltzmann solver Musubi_, Journal of Computational Science, Vol. 5(5), pp. 784-794, 2014. [Link](https://doi.org/10.1016/j.jocs.2013.11.001)

1. <a name="openlb_1"></a> Heuveline, V. and Latt, J., _The OpenLB project: an open source and object oriented implementation of lattice boltzmann methods_, International Journal of Modern Physics C, Vol. 18 (4), pp. 627-634,2007. [Link](https://doi.org/10.1142/S0129183107010875)

1. <a name="openlb_2"></a> Heuveline, V. and Krause, M.J., _OpenLB: Towards an Efficient Parallel Open Source Library for Lattice Boltzmann Fluid Flow Simulations_, PARA'08 Workshop on State-of-the-Art in Scientific and Parallel Computing, May 13-16, 2008. [Link](https://www.researchgate.net/publication/256120850_OpenLB_Towards_an_Efficient_Parallel_Open_Source_Library_for_Lattice_Boltzmann_Fluid_Flow_Simulations)

1. <a name="Palabos"></a> Latt, J. et al., _Palabos: Parallel Lattice Boltzmann Solver_, arXiv, 2019. [Link](https://www.researchgate.net/publication/336349704_Palabos_Parallel_Lattice_Boltzmann_Solver)

1. <a name="Sailfish"></a> Januszewski, M. and Kostur, M., _Sailfish: A flexible multi-GPU implementation of the lattice Boltzmann method_, Computer Physics Communications, Vol. 185 (9), pp. 2350-2368, 2014. [Link](https://doi.org/10.1016/j.cpc.2014.04.018)

1. <a name="TaxilaLBM_1"></a> Coon, E.T., Porter, M.L. and Kang, Q, _Taxila LBM: a parallel, modular lattice Boltzmann framework for simulating pore-scale flow in porous media_. Computational Geosciences, Vol. 18, pp. 17–27, 2014. [Link](https://doi.org/10.1007/s10596-013-9379-6)

1. <a name="TaxilaLBM_2"></a> Porter, M.L. et al., _Multicomponent interparticle-potential lattice Boltzmann model for fluids with large viscosity ratios_, Physical Review E, Vol. 86 (3), 036701, 2012. [Link](https://doi.org/10.1103/PhysRevE.86.036701)

1. <a name="waLBerla"></a> Bauer, M. et al., _waLBerla: A block-structured high-performance framework for multiphysics simulations_, To appear in Computers & Mathematics with Applications, 2020. [Link](https://doi.org/10.1016/j.camwa.2020.01.007)

1. <a name="wlb"></a> Bülling, A., _Modelling of electrokinetic flow using the lattice-Boltzmann method_, Master thesis, Chalmers University of Science and Technology, 2012. [Link](http://publications.lib.chalmers.se/records/fulltext/170015/170015.pdf)

1. <a name="lbmbook"></a> Krüger, T. et al., _The Lattice Boltzmann Method: Principles and Practice_, Springer International Publishing, ISBN 978-3-319-44647-9, 2017. [Link](https://www.springer.com/gp/book/9783319446479)

