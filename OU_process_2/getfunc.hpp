#include <iostream>
#include <memory> 
#include <cstring>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <cln/cln.h>
#include <random>
#include <new>
#include <cmath>
#include <cassert>
#include <typeinfo>
#include "get_vol.hpp"
#include <omp.h>
#include <stack>
#include <ginac/ginac.h>
#include <vector>
#include <utility>

GiNaC::symbol get_x0();

GiNaC::symbol get_x1();

GiNaC::symbol get_x2();


std::vector<GiNaC::ex> getfunc();
