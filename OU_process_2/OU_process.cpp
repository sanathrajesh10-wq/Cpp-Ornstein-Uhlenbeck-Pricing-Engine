#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <random>
#include <new>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "OU_optim.hpp"
#include <coin-or/IpIpoptApplication.hpp>
#include <coin-or/IpSmartPtr.hpp>
#include <cfenv>
using namespace Ipopt;
int main( 
    int,
   char**
)
{
TNLP* mynlp = new OU_optim();
SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
app->Options()->SetNumericValue("tol", 3.82e-6);
app->Options()->SetStringValue("derivative_test","second-order");
app->Options()->SetStringValue("derivative_test_print_all", "yes");
app->Options()->SetNumericValue("derivative_test_tol",1e-4);
app->Options()->SetStringValue("hessian_approximation", "limited-memory");
//app->Options()->SetNumericValue("bound_push",1e-2);
//app->Options()->SetNumericValue("bound_relax_factor",1e-8);
//app->Options()->SetIntegerValue("print_level",12);
app->Options()->SetStringValue("check_derivatives_for_naninf","yes");
app->Options()->SetStringValue("mu_strategy", "adaptive");
app->Options()->SetStringValue("output_file", "ipopt.out");
ApplicationReturnStatus status;
//feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
status = app->Initialize();
if( status != Solve_Succeeded )
{
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
}
status = app->OptimizeTNLP(mynlp);
 
if( status == Solve_Succeeded )
{
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
}
else
{
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
}


return (int) status;





















}









    
