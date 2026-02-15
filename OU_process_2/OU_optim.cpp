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
#include "OU_optim.hpp"
#include <coin-or/IpReferenced.hpp>
#include <typeinfo>
#include "get_vol.hpp"
#include <omp.h>
#include <gmpxx.h>
//#include <mpfr.h>
#include <mpreal.h>
#include <stack>
#include <ginac/ginac.h>
#include "getfunc.hpp"
#include <string>
#include <vector>
#include <map>
#include "map_ex.hpp"
#include <utility>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


std::vector<GiNaC::ex> hold=getfunc();
//GiNaC::ex hold_f=hold[0];
//GiNaC::ex hold_f_grad_1=hold[1];
//GiNaC::ex hold_f_grad_2=hold[2];
//GiNaC::ex hold_f_grad_3=hold[3];
//GiNaC::ex hold_h_1=hold[4];
//GiNaC::ex hold_h_2=hold[5];
//GiNaC::ex hold_h_3=hold[6];
//GiNaC::ex hold_h_4=hold[7];
//GiNaC::ex hold_h_5=hold[8];
//GiNaC::ex hold_h_6=hold[9];
//extern GiNaC::symbol x0;
//extern GiNaC::symbol x1;
//extern GiNaC::symbol x2;





//std::map<int, GiNaC::ex> M=map_ex(hold_f_grad_1);

//


//std::cout<<Z<<'\n';
using namespace Ipopt;

OU_optim::OU_optim(){ }

OU_optim::~OU_optim(){ }

bool OU_optim::get_nlp_info(
            Index&          n,
            Index&          m,
            Index&          nnz_jac_g,
            Index&          nnz_h_lag,
            IndexStyleEnum& index_style
        )
        {
            n=3;
            
            m=0;

            nnz_jac_g=0;

            nnz_h_lag = 6;

            index_style=TNLP::C_STYLE;

            return true;
        };
bool OU_optim::get_bounds_info(
            Index   n,
            Number* x_lower,
            Number* x_upper,
            Index   m,
            Number* g_l,
            Number* g_u
        )
        {
            assert(n==3);
            assert(m==0);

	    //theta is x0 sigma is x1 and kappa is x2
            x_lower[0] = -2e19;

            x_lower[1]=1e-12;

	    x_lower[2]=1e-12;

            x_upper[0] = 2e19;

            x_upper[1] = 2e19;

	    x_upper[2]=2e19;




            return true;
        };
bool OU_optim::get_starting_point(
            Index   n,
            bool    init_x,
            Number* x,
            bool    init_z,
            Number* z_L,
            Number* z_U,
            Index   m,
            bool    init_lambda,
            Number* lambda
        )
        {
            assert(init_x==true);
            assert(init_z==false);
            assert(init_lambda==false);
            x[0]=5;
            x[1]=1;
	    x[2]=1;
            return true;
        };
bool OU_optim::eval_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number&       obj_value
        )
        {
            assert(n==3);
	    
	   // std::vector<GiNaC::ex> hold=getfunc();
	    GiNaC::ex hold_f=hold[0];
            obj_value=GiNaC::ex_to<GiNaC::numeric>(hold_f.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	    //std::cout<<obj_value<<'\n';
	    return true;
        };
bool OU_optim::eval_grad_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number*       grad_f
        ){
            assert(n==3); 
            
	    //std::vector<GiNaC::ex> hold=getfunc();
	    GiNaC::ex hold_f_grad_1=hold[1];  
	    GiNaC::ex hold_f_grad_2=hold[2]; 
	    GiNaC::ex hold_f_grad_3=hold[3];  
	    //std::cout<<hold_f_grad_1.is_zero()<<'\n';
	    grad_f[2]=GiNaC::ex_to<GiNaC::numeric>(hold_f_grad_3.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();

            grad_f[1]=GiNaC::ex_to<GiNaC::numeric>(hold_f_grad_2.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();

            grad_f[0]=GiNaC::ex_to<GiNaC::numeric>(hold_f_grad_1.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();		
	    //std::cout<<hold_f_grad_1.has(get_x0())<<'\n';
	    //std::cout<<hold_f_grad_1.is_zero()<<'\n';
	    //std::cout<<"grad_f_1 "<<grad_f[0]<<'\n';
	   // std::cout<<"grad_f_1 "<<GiNaC::ex_to<GiNaC::numeric>(hold_f_grad_1.subs(GiNaC::lst{(get_x0())==x[1], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double()<<'\n';
	    return true;
        };
bool OU_optim::eval_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Number*       g
        )
        {
            assert(n==3);
            assert(m==0);

            return true;
        };

bool OU_optim::eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
        )
        {
            assert(n==3);
            assert(m==0);

            return true;
        };
bool OU_optim::eval_h(
        Index         n,
        const Number* x,
        bool          new_x,
        Number        obj_factor,
        Index         m,
        const Number* lambda,
        bool          new_lambda,
        Index         nele_hess,
        Index*        iRow,
        Index*        jCol,
        Number*       values
        )
        {
        assert(n == 3);
        assert(m == 0);
        
        if( values == NULL )
        {
            Index idx = 0;
            for( Index row = 0; row < 3; row++ )
            {
                for( Index col = 0; col <= row; col++ )
                {
                    iRow[idx] = row;
                    jCol[idx] = col;
                    idx++;
                }
            }
        
            assert(idx == nele_hess);
        }
        else{
  	   
	    //std::vector<GiNaC::ex> hold=getfunc();
	    GiNaC::ex hold_h_1=hold[4];
	    GiNaC::ex hold_h_2=hold[5];
	    GiNaC::ex hold_h_3=hold[6];
	    GiNaC::ex hold_h_4=hold[7];
	    GiNaC::ex hold_h_5=hold[8];
	    GiNaC::ex hold_h_6=hold[9];


	    double hold_0=GiNaC::ex_to<GiNaC::numeric>(hold_h_1.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	    values[0]=obj_factor*hold_0;

	    double hold_1=GiNaC::ex_to<GiNaC::numeric>(hold_h_2.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	    values[1]=obj_factor*hold_1;
	    double hold_2=GiNaC::ex_to<GiNaC::numeric>(hold_h_3.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	    values[2]=obj_factor*hold_2;
	   double hold_3=GiNaC::ex_to<GiNaC::numeric>(hold_h_4.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	   values[3]=obj_factor*hold_3;
	   double hold_4=GiNaC::ex_to<GiNaC::numeric>(hold_h_5.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
  	   values[4]=obj_factor*hold_4;
	   double hold_5=GiNaC::ex_to<GiNaC::numeric>(hold_h_6.subs(GiNaC::lst{(get_x0())==x[0], (get_x1())==x[1], (get_x2())==x[2]}).evalf()).to_double();
	   values[5]=obj_factor*hold_5;

     
	   //for(int i=0; i<6; ++i){
	 	//if(!std::isfinite(values[i])){
			//if(values[i]>0){
				//values[i]=0*obj_factor;
			//}
			//else{
		//		values[i]=0*obj_factor;
		//	}

		//	};
	//	};
	    };
	    return true;
        
	};

void OU_optim::finalize_solution(
        SolverReturn               status,
        Index                      n,
        const Number*              x,
        const Number*              z_L,
        const Number*              z_U,
        Index                      m,
        const Number*              g,
        const Number*              lambda,
        Number                     obj_value,
        const IpoptData*           ip_data,
        IpoptCalculatedQuantities* ip_cq
        )
        {
        std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
        for( Index i = 0; i < n; i++ )
        {
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }
        };

