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
#include <RInside.h>
#include <Rcpp.h>


//int argc=0;
//char **argv=NULL;
//RInside V(argc,argv);
//Eigen::VectorXd Z=get_vol(&V);
//int run=Z.size()-1;

GiNaC::symbol get_x0(){
	GiNaC::symbol x0("x0");
	return x0;
};

GiNaC::symbol get_x1(){
	GiNaC::symbol x1("x1");
	return x1;
};

GiNaC::symbol get_x2(){
	GiNaC::symbol x2("x2");
	return x2;
};

std::vector<GiNaC::ex>  getfunc(){
	GiNaC::symbol x0=get_x0();
	GiNaC::symbol x1=get_x1(); 
	GiNaC::symbol x2=get_x2();
	int argc=0;
	char **argv=NULL;
	RInside V(argc,argv);
	Eigen::VectorXd Z=get_vol(&V);
	//Z=Z.head(2);
	int run=Z.size()-1;


	GiNaC::ex hold_f=0;
	GiNaC::ex hold_grad_f_1=0;
	GiNaC::ex hold_grad_f_2=0;
	GiNaC::ex hold_grad_f_3=0;
 	GiNaC::ex hold_h_1=0;
	GiNaC::ex hold_h_2=0;
	GiNaC::ex hold_h_3=0;
       	GiNaC::ex hold_h_4=0;
	GiNaC::ex hold_h_5=0;
	GiNaC::ex hold_h_6=0;
	double dt=1.0/252.0;
	GiNaC::ex  e2kt=exp(-2*(x2)*dt);
	GiNaC::ex ekt=exp(-1*(x2)*dt);
	GiNaC::ex sigma_bar=(x1)*(GiNaC::sqrt((1-e2kt)/(2*(x2))));
	GiNaC::ex sqrt_term=GiNaC::sqrt((1-e2kt)/(2*(x2)));
	GiNaC::ex dsigma_bar_dk=((x1)/2)*(1/(GiNaC::sqrt((1-e2kt)/(2*(x2))))*(4*(x2)*dt*e2kt-2*(1-e2kt)));
	 GiNaC::ex dsigma_bar2_dk2=-1*((x1)/4)*(GiNaC::pow((1-e2kt)/(2*(x2)),-1.5))*(GiNaC::pow((((4*(x2)*dt*e2kt)-2*(1-e2kt))/(4*GiNaC::pow((x2),2))),2))+(((4*GiNaC::pow((x1),2))*(4*dt*e2kt-8*(x2)*(GiNaC::pow(dt,2))*e2kt+4*dt*e2kt)-(8*(x2))*(4*(x2)*dt*e2kt-2*(1-e2kt)))/(16*GiNaC::pow((x2),4))*(GiNaC::pow((1-e2kt)/(2*(x2)),-0.5))*((x1)/2));
	  GiNaC::ex std=sigma_bar;
	  GiNaC::ex dstdk=dsigma_bar_dk;
	 GiNaC::ex dstds=GiNaC::sqrt((1-e2kt)/(2*(x2)));


	  hold_h_6+=(run)*(-1*GiNaC::pow(dsigma_bar_dk/sigma_bar,2)+(dsigma_bar2_dk2/sigma_bar));

	  #pragma omp paralell for reudction(+:hold_f, hold_grad_f_1, hold_grad_f_2, hold_grad_f_3, hold_h_1, hold_h_2, hold_h_3, hold_h_4, hold_h_5, hold_h_6)
	  for(int i=0; i<run; ++i){
		GiNaC::ex u=Z(i)*(x0)*ekt*(1-ekt);
		GiNaC::ex duk=Z(i)*(x0)*dt*(ekt+2*e2kt);
		GiNaC::ex dut=Z(i)*ekt*(1-ekt);
		GiNaC::ex mu=u;	
		GiNaC::ex dmu_dk=-1*Z(i)*(x0)*dt*(ekt+0.5*e2kt);
		GiNaC::ex dmu_dt=Z(i)*ekt*(1-ekt);
		GiNaC::ex dmu2_dtdk=-1*(Z(i))*dt*(ekt+2*e2kt);
		GiNaC::ex dmu2_dk2=Z(i)*(x0)*dt*(ekt+0.25*e2kt);
		//-= due to -1*hold in orginal
		hold_f+=0.5*(GiNaC::log(2*GiNaC::Pi))+GiNaC::log(std)+(GiNaC::pow(Z(i+1)-u,2)/(2*GiNaC::pow(std,2)));
		hold_grad_f_3-=(-1/(std))*dstdk+(duk*(Z(i+1)-u))/(GiNaC::pow(std,2))+(dstdk*(GiNaC::pow((Z(i+1)-u),2)))/(GiNaC::pow(std,3));
		hold_grad_f_2-=(-1/std)+dstds*((GiNaC::pow(Z(i+1),2)*dstds)/(GiNaC::pow(std,3)));
	 	hold_grad_f_1+=((Z(i+1)-u)/(GiNaC::pow(std,2)))*dut;
		hold_h_1+=GiNaC::pow((Z(i)*ekt*(1-ekt))/(sigma_bar),2);
		hold_h_2+=(2*(Z(i+1)-mu)*ekt*(1-ekt)*sqrt_term*Z(i))/(GiNaC::pow(sigma_bar,3));
		hold_h_3+=(-1/(GiNaC::pow((x1),2)))+3*(GiNaC::pow((Z(i+1)-mu)*sqrt_term/(GiNaC::pow(sigma_bar,2)),2));
		hold_h_4+=dmu_dk*dmu_dt+dmu2_dtdk+(2/(x1))*dsigma_bar_dk*(Z(i+1)-mu)*dmu_dt;
		hold_h_5+=(((2*GiNaC::pow(sigma_bar,3)*(Z(i+1)-mu)*dmu_dk+3*(GiNaC::pow(sigma_bar*(Z(i+1)-mu),2)*dsigma_bar_dk)/(GiNaC::pow(sigma_bar,6)))))*sqrt_term-(1/(4*GiNaC::pow((x2),2)*sqrt_term))*(4*(x2)*dt*e2kt-2*(1-e2kt))*GiNaC::pow((Z(i+1)-mu)/(GiNaC::pow(sigma_bar,2)),2);
		hold_h_6-=GiNaC::pow(dmu_dk,2)+dmu2_dk2*(Z(i+1)-mu);
		hold_h_6-=-1*(GiNaC::pow((dsigma_bar_dk*(Z(i+1)-mu))/(sigma_bar),2))+(dsigma_bar2_dk2/sigma_bar)*(GiNaC::pow((Z(i+1)-mu),2))-(2/sigma_bar)*(Z(i+1)-mu)*dsigma_bar_dk*dmu_dk;
		//hold_f.print(GiNaC::print_tree(std::cout));	
	  };	
	  std::vector<GiNaC::ex> hold_1{hold_f, hold_grad_f_1, hold_grad_f_2, hold_grad_f_3, hold_h_1, hold_h_2, hold_h_3, hold_h_4, hold_h_5, hold_h_6};
//	std::vector<GiNaC::symbol> hold_2{(x0),(x1),(x2)};
//	std::pair<std::vector<GiNaC::ex>,std::vector<GiNaC::symbol>> hold(hold_1,hold_2);
	//std::cout<<"check "<<hold_f.normal().is_zero()<<'\n';
	//std::cout<<"Z len " <<run<<'\n';    
	hold_f.print(GiNaC::print_latex(std::cout));	
	return hold_1;
};

