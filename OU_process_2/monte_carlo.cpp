#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <random>
#include <new>
#include <cmath>
#include <typeinfo> 
using namespace Eigen;
int main(){
    double S=101.5;
    double K=99.1;
    double vol=0.0991;
    double r=0.03;
    double N=10;
    double M=100;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::normal_distribution<> d(0.0, 1.0);
    auto normal = [&] (double) {return d(mt);};
    std::cout << typeid(rd).name()<<'\n';
    Eigen::MatrixXd Z =Eigen::MatrixXd::NullaryExpr(M, N, [&] (double) {return d(mt);});
    double T=60.0/365;
    double dt=T/N;
    auto nudt=[&] (double) {return (r-std::pow(vol,2))*dt;};
    double voldt=vol*std::sqrt(dt);
    Eigen::MatrixXd U =Eigen::MatrixXd::NullaryExpr(M, N, nudt);
    Eigen::MatrixXd delta=U+voldt*Z; 
    
    
    
    
    Eigen::VectorXd hold=Eigen::VectorXd::Zero(M);
    for(int i=0; i<N; ++i){
        hold = hold+delta.col(i);    
        };
    Eigen::ArrayXd hold_array=((S*(hold.array()).exp())-K).max(0);
    double C0 = std::exp(-r*T)*(hold_array.sum())/M;
    double sigma= std::sqrt((((hold.array())-C0).pow(2)).sum()/(M-1));
    //std::default_random_engine generator;
    //std::poisson_distribution<int> distribution(4.1);
    //auto poisson = [&] (int) {return distribution(generator);};
    //auto hold= [&] (int) {return 1;};;

    //RowVectorXi v = RowVectorXi::NullaryExpr(10, poisson);
    std::cout << '$' << C0 << '\n';
    std::cout<<sigma<<'\n';    
    return 0;
}