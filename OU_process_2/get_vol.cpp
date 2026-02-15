#include <RInside.h>
#include <Rcpp.h>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <iostream>

Eigen::MatrixXd get_vol(RInside* R){
	RInside R_env=*R;
	std::string r_code="library(quantmod);\ngetSymbols('^GSPC', src='yahoo', from='2003-01-07', to='2022-03-07');\nlog_returns <- dailyReturn(GSPC, type='log');\nsquared_log_returns <- log_returns^2;\nrolling_sd <- rollapply(squared_log_returns, width=40, FUN=sd, align='right');\noutput <-na.omit(rolling_sd);\noutput_vector=as.vector(output);\nrows <-length(output_vector);\n";
	//std::string r_code_test=".libPaths();\n";
	
	R_env.parseEvalQ(r_code);

	
	Rcpp::NumericVector r_vector=Rcpp::as<Rcpp::NumericVector>(R_env["output_vector"]);
	//int cols= Rcpp::as<int>(R_env["cols"]);
	int rows= Rcpp::as<int>(R_env["rows"]);
	Eigen::VectorXd hold(rows);
	//std::cout<<cols<<" "<<rows<<'\n';
	for(int i=0; i<rows; i++)
	{
		hold(i)=r_vector(i);
	}
//	std::cout<<hold<<'\n';
	return hold;
};
		
		
				

