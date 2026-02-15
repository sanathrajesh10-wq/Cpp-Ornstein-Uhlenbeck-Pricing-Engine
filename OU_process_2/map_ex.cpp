#include <cln/cln.h>
#include <ginac/ginac.h>
#include <string>
#include <map>


std::map<int, GiNaC::ex>  map_ex(GiNaC::ex hold){
	std::map<int, GiNaC::ex> m;

	GiNaC::exset symset;

	hold.find(GiNaC::symbol(), symset);

	for(auto const &s: symset){
		std::string name=GiNaC::ex_to<GiNaC::symbol>(s).get_name();
		if(name=="x0"){
			m[0]=s;
		};
		if(name=="x1"){
			m[1]=s;
		};
		if(name=="x3"){
			m[2]=s;
		};
	};
	return m;
};	
