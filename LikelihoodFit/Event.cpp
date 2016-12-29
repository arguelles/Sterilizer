#include "Event.h"

double differenceAngle(double zenith1, double azimuth1, double zenith2, double azimuth2){
	using boost::math::constants::pi;
	
	double cad=cos(azimuth1-azimuth2);
	double czd=cos(zenith1-zenith2);
	double czs=cos(zenith1+zenith2);
	double dot=(cad*(czd-czs)+czd+czs)/2;
	if(dot>1.)
		return(0.0);
	if(dot<-1.)
		return(pi<double>());
	return(acos(dot));
}

herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data){
	std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
	items->insert(member_name);
	return(0);
}
