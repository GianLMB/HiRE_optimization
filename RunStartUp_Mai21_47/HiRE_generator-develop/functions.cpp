#include "functions.h"

namespace OpepFunctions
{

	#ifndef M_PI
	   #define M_PI 3.14159265358979323846
	#endif

	using namespace std;

	size_t index_of(vector<string> vec, string elementToSearch)
	{
		for(size_t i=0; i<vec.size(); ++i) {
			if(vec[i] == elementToSearch) {
				return i;
			}
		}
		return -1;
	}

	string getFortranFloat(double n)
	{

		int p;
		if(n == 0) {
			p = 1;
		} else {
			p = floor(log10(fabs(n))) + 1; }

		double number = n/pow(10, p);

		ostringstream stringStream;
		stringStream << " ";
		stringStream << ((number<0)?"-":" ");
		stringStream << fixed << setprecision(8) << fabs(number) ;
		stringStream << "E";
		stringStream << ((p<0)?"-":"+");
		stringStream << setfill('0') << setw(2)<< abs(p);

		return stringStream.str();
	}

}
