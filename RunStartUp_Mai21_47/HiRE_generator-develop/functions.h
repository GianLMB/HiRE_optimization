#pragma once

#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>

namespace OpepFunctions
{

	size_t index_of(std::vector<std::string> vec, std::string elementToSearch);

	std::string getFortranFloat(double n);

	// template function must be implementd in header
	template<typename Type>
	Type deg2rad(Type angleDeg)
	{
		return angleDeg / 180 * M_PI;
	}

	template<typename Type>
	Type rad2deg(Type angleRad)
	{
		return angleRad * 180 / M_PI;
	}

}
