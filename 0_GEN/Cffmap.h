

#pragma once


#include	"CTForm.h"

#include	<string>
using namespace std;


class ffmap {  // set of maps from one file to another

public:
	string			fname;		// name of image mapped to
	string			mname;		// name of mapping file
	vector<Point>	centers;	// centers of sub-maps
	vector<TForm>	transforms;	// transform for each portion
};
