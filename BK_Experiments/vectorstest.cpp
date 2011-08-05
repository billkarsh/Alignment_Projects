

#include	<stdio.h>

#include	<vector>
using namespace std;






int main(int argc, char **argv)
{
	vector<int>	v0(50);
	printf( "v0: size=%d, cap=%d.\n", v0.size(), v0.capacity() );

	v0.clear();
	v0.reserve(55);
	v0.push_back( 5 );
	printf( "v0: size=%d, cap=%d.\n", v0.size(), v0.capacity() );

	vector<int>	v = v0;
	printf( "v: size=%d, cap=%d.\n", v.size(), v.capacity() );

	return 0;
}


