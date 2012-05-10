

#include	"File.h"

#include	<stdio.h>

#include	<set>
#include	<string>
using namespace std;






int main(int argc, char **argv)
{
	set<string>	used;
	FILE		*f = FileOpenOrDie( "usedseries.log", "r" );
	CLineScan	LS;

	for(;;) {

		char	buf[512];

		if( LS.Get( f ) <= 0 )
			break;

		sscanf( LS.line, "%s", buf );
		used.insert( string(buf) );
	}

	fclose( f );



	set<string>	notused;

	f = FileOpenOrDie( "allseries.txt", "r" );

	for(;;) {

		char	buf[512];

		if( LS.Get( f ) <= 0 )
			break;

		sscanf( LS.line, "%s", buf );
		string	s( buf );

		if( used.find( s ) != used.end() )
			;
		else
			notused.insert( s );
	}

	fclose( f );



	set<string>::iterator	it = notused.begin();
	int						nd = notused.size();

	f = FileOpenOrDie( "notused.txt", "w" );

	for( int i = 0; i < nd; ++i, ++it )
		fprintf( f, "%s\n", it->c_str() );

	fclose( f );




	return 0;
}


