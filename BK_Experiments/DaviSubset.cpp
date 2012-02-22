#include	"Disk.h"
#include	"File.h"
#include	"CTForm.h"


typedef struct {
// entry: TileToImage.txt
	TForm	T;
	int		tile;
	string	path;
} Til2Img;


static void ReadAll( vector<Til2Img> &t2i )
{
	char	name[2048];
	FILE	*f;

	t2i.clear();

	sprintf( name,
	"/groups/bock/home/bockd/bocklab/karsh/idb0/"
	"0/TileToImage.txt" );

	if( f = fopen( name, "r" ) ) {

		CLineScan	LS;

		if( LS.Get( f ) <= 0 )
			goto exit;

		while( LS.Get( f ) > 0 ) {

			Til2Img	E;
			char	buf[2048];

			sscanf( LS.line,
			"%d\t%lf\t%lf\t%lf"
			"\t%lf\t%lf\t%lf\t%[^\t\n]",
			&E.tile,
			&E.T.t[0], &E.T.t[1], &E.T.t[2],
			&E.T.t[3], &E.T.t[4], &E.T.t[5],
			buf );

			E.path = buf;

			t2i.push_back( E );
		}
	}

exit:
	if( f )
		fclose( f );
}


static void ColRow( int &col, int &row, const string &path )
{
	const char *c = path.c_str();
	const char *s;

	s = strstr( c, "col" );
	col = atoi( s + 3 );
	s = strstr( c, "row" );
	row = atoi( s + 3 );
}


static void Run()
{
	vector<Til2Img>	t2i;

	ReadAll( t2i );

	FILE	*f = FileOpenOrDie( "TileToImage.txt", "w" );
	int		nt = t2i.size();

	fprintf( f, "Tile\tT0\tT1\tX\tT3\tT4\tY\tPath\n" );

	for( int i = 0; i < nt; ++i ) {

		int	col, row;

		ColRow( col, row, t2i[i].path );

		if( col >= 38 && col <= 43 && row >= 86 && row <= 91 ) {
			fprintf( f,
				"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
				t2i[i].tile,
				t2i[i].T.t[0], t2i[i].T.t[1], t2i[i].T.t[2],
				t2i[i].T.t[3], t2i[i].T.t[4], t2i[i].T.t[5],
				t2i[i].path.c_str() );
		}
	}

	fclose( f );
}


int main( int argc, char **argv )
{
	Run();

	return 0;
}
