

#include	"PipeFiles.h"
#include	"File.h"






static void Run()
{
    vector<Til2Img>	t2i;

    IDBT2IGetAll( t2i,
    string( "/groups/bock/home/bockd/bocklab/karsh/idb0" ), 0 );

    FILE	*f = FileOpenOrDie( "TileToImage.txt", "w" );
    int		nt = t2i.size();

    fprintf( f, "ID\tT0\tT1\tX\tT3\tT4\tY\tCol\tRow\tCam\tPath\n" );

    for( int i = 0; i < nt; ++i ) {

        const Til2Img	&I = t2i[i];

        if( I.col >= 38 && I.col <= 43 &&
            I.row >= 86 && I.row <= 91 ) {

            fprintf( f,
                "%d"
                "\t%f\t%f\t%f\t%f\t%f\t%f"
                "\t%d\t%d\t%d"
                "\t%s\n",
                I.id,
                I.T.t[0], I.T.t[1], I.T.t[2],
                I.T.t[3], I.T.t[4], I.T.t[5],
                I.col, I.row, I.cam,
                I.path.c_str() );
        }
    }

    fclose( f );
}


int main( int argc, char **argv )
{
    Run();

    return 0;
}


