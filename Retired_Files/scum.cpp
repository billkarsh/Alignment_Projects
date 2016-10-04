

#include	"File.h"






void MakeMakePts(FILE *fp, int N)
{
    fprintf(fp, "all: pts.down pts.same pts.up\n");
    fprintf(fp, "\n");
    fprintf(fp, "pts.same: */%d*\n", N);
    fprintf(fp, "\talign -fpts.same make.same\n");
    fprintf(fp, "\n");
    fprintf(fp, "pts.up: */%d*\n", N+1);
    fprintf(fp, "\talign -fpts.up make.up\n");
    fprintf(fp, "\n");
    fprintf(fp, "pts.down: */%d*\n",N-1);
    fprintf(fp, "\talign -fpts.down make.down\n");
}


int main( int arc, char **argv )
{
    for( int i = 161; i <= 1978; ++i ) {

        char name[128];

        sprintf( name, "temp/%d/make.pts", i );

        FILE	*fp = FileOpenOrDie( name, "w" );

        MakeMakePts( fp, i );

        fclose( fp );
    }

    return 0;
}


