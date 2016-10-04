

#include	"CAffineLens.h"
#include	"GenDefs.h"
#include	"File.h"






/* --------------------------------------------------------------- */
/* ReadFile ------------------------------------------------------ */
/* --------------------------------------------------------------- */

bool CAffineLens::ReadFile( const char *path, FILE* flog )
{
    this->flog = flog;

    FILE		*f = fopen( path, "r" );
    CLineScan	LS;
    uint8		is[4] = {0,0,0,0};
    bool		ok = false;

// Header

    if( LS.Get( f ) <= 0 ) {

        fprintf( flog,
        "Lens::Read: Empty file [%s].\n", path );
        goto exit;
    }

// Entries

    while( LS.Get( f ) > 0 ) {

        TAffine	T;
        int		cam;

        if( 7 != sscanf( LS.line,
                "%d%lf%lf%lf%lf%lf%lf",
                &cam,
                &T.t[0], &T.t[3], &T.t[1],
                &T.t[4], &T.t[2], &T.t[5] ) ) {

            fprintf( flog,
            "Lens::Read: Bad line [%s].\n", LS.line );
            goto exit;
        }

        if( cam < 0 || cam > 3 ) {

            fprintf( flog,
            "Lens::Read: Bad cam index [%s].\n", LS.line );
            goto exit;
        }

        Tf[cam] = T;
        is[cam] = 1;
    }

// Got all four?

    if( 4 != is[0] + is[1] + is[2] + is[3] ) {

        fprintf( flog,
        "Lens::Read: Missing entry [%s].\n", path );
        goto exit;
    }

// Make inverses

    for( int i = 0; i < 4; ++i ) {

//		Tf[i].SetXY( 0, 0 );
        Ti[i].InverseOf( Tf[i] );
    }

    ok = true;

exit:
    if( f )
        fclose( f );

    return ok;
}

/* --------------------------------------------------------------- */
/* ReadIDB ------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CAffineLens::ReadIDB( const string &idb, FILE* flog )
{
    char	path[2048];

    sprintf( path, "%s/lens.txt", idb.c_str() );

    return ReadFile( path, flog );
}

/* --------------------------------------------------------------- */
/* UpdateDoublesRHS ---------------------------------------------- */
/* --------------------------------------------------------------- */

void CAffineLens::UpdateDoublesRHS( double *D, int cam, bool inv )
{
    TAffine	T = TAffine( D ) * (inv ? Ti : Tf)[cam];

    T.CopyOut( D );
}

/* --------------------------------------------------------------- */
/* UpdateTFormRHS ------------------------------------------------ */
/* --------------------------------------------------------------- */

void CAffineLens::UpdateTFormRHS( TAffine &T, int cam, bool inv )
{
    T = T * (inv ? Ti : Tf)[cam];
}

/* --------------------------------------------------------------- */
/* UpdateTFormLHS ------------------------------------------------ */
/* --------------------------------------------------------------- */

void CAffineLens::UpdateTFormLHS( TAffine &T, int cam, bool inv )
{
    T = (inv ? Ti : Tf)[cam] * T;
}


