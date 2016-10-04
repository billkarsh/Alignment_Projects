

#include	"CRegexID.h"

#include	<string.h>






/* --------------------------------------------------------------- */
/* Set ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

void CRegexID::Set( const char* _pat )
{
    sprintf( pat, "%.*s", int(sizeof(pat)), _pat );
}

/* --------------------------------------------------------------- */
/* Compile ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Compile pattern into a regexp and matching format string.
//
// Each instance of 'N' in a pattern expands to one or more
// consecutive digits in the regexp--in grep speak "[0-9]+".
// The format string for scanf gets '%d' in the right place.
//
void CRegexID::Compile( FILE* flog )
{
    char	pat_expanded[2048]	= "";	// expanded pattern
    char	tmp[8]				= "a";
    char	*pf					= pat_format;

    for( char *p = pat; *p; ++p ) {

        if( *p == 'N' ) {
            strcat( pat_expanded, "[0-9]+" );
            *pf++ = '%';
            *pf++ = 'd';
        }
        else {
            tmp[0] = *p;
            strcat( pat_expanded, tmp );
            *pf++ = *p;
        }
    }

    regcomp( &pat_compiled, pat_expanded, REG_EXTENDED );
    *pf = 0;

// Log compilation

    if( flog ) {

        fprintf( flog,
        "Pattern '%s' regex '%s' format '%s'.\n\n",
        pat, pat_expanded, pat_format );
    }
}

/* --------------------------------------------------------------- */
/* Decode -------------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CRegexID::Decode( int &i, const char *s ) const
{
    regmatch_t	matches[12];

// regmatch_t.rm_so is 'regmatch start offset'

    return
    (REG_NOMATCH != regexec( &pat_compiled, s, 12, matches, 0 ))
    &&
    (1 == sscanf( s + matches[0].rm_so, pat_format, &i ));
}


