

#pragma once


#include	<regex.h>
#include	<stdio.h>


/* --------------------------------------------------------------- */
/* class CRegexID ------------------------------------------------ */
/* --------------------------------------------------------------- */

class CRegexID {

private:
    char		pat[16],
                pat_format[16];
    regex_t		pat_compiled;

public:
    CRegexID()						{pat[0] = 0; pat_format[0] = 0;};
    CRegexID( const char *_pat )	{Set( pat ); pat_format[0] = 0;};
    void Set( const char* _pat );
    void Compile( FILE* flog );
    bool Decode( int &i, const char *s ) const;
};


