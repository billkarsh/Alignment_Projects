

#include	"TrakEM2_UTL.h"

#include	<stdio.h>
#include	<string.h>






/* --------------------------------------------------------------- */
/* NextOID ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Recursively search all children of par node for either
// 'id' or 'oid' tags.
//
// Return highest found.
//
static int NextOID_IntoNode( const TiXmlNode* par )
{
	int	highest = 0;

	const TiXmlNode*	ch = par->FirstChild();

	for( ; ch; ch = par->IterateChildren( ch ) ) {

		const char*	s;
		int			id;

		s = ((TiXmlElement*)ch)->Attribute( "oid" );

		if( s ) {

			id = atoi( s );

			if( id > highest )
				highest = id;
		}

		s = ((TiXmlElement*)ch)->Attribute( "id" );

		if( s ) {

			id = atoi( s );

			if( id > highest )
				highest = id;
		}

		id = NextOID_IntoNode( ch );

		if( id > highest )
			highest = id;
	}

	return highest;
}


// Recursively search all children of "trakem2" node for either
// 'id' or 'oid' tags.
//
// Return highest found + 1.
//
int NextOID( const TiXmlHandle hDoc )
{
	int	highest = 0;

	TiXmlNode*	par = hDoc.FirstChild( "trakem2" ).ToNode();

	if( par )
		highest = NextOID_IntoNode( par );

	return highest + 1;
}

/* --------------------------------------------------------------- */
/* SetOID -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// If par object has an 'oid' set it to nextoid and recursively
// do the same for all children of par.
//
// Return highest applied oid + 1.
//
int SetOID( TiXmlNode* par, int nextoid )
{
	const char*	s;

	s = ((TiXmlElement*)par)->Attribute( "oid" );

	if( s )
		((TiXmlElement*)par)->SetAttribute( "oid", nextoid++ );

	TiXmlNode*	ch = par->FirstChild();

	for( ; ch; ch = par->IterateChildren( ch ) )
		nextoid = SetOID( ch, nextoid );

	return nextoid;
}

/* --------------------------------------------------------------- */
/* CopyDTD ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// TrakEM2 requires input xml files to begin with DTD <!DOCTYPE>
// block, but TinyXML ignores DTD completely. Therefore, we copy
// those data from the original input file.
//
void CopyDTD( const char *dtdsrc, const char *bodysrc )
{
	char	line[2048];
	FILE	*fi, *fo;
	int		len;

// Open output file
	len = sprintf( line, "%s", dtdsrc );
	strcpy( line + len - 4, "_v2.xml" );

	if( fo = fopen( line, "w" ) ) {

		// Get DTD from input file
		if( fi = fopen( dtdsrc, "r" ) ) {

			while( fgets( line, sizeof(line), fi ) ) {

				fprintf( fo, line );

				if( !strncmp( line, "] >", 3 ) ) {
					fprintf( fo, "\n" );
					break;
				}
			}

			fclose( fi );
		}

		// Get all but first header line from bodysrc
		if( fi = fopen( bodysrc, "r" ) ) {

			fgets( line, sizeof(line), fi );

			while( fgets( line, sizeof(line), fi ) )
				fprintf( fo, line );

			fclose( fi );
		}

		fclose( fo );
	}

	remove( bodysrc );
}


