

#include	"CCorrImages.h"

#include	"tinyxml.h"


// ---------------------------------------------------------------
// !! Compilation requires:
//	- g++ switch -DTIXML_USE_STL
//	- modules: {tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp}
// ---------------------------------------------------------------






/* --------------------------------------------------------------- */
/* CCorrImages::Read --------------------------------------------- */
/* --------------------------------------------------------------- */

CCorrImages* CCorrImages::Read( const char *fname, FILE *flog )
{
	CCorrImages*	CI = new CCorrImages;
	TiXmlDocument	doc(fname);
	bool			loadOK = doc.LoadFile();

	printf( "CorrImage: XML load of corr file gives %d\n", loadOK );

	if( !loadOK ) {
		printf(
		"CorrImage: Could not open XML file [%s].\n", fname );
		fprintf( flog,
		"CorrImage: Could not open XML file [%s].\n", fname );
		exit( 42 );
	}

	TiXmlHandle		hDoc(&doc);
	TiXmlElement*	child;
	TiXmlHandle		hRoot(0);

// block: should be <image_correlations>
	TiXmlNode* node = doc.FirstChild();

	if( node == NULL )
		printf( "CorrImage: No first node??\n" );

	child = hDoc.FirstChild( "image_correlations" )
				.FirstChild().ToElement();

// should always have a valid root but
// handle gracefully if it doesn't

	if( !child ) {
		printf(
		"CorrImage: File exists, but no correlations - OK...\n" );
		return CI;
	}

	for( child; child; child=child->NextSiblingElement() ) {

		const char	*what = child->Value();

		if( !strcmp( what, "map" ) ) {

			const char	*id	= child->Attribute( "id" );
			const char	*nm	= child->Attribute( "name" );
			map<string,int>::iterator	i1;

			//printf( "CorrImage: id=%s, name=[%s].\n", id, nm );

			i1 = CI->names.find( nm );	// add if not there

			if( i1 == CI->names.end() )
				CI->names.insert( pair<string,int>( nm, atoi( id ) ) );
		}
		else if( !strcmp( what, "pair" ) ) {

			const char	*id	= child->Attribute( "idpair" );
			int			id1, id2;

			sscanf( id,"%d %d", &id1, &id2 );

			CorrPair	c( id1, id2 );

			//printf( "CorrImage: Inserting pair %d %d\n", id1, id2 );

			TiXmlElement*	pa = child->FirstChildElement( "pts" );

			for( pa; pa; pa = pa->NextSiblingElement() ) {

				const char *xy1 = pa->Attribute( "xy1" );
				const char *xy2 = pa->Attribute( "xy2" );
				Point		p1, p2;

				//printf( "CorrImage: xypair : %s %s.\n", xy1, xy2 );

				sscanf( xy1, "%lf %lf", &p1.x, &p1.y );
				sscanf( xy2, "%lf %lf", &p2.x, &p2.y );

				c.p1s.push_back( p1 );
				c.p2s.push_back( p2 );
			}

			CI->corrs.insert(c);
		}
	}

	return CI;
}

/* --------------------------------------------------------------- */
/* CCorrImages::Write -------------------------------------------- */
/* --------------------------------------------------------------- */

int CCorrImages::Write( const char *fname )
{
	FILE	*fp = fopen( fname, "w" );

	if( fp == NULL ) {
		printf( "CorrImage: Could not open [%s] for write.\n", fname );
		return 0;
	}

	fprintf( fp, "<image_correlations>\n" );

// write the map
	map<string,int>::iterator	i;

	for( i = names.begin(); i != names.end(); ++i ) {

		fprintf( fp, "<map id=\"%d\" name=\"%s\" />\n",
			i->second, i->first.c_str() );
	}

// write the correlations
	set<CorrPair>::iterator	it;

	for( it = corrs.begin(); it != corrs.end(); ++it ) {

		fprintf( fp, "<pair idpair=\"%d %d\">\n", it->i1, it->i2 );

		for( int j=0; j<it->p1s.size(); ++j ) {

			fprintf( fp, " <pts xy1=\"%f %f\" xy2=\"%f %f\" />\n",
			it->p1s[j].x, it->p1s[j].y, it->p2s[j].x, it->p2s[j].y );
		}

		fprintf( fp,"</pair>\n" );
	}

	fprintf( fp, "</image_correlations>\n" );

	fclose( fp );
	return 1;
}

/* --------------------------------------------------------------- */
/* CCorrImages::Find --------------------------------------------- */
/* --------------------------------------------------------------- */

int CCorrImages::Find(
	string			name1,
	string			name2,
	vector<Point>	&p1s,
	vector<Point>	&p2s )
{
	map<string,int>::iterator	i1, i2;

	p1s.clear();
	p2s.clear();

	i1 = names.find( name1 );
	i2 = names.find( name2 );

	if( i1 == names.end() || i2 == names.end() )
		return 0;

	printf( "CorrImage: Found names...\n" );

	set<CorrPair>::iterator	s1;

	printf( "CorrImage: Looking for pair %d %d\n",
		i1->second, i2->second );

	if( i1->second < i2->second ) { // already in the right order
		s1 = corrs.find( CorrPair( i1->second, i2->second ) );

		if( s1 == corrs.end() )
			return 0;

		printf( "CorrImage: Found correlated points (same order).\n" );

		p1s = s1->p1s;
		p2s = s1->p2s;
	}
	else {
		s1 = corrs.find( CorrPair( i2->second, i1->second ) );

		if( s1 == corrs.end() )
			return 0;

		printf(
		"CorrImage: Found correlated points (reversed order).\n" );

		p1s = s1->p2s;
		p2s = s1->p1s;
	}

	return p1s.size();
}

/* --------------------------------------------------------------- */
/* CCorrImages::Add ---------------------------------------------- */
/* --------------------------------------------------------------- */

int CCorrImages::Add(
	string			name1,
	string			name2,
	vector<Point>	&p1s,
	vector<Point>	&p2s )
{
	map<string,int>::iterator	it1, it2;

	it1 = names.find( name1 );	// add if not there

	if( it1 == names.end() ) {
		int j = names.size();
		names.insert( pair<string,int>( name1, j ) );
		it1 = names.find( string( name1 ) );	// Now cannot fail
	}

	it2 = names.find(name2);

	if( it2 == names.end() ) {
		int j = names.size();
		names.insert( pair<string,int>( name2, j ) );
		it2 = names.find( string( name2 ) );	// Now cannot fail
	}

	CorrPair c( it1->second, it2->second );

	if( it1->second < it2->second ) {	// normal order
		c.p1s	= p1s;
		c.p2s	= p2s;
	}
	else {	// swap order - always want lowest index first
		c.i1	= it2->second;
		c.i2	= it1->second;
		c.p1s	= p2s;
		c.p2s	= p1s;
	}

	corrs.insert( c );
	return 1;
}


