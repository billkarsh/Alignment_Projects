

#include	"lsq_CNX.h"

#include	"File.h"

#include	<stack>
using namespace std;


/* --------------------------------------------------------------- */
/* ListWeakConnections ------------------------------------------- */
/* --------------------------------------------------------------- */

// Add to r12Bad those region pairs having fewer than minLinks
// corr. points between them.
//
void CNX::ListWeakConnections( set<CRPair> &r12Bad )
{
	if( minLinks <= 0 )
		return;

// Scan for weak pairs

	printf( "---- Scan for weakly connected pairs ----\n" );

	int	nct = cnxTbl.size();

	for( int i = 0; i < nct; ++i ) {

		const CnxEntry&	Ci = cnxTbl[i];
		int				np = Ci.nlinks.size();

		for( int j = 0; j < np; ++j ) {

			if( Ci.nlinks[j] < minLinks )
				r12Bad.insert( CRPair( i, Ci.linkto[j] ) );
		}
	}

// Print weak pairs

	if( r12Bad.size() ) {

		printf(
		"%d weak pair entries reported in 'FewCorr'.\n",
		r12Bad.size() );

		FILE	*f = FileOpenOrDie( "FewCorr", "w" );

		fprintf( f, "Lyr.til:rgn pairs with few ctl pts between\n" );

		set<CRPair>::iterator	pi;

		for( pi = r12Bad.begin(); pi != r12Bad.end(); ++pi ) {

			const RGN	&A = vRgn[pi->a];
			const RGN	&B = vRgn[pi->b];

			fprintf( f, "%d.%d:%d - %d.%d:%d\n",
			A.z, A.id, A.rgn, B.z, B.id, B.rgn );
		}

		fclose ( f );
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* MaxConnectedSet ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Find the longest interconnected set of regions, and
// fill list (ignore) with those that do not belong.
//
void CNX::MaxConnectedSet( set<int> &ignore )
{
	printf( "---- Find maximum set ----\n" );

	int	nct			= cnxTbl.size(),
		bestSeed	= -1,
		bestLen		= -1;

/* --------------------------------- */
/* Init all seeds to -1 'unassigned' */
/* --------------------------------- */

	for( int i = 0; i < nct; ++i )
		cnxTbl[i].seed = -1;

/* ------------------------------------------ */
/* Consider each entry as a possible new seed */
/* ------------------------------------------ */

	for( int i = 0; i < nct; ++i ) {

		/* ----------------- */
		/* Already assigned? */
		/* ----------------- */

		if( cnxTbl[i].seed >= 0 )
			continue;

		/* -------------- */
		/* Got a new seed */
		/* -------------- */

		printf( "Next seed %d...\n", i );

		// Push the new seed...
		// ...and, as we go, its connected
		// neighbors, for further processing.

		stack<int>	s;
		int			len = 0;

		s.push( i );

		/* --------------------------------- */
		/* Pop and process connected entries */
		/* --------------------------------- */

		while( !s.empty() ) {

			/* ------- */
			/* Get one */
			/* ------- */

			int			j   = s.top();
			CnxEntry	&Cj = cnxTbl[j];

			s.pop();

			/* ------------------ */
			/* Already processed? */
			/* ------------------ */

			if( Cj.seed >= 0 )
				continue;

			/* ----------------- */
			/* Belongs to seed i */
			/* ----------------- */

			Cj.seed = i;
			++len;

			/* --------------------------------- */
			/* Push all well-connected neighbors */
			/* --------------------------------- */

			int	nneib = Cj.linkto.size();

			for( int k = 0; k < nneib; ++k ) {

				// require minLinks corr pts
				if( Cj.nlinks[k] >= minLinks )
					s.push( Cj.linkto[k] );
			}
		}

		/* ------------------- */
		/* This seed completed */
		/* ------------------- */

		printf( "...Got %d connections.\n", len );

		if( len > bestLen ) {
			bestSeed	= i;
			bestLen		= len;
		}

		/* --------------------------------- */
		/* If > nct/2 can't get a longer one */
		/* --------------------------------- */

		if( len * 2 >= nct )
			break;
	}

	printf( "Maximum length %d of %d candidates.\n", bestLen, nct );

/* -------------------------------- */
/* Create and print the ignore list */
/* -------------------------------- */

	if( bestLen < nct ) {

		printf( "\nIgnoring:\n" );

		for( int i = 0; i < nct; ++i ) {

			if( cnxTbl[i].seed != bestSeed ) {

				const RGN	&A = vRgn[i];

				ignore.insert( i );
				printf( "\t%d.%d:%d\n", A.z, A.id, A.rgn );
			}
		}
	}

	printf( "\n" );
}

/* --------------------------------------------------------------- */
/* Set_itr_set_used ---------------------------------------------- */
/* --------------------------------------------------------------- */

// Assign each regions's global-space transform index (itr).
// The itrs will be: {-1=invalid, 0+=valid} and a region's
// k transform values start at location X[itr * k].
//
// Set the Constraint.used flags selecting those that will be fit.
//
// Return valid transform count.
//
int CNX::Set_itr_set_used( set<CRPair> &r12Bad, set<int> &ignore )
{
	printf( "---- Count valid regions ----\n" );

	int	nGoodC = 0, nTrans = 0, nc = vAllC.size();

	for( int i = 0; i < nc; ++i ) {

		int	r1 = vAllC[i].r1,
			r2 = vAllC[i].r2;

		vAllC[i].used = false;

		if( ignore.find( r1 ) != ignore.end() ||
			ignore.find( r2 ) != ignore.end() ||
			r12Bad.find( CRPair( r1, r2 ) ) != r12Bad.end() ) {

			continue;
		}

		if( vRgn[r1].itr == -1 )
			vRgn[r1].itr = nTrans++;

		if( vRgn[r2].itr == -1 )
			vRgn[r2].itr = nTrans++;

		vAllC[i].used = true;
		++nGoodC;
	}

	printf(
	"%d transforms to be determined, %d point correspondences.\n",
	nTrans, nGoodC );

	if( nTrans == 0 || nGoodC < minLinks ) {
		printf( "Too few transforms or constraints.\n" );
		exit( 42 );
	}

	printf( "\n" );

	return nTrans;
}

/* --------------------------------------------------------------- */
/* AddCorrespondence --------------------------------------------- */
/* --------------------------------------------------------------- */

// Update cnxTbl whose entries track which regions connect
// to this one and how many points compose that connection.
//
// Note cnxTbl.size() kept same as vRgn.size().
//
void CNX::AddCorrespondence( int r1, int r2 )
{
// Make room

	if( r1 >= cnxTbl.size() || r2 >= cnxTbl.size() )
		cnxTbl.resize( max( r1, r2 ) + 1 );

	CnxEntry	&C1 = cnxTbl[r1],
				&C2 = cnxTbl[r2];
	int			j, n;

// Counts for r1

	n = C1.linkto.size();

	for( j = 0; j < n; ++j ) {

		if( C1.linkto[j] == r2 ) {
			++C1.nlinks[j];
			goto do_n2;
		}
	}

	C1.linkto.push_back( r2 );
	C1.nlinks.push_back( 1 );

// Counts for r2

do_n2:
	n = C2.linkto.size();

	for( j = 0; j < n; ++j ) {

		if( C2.linkto[j] == r1 ) {
			++C2.nlinks[j];
			return;
		}
	}

	C2.linkto.push_back( r1 );
	C2.nlinks.push_back( 1 );
}

/* --------------------------------------------------------------- */
/* SelectIncludedImages ------------------------------------------ */
/* --------------------------------------------------------------- */

// Decide which regions to include based on connection strength.
//
// The set<> container holds sorted, unique key-values.
//
// Return valid transform count.
//
int CNX::SelectIncludedImages( int _minLinks )
{
	set<CRPair>	r12Bad;		// suspicious region pairs
	set<int>	ignore;		// unconnected regions

	minLinks = _minLinks;

	ListWeakConnections( r12Bad );
	MaxConnectedSet( ignore );
	return Set_itr_set_used( r12Bad, ignore );
}


