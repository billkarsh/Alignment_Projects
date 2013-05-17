

#pragma once


#include	"lsq_Types.h"

#include	<set>
using namespace std;


/* --------------------------------------------------------------- */
/* CNX ----------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Connection table operations

class CNX {

private:
	int	minLinks,
		minPrsPerLink;

	typedef struct {
		vector<int>	linkto;	// which regions it connects to
		vector<int>	npairs;	// # points in that connection
		int			seed;	// starting entry of conn graph
	} CnxEntry;

private:
	vector<CnxEntry>	cnxTbl;

private:
	void ListWeakConnections( set<CRPair> &r12Bad );
	bool IsWeakLink( int irgn );
	void MaxConnectedSet( set<int> &ignore );
	int  Set_itr_set_used( set<CRPair> &r12Bad, set<int> &ignore );

public:
	void AddCorrespondence( int r1, int r2 );
	int  SelectIncludedImages( int _minLinks, int _minPrsPerLink );
};


