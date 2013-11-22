/* --------------------------------------------------------------- */
/* DecimateVector ------------------------------------------------ */
/* --------------------------------------------------------------- */

// This is Lou's original (cleaned up) implementation.

class Decimate {  // used for decimating Point lists
public:
	double	v;
	int		x, y;
public:
	inline bool operator < (const Decimate &rhs) const
		{return y < rhs.y || (y == rhs.y && x < rhs.x);};
};


// Bin point and value lists by linear dimension lbin.
//
void DecimateVector(
	vector<Point>	&p,
	vector<double>	&v,
	int				lbin )
{
// First, calculate a list of all (integer) coordinates that will
// appear in the output. There will be several input points that
// map to the same output point.

	int					np = p.size();
	vector<Decimate>	dv( np );

	for( int i = 0; i < np; ++i ) {

		dv[i].v	= v[i];
		dv[i].x	= int(p[i].x / lbin);
		dv[i].y	= int(p[i].y / lbin);
	}

// Sort by coordinates it so the points that
// map together are adjacent (grouped).

	sort( dv.begin(), dv.end() );

// Average the point groups.

	int	j, k = 0;

	for( int i = 0; i < np; i = j ) {

		double	sum = dv[i].v;

		// Advance j until a point goes off the end,
		// or is a different point.

		for( j = i + 1; j < np && !(dv[i] < dv[j]); ++j )
			sum += dv[j].v;

		// j is the first one that failed,
		// so (j-i) is number that succeeded.

		v[k]	= sum / (j-i);
		p[k]	= Point( dv[i].x, dv[i].y );
		++k;
	}

	v.resize( k );
	p.resize( k );
}


