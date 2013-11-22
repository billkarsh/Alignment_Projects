

#pragma once


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CLoadPoints {
	friend void* _Loader( void* ithr );
private:
	class PntPair {
	public:
		double	x1, y1,
				x2, y2;
		int		z1, d1, r1,
				z2, d2, r2;
	};
	class CJob {
	public:
		vector<PntPair>	vP;
		int				z, SorD, x, y, done;
	public:
		CJob( int z, int SorD, int x, int y )
		: z(z), SorD(SorD), x(x), y(y), done(false) {};
	};
private:
	const char		*tempdir;
	vector<CJob>	vJ;
	int				nj;
private:
	void AppendJobs(
		int			z,
		int			SorD,
		int			xhi,
		int			yhi );
public:
	void Load( const char *tempdir, bool isstack );
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */


