

#pragma once


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CLoadPoints {
	friend void* _Loader( void* ithr );
private:
	class CJob {
	public:
		int	z, SorD, x, y;
	public:
		CJob( int z, int SorD, int x, int y )
		: z(z), SorD(SorD), x(x), y(y) {};
	};
private:
	const char		*tempdir;
	vector<CJob>	vJ;
	int				njob,
					nthr;
private:
	void AppendJobs(
		int			z,
		int			SorD,
		int			xhi,
		int			yhi );
	void Remap();
public:
	void Load( const char *tempdir, bool isstack );
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */


