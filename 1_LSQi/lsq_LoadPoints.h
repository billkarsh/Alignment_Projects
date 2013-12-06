

#pragma once


#include	<stdio.h>


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CLoadPoints {
	friend void* _Gather( void* ithr );
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
	FILE			*fpnts;
	vector<CJob>	vJ;
	int				njob,
					nthr,
					wkid;
private:
	char* NameBinary( char *buf );
	bool IsBinary();
	void AppendJobs(
		int			z,
		int			SorD,
		int			xhi,
		int			yhi );
	void MakeBinary();
	void LoadBinary();
	void Remap();
public:
	void Load( const char *tempdir, int wkid );
};

/* --------------------------------------------------------------- */
/* Functions ----------------------------------------------------- */
/* --------------------------------------------------------------- */


