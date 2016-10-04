

#pragma once


#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class XArray {
public:
    int						NE;
    vector<vector<double> >	X;
public:
    void Resize( int ne );
    void Load( const char *path );
    void Save() const;
    void UpdtFS();
private:
    bool Send( int zlo, int zhi, int XorF, int toLorR );
    bool Recv( int zlo, int zhi, int XorF, int fmLorR );
public:
    bool Updt();
};


