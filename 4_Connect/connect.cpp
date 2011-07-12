

#include	"LinEqu.h"
#include	"File.h"
#include	"CPicBase.h"

#include	<math.h>
#include	<string.h>

#include	<map>
#include	<stack>
#include	<set>
#include	<algorithm>
using namespace std;


#define	NM( pic )	pic.fname.c_str()


class Picture : public PicBase {

public:
	double	dist_from_center;	// distance from center if exists

public:
	Picture()	{dist_from_center = 0.0;};

	bool operator < (const Picture &rhs) const
		{return dist_from_center < rhs.dist_from_center;};
};


class ConnEntry {
  public:
    bool used;           // for use in connected components
    vector<int> Which;    // which pictures does it connect to
    vector<int> HowMany;  // how many points
    };

class Error {
  public:
    int which;  // which constraint
    double amt;
    bool operator<(const Error &rhs) const {return amt < rhs.amt;};
    };

class CorrPts {
    public:
    vector<Point> pa;
    vector<Point> pb;
    };


/* --------------------------------------------------------------- */
/* PrintMagnitude ------------------------------------------------ */
/* --------------------------------------------------------------- */

static void PrintMagnitude( const vector<double> &X )
{
	int	k = X.size() - 6;

	if( k >= 0 ) {

		double	mag	= sqrt( X[k]*X[k] + X[k+1]*X[k+1] );

		printf( "Final magnitude is %f = %.6e.\n", mag, mag );
	}
}


//----------------------Start stuff for pairwise testing ------------------------------------------
//tells if two vectors of points can be aligned by a single transform of the form
//   a x - b y + c  = x'
//   b x + a y + d  = y'   (4 free parameters, a, b, c, d)
double CanAlign(vector<Point> &p1, vector<Point> &p2, bool print)
{

// create the system of normal equations
vector<double> RHS(4, 0.0);   // create a vector for the right hand side
vector<LHSCol> LHS(4);        // a vector of columns
for(int i=0; i<p1.size(); i++) {
    int indices[3] = {0, 1, 2};
    double vals[3] = {p1[i].x,  -p1[i].y, 1.0};
    AddConstraint(LHS, RHS, 3, indices, vals, p2[i].x);
    int i2[6] =      {0, 1, 3};
    double v2[3]   = {p1[i].y,   p1[i].x, 1.0};
    AddConstraint(LHS, RHS, 3, i2     ,v2, p2[i].y);   // 0.0 is the right hand side of the constraint
    }
vector<double> X(4);
WriteSolveRead( X, LHS, RHS, false );
PrintMagnitude( X );
if( print ) {
    double mag = sqrt(X[0]*X[0] + X[1]*X[1]);
    printf(" a=%f b=%f x0=%f y0=%f mag=%f\n", X[0], X[1], X[2], X[3], mag);
    }

// find the error
double sum=0.0;
for(int i=0; i<p1.size(); i++) {
    double px = p1[i].x*X[0] - p1[i].y*X[1] + X[2];
    double py = p1[i].x*X[1] + p1[i].y*X[0] + X[3];
    double dx = px - p2[i].x;
    double dy = py - p2[i].y;
    sum += dx*dx + dy*dy;
    if( print ) {
        printf("(%8.2f %8.2f) -> (%8.2f %8.2f) =? (%8.2f %8.2f) d=%8.2f\n",
	p1[i].x, p1[i].y, px, py, p2[i].x, p2[i].y, sqrt(dx*dx + dy*dy) );
	}
    }
double rms = sqrt(sum/p1.size());
if( print ) {
    printf(" RMS error is %f pixels\n", rms);
    }
return rms;
}

map<string, int> PairIndex;  // maps the name to an index in next table
vector<CorrPts> AllPts;      // keeps track of corresponding points for each pair

void AddPairwisePoints(const char *name1, double x1, double y1, const char *name2, double x2, double y2)
{
bool InOrder = strcmp(name1, name2) < 0;
string name;  // canonical name, first (alphabetically) name, then second
if( InOrder )
    name = string(name1) + string(",") + string(name2);
else
    name = string(name2) + string(",") + string(name1);

// Now look up the string in the table, adding if necessary
int loc; // will be the location in the table when this is done
map<string, int>::iterator mi = PairIndex.find(name);
if( mi == PairIndex.end() ) {
    // did not exist; add it
    AllPts.push_back(CorrPts());
    loc = AllPts.size() - 1;
    PairIndex[name] = loc;
    }
else
    loc = mi->second;

// Now add the points in the correct order
if( InOrder ) {
    AllPts[loc].pa.push_back(Point(x1,y1));
    AllPts[loc].pb.push_back(Point(x2,y2));
    }
else {
    AllPts[loc].pa.push_back(Point(x2,y2));
    AllPts[loc].pb.push_back(Point(x1,y1));
    }
}

//------------------------------
// Given name1, name2, tells whether to ignore this line
bool InIgnoreNames(set<string> &ign, set<string> &BadPairs, const char *name1, const char *name2)
{
set<string>::iterator si = BadPairs.find(string(name1) + string(",") + string(name2));
if( si != BadPairs.end() ) {
    //printf("Reject %s,%s\n", name1, name2);
    return true;
    }
si = ign.find(string(name1));
if( si != ign.end() )
     return true;
si = ign.find(string(name2));
return si != ign.end();
}

// Converts a name to a zero based index.  Looks for the string in the mapping, and returns that
// if found.  If not, adds to the table, assigning the value nt, then increments nt.
int LookupAddIfNeeded(const char *name, map<string, int> &m, int &nt)
{
map<string,int>::iterator it = m.find(string(name));
if( it != m.end() )
    return it->second;
// else it was not there, need to add it
string *s = new string; *s = name;
m[*s] = nt;
return nt++;
}

vector<ConnEntry> ConnTable;
// keeps track of how many correspondences between pictures, to compute connected components
void AddCorrespondence(int n1, int n2)
{
if( n1 >= ConnTable.size()  )
    ConnTable.push_back(ConnEntry());
if( n2 >= ConnTable.size()  )
    ConnTable.push_back(ConnEntry());
// Now both should be within the table.  double check
if( n1 < 0 || n1 >= ConnTable.size() || n2 < 0 || n2 >= ConnTable.size()  ) {
    printf("Bad table lookup! %d %d size %d\n", n1, n2, ConnTable.size() );
    exit( 42 );
    }

// Now see if this connection is already known.  If so, increment the count.  If not,
// add the connection to the table.
int j;
for(j=0; j<ConnTable[n1].Which.size(); j++)
    if( ConnTable[n1].Which[j] == n2 )
	break;
if( j >= ConnTable[n1].Which.size() ) {
    ConnTable[n1].Which.push_back(n2);
    ConnTable[n1].HowMany.push_back(0);
    }
ConnTable[n1].HowMany[j]++;

//Now the same for n2
for(j=0; j<ConnTable[n2].Which.size(); j++)
    if( ConnTable[n2].Which[j] == n1 )
	break;
if( j >= ConnTable[n2].Which.size() ) {
    ConnTable[n2].Which.push_back(n1);
    ConnTable[n2].HowMany.push_back(0);
    }
ConnTable[n2].HowMany[j]++;
}

// Lookup by number.  Very inefficient; could (should) keep tables both ways
string LookupByNumber(map<string,int> &mapping, int k)
{
map<string,int>::iterator it;
for(it=mapping.begin(); it != mapping.end(); it++) {
    if( it->second == k )
	return it->first;
    }
printf("looking for %d, not there?\n", k);
exit( 42 );
return string("oops");
}

// Find the layer number by searching through the substrings.  Inefficient, but should be
// OK here
int FindZ(const char *name, vector<char *> &dnames, vector<int> &lnums)
{
for(int i=0; i<dnames.size(); i++) {
    if( strstr(name, dnames[i]) != NULL )
	return lnums[i];
    }
return 0;  // not found, return layer 0
}

// returns the number part after the last '/' in the name
int GetImageNumber(const char *name)
{
char *p = strrchr(name,'/');  // find the last '/'
if( p == NULL ) {
    printf("Odd - no slash in '%s'\n", name);
    exit( 42 );
    }
int n;
if( sscanf(p+1,"%d", &n) != 1 ) {
    printf("Odd - no number after last slash in '%s'\n", name);
    exit( 42 );
    }
return n;
}

//Reads a file produced by lsq.  only 4 kinds of LINES
//  TRANSFORM - transform for a given image
//  IMAGESIZE - size of each image (all are the same)
//  Compiles a list of images, sorted by their distance from the center
int ReadAFile(char *name, vector<Picture> &ps)
{
double xmin, xmax, ymin, ymax;  // BBOX of layer
double w, h;                    // W and H of images

FILE		*fdir = FileOpenOrDie( name, "r" );
CLineScan	LS;

for(;;) {
	if( LS.Get( fdir ) <= 0 )
		break;
   if( strncmp(LS.line,"TRANSFORM",9) == 0 ) {
        char fname[10240];
        double a,b,c,d,e,f;
        if( sscanf(LS.line+9,"%s %lf %lf %lf %lf %lf %lf", fname, &a, &b, &c, &d, &e, &f) != 7 ) {
	    printf("Not expecting this IMAGESIZE: '%s'\n", LS.line);
	    exit( 42 );
	    }
        Picture p;
        p.fname = strtok(fname," ':\n");
        p.tr = TForm(a,b,c,d,e,f);
        ps.push_back(p);
	}
   if( strncmp(LS.line,"IMAGESIZE",9) == 0 ) {
        if( sscanf(LS.line+9,"%lf %lf", &w, &h) != 2 ) {
	    printf("Not expecting this IMAGESIZE: '%s'\n", LS.line);
	    exit( 42 );
	    }
	}
   if( strncmp(LS.line,"BBOX",4) == 0 ) {
        if( sscanf(LS.line+4,"%lf %lf %lf %lf", &xmin, &ymin, &xmax, &ymax) != 4 ) {
	    printf("Not expecting this BBOX: '%s'\n", LS.line);
	    exit( 42 );
	    }
	}
}

printf("read %d transform statements\n", ps.size());
fclose(fdir);

// Compute the distance from the center for each image read.
for(int i=0; i<ps.size(); i++) {
    Point p(w/2, h/2);  // midpoint of the picture
    ps[i].tr.Transform( p );
    double dx = p.x - (xmin+xmax)/2.0;
    double dy = p.y - (ymin+ymax)/2.0;
    ps[i].dist_from_center = sqrt(dx*dx + dy*dy);  // sqrt not strictly needed here
    }
// Sort by distance from the center
sort(ps.begin(), ps.end());

// write in order of distance
//for(int i=0; i<ps.size(); i++)
    //printf("%12.2f %s\n", ps[i].dist_from_center, ps[i].name);
}

// Connects two layers

int main(int argc, char **argv)
{
if( argc < 4 ) {
    printf("Usage:  connect <file-1-from-lsq> <file-2-from-lsq> <ldir>\n");
    return 42;
    }

FILE		*FP = FileOpenOrDie( argv[1], "r" );
CLineScan	LS;

// How many close to the center images should we try.  Time will be N^2
int HowMany = 4;  // How many to try
for(int i=0; i<argc; i++) {
    if( strncmp(argv[i],"-N=",3) == 0 ) {
	HowMany = atoi(argv[i]+3);
	printf("Will try %d central images\n", HowMany);
	}
    }

vector<Picture> l1, l2;

ReadAFile(argv[1], l1);
ReadAFile(argv[2], l2);

int w=4056, h=4056; // size of the images (all assumed to be the same)

// (d) read the layer directory to map layer numbers and directory names, if specified
vector<int> lnums;    // layer numbers
vector<char*> dnames; // directory names where these layers are found
if( argc > 3 ) {

	FILE	*fdir = FileOpenOrDie( argv[3], "r" );

	for(;;) {

		if( LS.Get( fdir ) <= 0 )
			break;
		if( strncmp(LS.line,"DIR",3) == 0 ) {      // simply store directory string
			char *anum = strtok(LS.line+3," ");  // and layer number in parallel
			char *dname = strtok(NULL," \n");     // vectors
			dnames.push_back(strdup(dname));
			lnums.push_back(atoi(anum));
			//printf("got %s = %s\n", anum, dname);
		}
		else {
			printf("Bad line '%s' in dir file\n", LS.line);
			exit( 42 );
		}
	}

	fclose(fdir);
	printf("Read %d name:layer-number pairs.\n", dnames.size());
}

if( l1.size() < 1 || l2.size() < 1 ) {
    printf("Need at least one image on each layer - got %d %d\n", l1.size(), l2.size());
    return 42;
    }

// Find the layers involved.  We want the highest of the first set, lowest of the second.
printf("%s\n%s\n", NM(l1[0]), NM(l2[0]));
int z1 = -1;
for(int i=0; i<l1.size(); i++)
    z1 = max(z1, FindZ( NM(l1[i]), dnames, lnums ));
int z2 = 0x7FFFFFFF;
for(int i=0; i<l2.size(); i++)
    z2 = min(z2, FindZ( NM(l2[i]), dnames, lnums ));
printf("z1 = %d, z2 = %d\n", z1, z2);
if( z1 + 1 != z2 ) {
    printf("layers should be adjacent!\n");
    return 42;
    }

// Keep only the highest of the first set, lowest of the second set.
// This will result in a memory leak, but that's not important here
int k=0;
for(int i=0; i<l1.size(); i++) {
    if( FindZ( NM(l1[i]), dnames, lnums ) == z1 )
	l1[k++]  = l1[i];
    }
l1.erase(l1.begin()+k, l1.end());

k=0;
for(int i=0; i<l2.size(); i++) {
    if( FindZ( NM(l2[i]), dnames, lnums) == z2 )
	l2[k++]  = l2[i];
    }
l2.erase(l2.begin()+k, l2.end());
printf("After layer pruning: %d bottom, %d top.\n", l1.size(), l2.size());

// start make files make.up for the lower layer, and make.down for the upper.
char fname[1024];
sprintf( fname, "%d/make.up", z1 );
FILE	*fp_up = FileOpenOrDie( fname, "w" );

sprintf( fname, "%d/make.down", z2 );
FILE	*fp_down = FileOpenOrDie( fname, "w" );

// OK, write the target statements
fprintf( fp_up, "all: " );
// Now write the targets
for(int i=0; i<HowMany && i < l1.size(); i++) {
    int num1 = GetImageNumber( NM(l1[i]) );
    for(int j=0; j<HowMany && j < l2.size(); j++) {
        int num2 = GetImageNumber( NM(l2[j]) );
	fprintf(fp_up, "%d/%d.%d.map.tif ", num1, z2, num2);
	}
    }
fprintf(fp_up, "\n\n");

// Now write the make lines
for(int i=0; i<HowMany && i < l1.size(); i++) {
    int num1 = GetImageNumber( NM(l1[i]) );
    for(int j=0; j<HowMany && j < l2.size(); j++) {
        int num2 = GetImageNumber( NM(l2[j]) );
	fprintf(fp_up, "%d/%d.%d.map.tif: %s %s\n", num1, z2, num2, NM(l2[j]), NM(l1[i]) );
	fprintf(fp_up, "	ptest '%s' '%s' ../%d/%d/fm.tif %d/fm.tif ../%d/%d/%d.%d.tf.txt ../%d/%d/%d.%d.map.tif ${EXTRA}\n",
	 NM(l2[j]), NM(l1[i]), z2, num2, num1,
         z1, num1, z2, num2,
         z1, num1, z2, num2);
        fprintf(fp_up, "\n");
	}
    }
fclose(fp_up);

// OK, same for the make.down from the upper layer
// OK, write the target statements
fprintf(fp_down,"all: ");
// Now write the targets
for(int i=0; i<HowMany && i < l2.size(); i++) {
    int num1 = GetImageNumber( NM(l2[i]) );
    for(int j=0; j<HowMany && j < l1.size(); j++) {
        int num2 = GetImageNumber( NM(l1[j]) );
	fprintf(fp_down, "%d/%d.%d.map.tif ", num1, z1, num2);
	}
    }
fprintf(fp_down, "\n\n");

// Now write the make lines
for(int i=0; i<HowMany && i < l2.size(); i++) {
    int num1 = GetImageNumber( NM(l2[i]) );
    for(int j=0; j<HowMany && j < l1.size(); j++) {
        int num2 = GetImageNumber( NM(l1[j]) );
	fprintf(fp_down, "%d/%d.%d.map.tif: %s %s\n", num1, z1, num2, NM(l1[j]), NM(l2[i]) );
	fprintf(fp_down, "	ptest '%s' '%s' ../%d/%d/fm.tif %d/fm.tif ../%d/%d/%d.%d.tf.txt ../%d/%d/%d.%d.map.tif ${EXTRA}\n",
	 NM(l1[j]), NM(l2[i]), z1, num2, num1,
         z2, num1, z1, num2,
         z2, num1, z1, num2);
        fprintf(fp_down, "\n");
	}
    }
fclose(fp_down);

exit( 42 );

FILE *FOUT;
// First pass through the file, to (a) determine how many entries.
// (b), create a graph for connected components
// (c) look for bad pairs (not enough points, inconsistent points)
int nt = 0;
int nc = 0;  // number of correspondences
map<string,int> mapping;   // maps the name to a zero based index
char name1[10240], name2[10240];  // make big enough to avoid obvious problems

for(;;) {
	if( LS.Get( FP ) <= 0 )
		break;
    if( strncmp(LS.line,"POINT",5) == 0 ) {
        double x1, y1, x2, y2;
        if( sscanf(LS.line+5, "%s %lf %lf %s %lf %lf", name1, &x1, &y1, &name2, &x2, &y2) != 6 )
	    break;
        int n1 = LookupAddIfNeeded(name1, mapping, nt);
        int n2 = LookupAddIfNeeded(name2, mapping, nt);
        nc++;
        AddCorrespondence(n1, n2);
        AddPairwisePoints(name1, x1, y1, name2, x2, y2);
	}
}

// Look for suspiciouos connections. Save in a table of bad pairs
set<string> BadPairs;
for(int i=0; i<ConnTable.size(); i++) {
    for(int k=0; k<ConnTable[i].HowMany.size(); k++) {
	if( ConnTable[i].HowMany[k] < 3 ) {
	    printf("Odd... only one connection %d %d\n", i, ConnTable[i].Which[k]);
            BadPairs.insert(LookupByNumber(mapping, i) + string(",") +
	     LookupByNumber(mapping, ConnTable[i].Which[k]) );
            }
	}
    }
// Likewise, if three or more corresponding points, see if they are all consistent
map<string,int>::iterator pi;  // iterator over all pairs
for(pi=PairIndex.begin(); pi != PairIndex.end(); pi++) {
    if( AllPts[pi->second].pa.size() >= 3 ) {
        double RMS = CanAlign(AllPts[pi->second].pb, AllPts[pi->second].pa, false);
        if( RMS > 50.0 ) { // rough limit
            printf("Testing pair: '%s'\n", (pi->first).c_str());
            CanAlign(AllPts[pi->second].pa, AllPts[pi->second].pb, true); // this time print
	    }
	}
    }

// Find the connected components
vector<int> ignore;
for(int i=0; i<ConnTable.size(); i++)
    ConnTable[i].used = false;
for(int i=0; i<ConnTable.size(); i++) {
    if( ConnTable[i].used == false ) {
        int got = 0;
	printf("Starting connected components with entry %d\n", i);
        stack<int> s;
        s.push(i);
        while (!s.empty()) {
	    int j = s.top(); s.pop();  // remove top element
            if( !ConnTable[j].used ) {
		ConnTable[j].used = true;
                got++;
                if (i != 0) // not part of first connected component
		    ignore.push_back(j);
	        for(int k=0; k<ConnTable[j].Which.size(); k++)
                    if (ConnTable[j].HowMany[k] >= 3)  // only use images connected by at least 3 pts
		        s.push(ConnTable[j].Which[k]);
		}
	    }
        printf("Got %d components\n", got);
	}
    }
set<string> ignore_names;
for(int i=0; i<ignore.size(); i++) {
    printf("Ignoring %d\n", ignore[i]);
    ignore_names.insert(LookupByNumber(mapping, ignore[i]));
    }
set<string>::iterator is;
for(is=ignore_names.begin(); is != ignore_names.end(); is++)
    printf("%s\n", (*is).c_str());

// Restart, now we know which images to ignore
rewind(FP);
mapping.clear();
nt = 0;
nc = 0;

for(;;) {
	if( LS.Get( FP ) <= 0 )
		break;
    if( strncmp(LS.line,"POINT",5) == 0 ) {
        double x1, y1, x2, y2;
        if( sscanf(LS.line+5, "%s %lf %lf %s %lf %lf", name1, &x1, &y1, &name2, &x2, &y2) != 6 )
	    break;
        if(InIgnoreNames(ignore_names, BadPairs, name1, name2))
	    continue;
        int n1 = LookupAddIfNeeded(name1, mapping, nt);
        int n2 = LookupAddIfNeeded(name2, mapping, nt);
        nc++;
	}
    else if( strncmp(LS.line,"IMAGESIZE",9) == 0 ) {
        if( sscanf(LS.line+9, "%d %d", &w, &h) != 2 ) {
	    printf("Bad image size line '%s'\n", LS.line);
	    return 42;
	    }
	}
}

printf("Image size is w=%d h=%d\n", w, h);

printf("%d transforms to be determined, %d point correspondences\n", nt, nc);
int nvars = 6*nt;  // 6 variable for each transform
int neqs  = 2*nc + 6 + 2*nt;  // each point makes 2 constraints (x and y)
                              // plus 6 to define the fixed tile
                              // plus 2 for regularization for each transform

printf("%d equations in %d unknowns\n", neqs, nvars);
// create the system of normal equations
vector<double> RHS(nvars, 0.0);   // create a vector for the right hand side
vector<LHSCol> LHS(nvars);       // a vector of columns

// write directory out first.  Write the directory to layer number map
for(int i=0; i<dnames.size(); i++)
    fprintf(FOUT, "DIR %d %s\n", lnums[i], dnames[i] );

double scale = 10000.0;
rewind(FP);
vector<int> z(nt, -1);

for(;;) {
	if( LS.Get( FP ) <= 0 )
		break;
    if( strncmp(LS.line,"POINT",5) == 0 ) {
        double x1, y1, x2, y2;
        if( sscanf(LS.line+5, "%s %lf %lf %s %lf %lf", name1, &x1, &y1, name2, &x2, &y2) != 6 )
	    break;
        if(InIgnoreNames(ignore_names, BadPairs, name1, name2))
	    continue;
        int t1 = LookupAddIfNeeded(name1, mapping, nt);
        if( z[t1] == -1 )
	    z[t1] = FindZ(name1, dnames, lnums);
        int t2 = LookupAddIfNeeded(name2, mapping, nt);
        if( z[t2] == -1 )
	    z[t2] = FindZ(name2, dnames, lnums);
        x1 /= scale; y1 /= scale; x2 /= scale; y2 /= scale;
        int j = t1*6;  // variable number
        int k = t2*6;
        int indices[6] = {j, j+1, j+2, k, k+1,  k+2};
        double vals[6] = { x1,  y1, 1.0, -x2, -y2, -1.0};
        AddConstraint(LHS, RHS, 6, indices, vals, 0.0);   // 0.0 is the right hand side of the constraint
        int i2[6] =      {j+3, j+4, j+5, k+3, k+4,  k+5};
        AddConstraint(LHS, RHS, 6, i2     , vals, 0.0);   // 0.0 is the right hand side of the constraint
        }
    else if( strncmp(LS.line,"FOLDMAP",7) == 0 )
	fprintf(FOUT, "%s", LS.line);
    else if (strncmp(LS.line,"IMAGESIZE",9) != 0){
	printf("Unknown line type '%s'", LS.line);
	return 42;
	}
}

// OK, add constraints for the initial square
double InitStiff = 10000.0/scale;
int n;
n = 0;  AddConstraint(LHS, RHS, 1, &n, &InitStiff, InitStiff);
n = 1;  AddConstraint(LHS, RHS, 1, &n, &InitStiff,       0.0);
n = 2;  AddConstraint(LHS, RHS, 1, &n, &InitStiff,       0.0);
n = 3;  AddConstraint(LHS, RHS, 1, &n, &InitStiff,       0.0);
n = 4;  AddConstraint(LHS, RHS, 1, &n, &InitStiff, InitStiff);
n = 5;  AddConstraint(LHS, RHS, 1, &n, &InitStiff,       0.0);

// Add the constraints so the transforms stay close to square
double stiff_sim = 10000.0/scale;
for(int i=0; i<nt; i++) {
    int j = i*6;
    int i3[2] = {j, j+4};  // cosine terms
    double vs[2] = {stiff_sim, -stiff_sim};
    AddConstraint(LHS, RHS, 2, i3, vs, 0.0);  // cosine terms should be equal
    i3[0] = j+1; i3[1] = j+3; // sin() terms
    vs[1] = -vs[1];
    AddConstraint(LHS, RHS, 2, i3, vs, 0.0);  // sin() terms should be opposite
    }

vector<double> X(nvars);
WriteSolveRead( X, LHS, RHS, false );
PrintMagnitude( X );

// add some constraints so the left edge of the array needs to be the same as the right edge,
// and the same for top and bottom. Directions are for conventional xy coordinates.
const int MAXINT = 0x7FFFFFFF;
int lownum = MAXINT;
for(int i=0; i<nt; i++)
    lownum = min(lownum, z[i]);
printf("Lowest layer is %d\n",lownum);

int nwbest, nebest, swbest, sebest;
double nw = -MAXINT, ne = -MAXINT, se = -MAXINT, sw = -MAXINT;
for(int i=0; i<nt; i++) {
    if( z[i] != lownum )
	continue;
    int j = i*6;
    double t;
    t =  X[j+2] + X[j+5];
    if( t > ne ) { nebest = j; ne = t;}
    t =  X[j+2] - X[j+5];
    if( t > se ) { sebest = j; se = t;}
    t = -X[j+2] + X[j+5];
    if( t > nw ) { nwbest = j; nw = t;}
    t = -X[j+2] - X[j+5];
    if( t > sw ) { swbest = j; sw = t;}
    }
printf("Corners are: se %d (%f %f), ne %d (%f %f), nw %d (%f %f), sw %d (%f %f)\n",
    sebest, X[sebest+2], X[sebest+5],   nebest, X[nebest+2], X[nebest+5],
    nwbest, X[nwbest+2], X[nwbest+5],   swbest, X[swbest+2], X[swbest+5] );

// add the constraints
{
double oomph = 10000;
int indices[4] = {sebest+2, swbest+2, nebest+2, nwbest+2};
double   va[4] = {oomph,    -oomph,   -oomph,   oomph};
AddConstraint(LHS, RHS, 4, indices, va, 0.0);
int indice2[4] = {sebest+5, swbest+5, nebest+5, nwbest+5};
AddConstraint(LHS, RHS, 4, indice2, va, 0.0);
}

WriteSolveRead( X, LHS, RHS, false );
PrintMagnitude( X );

// Now add the new constraints to preserve the angle, but push the magnitude to 1
double scale_stiff = 1E4/scale;
for(int i=0; i<nt; i++) {
    int j = i*6;
    double a = X[j];
    double b = X[j+1];
    double m = sqrt(a*a+b*b);
    int ind[2] = {j,j+1};
    double va[2] = {a/m*scale_stiff, b/m*scale_stiff};
    AddConstraint(LHS, RHS, 2, ind, va, scale_stiff);
    }
WriteSolveRead( X, LHS, RHS, false );
PrintMagnitude( X );

//rescale the constant terms back to pixel space
map<string,int>::iterator it;
for(it=mapping.begin(); it != mapping.end(); it++) {
    int j = it->second * 6;
    X[j+2] *= scale;
    X[j+5] *= scale;
    }

// Now find the bounding box in global space
const double BIG = 1.0e30;
double xmin = BIG, ymin = BIG;
double xmax = -BIG, ymax = -BIG;
vector<Point> pts;
pts.push_back(Point(0.0,0.0));
pts.push_back(Point(w-1,0.0));
pts.push_back(Point(w-1, h-1));
pts.push_back(Point(0.0, h-1));
for(int i=0; i<nt; i++) {
    int j=i*6;
    for(int k=0; k<4; k++) {
        double x = X[j  ]*pts[k].x + X[j+1]*pts[k].y + X[j+2];
        double y = X[j+3]*pts[k].x + X[j+4]*pts[k].y + X[j+5];
	xmin = fmin(xmin, x);
	ymin = fmin(ymin, y);
	xmax = fmax(xmax, x);
	ymax = fmax(ymax, y);
	}
    }
printf("Bounds of global image are x=[%f %f] y=[%f %f]\n", xmin, xmax, ymin, ymax);

// Modify all transformations so the min is (0,0).  Always change by an integer amount, to make later
// layer-to-layer alignment easier
int xfl = int(floor(xmin));
int yfl = int(floor(ymin));
xmin -= xfl; xmax -= xfl;              // change the bounding box
ymin -= yfl; ymax -= yfl;
for(int i=0; i<nt; i++) {
    int j=i*6;
    X[j+2] -= xfl;
    X[j+5] -= yfl;
    }
printf("Bounds of global image are x=[%f %f] y=[%f %f]\n", xmin, xmax, ymin, ymax);

// Write the transforms
double smin = 100.0;
double smax = 0.0;
double smag = 0.0;
for(it=mapping.begin(); it != mapping.end(); it++) {
    int j = it->second * 6;
    fprintf(FOUT, "TRANSFORM %s %f %f %f %f %f %f\n", (it->first).c_str(), X[j], X[j+1], X[j+2], X[j+3], X[j+4], X[j+5]);
    double det = X[j]*X[j+4]-X[j+1]*X[j+3];
    //printf(" det= %f\n", det);
    smag += sqrt(det);
    smin = fmin(smin, sqrt(det));
    smax = fmax(smax, sqrt(det));
    }
printf("Average magnitude %f, min=%f max=%f ratio=%f\n", smag/mapping.size(), smin, smax, smax/smin );

// Find the residuals
rewind(FP);
vector<Error> Errors;
double sum = 0.0;
double biggest = 0.0;
int nn = 0;

for(;;) {
	if( LS.Get( FP ) <= 0 )
		break;
    if( strncmp(LS.line,"POINT",5) == 0 ) {
        double x1, y1, x2, y2;
        if( sscanf(LS.line+5, "%s %lf %lf %s %lf %lf", name1, &x1, &y1, name2, &x2, &y2) != 6 )
	    break;
        if(InIgnoreNames(ignore_names, BadPairs, name1, name2))
	    continue;
        int t1 = LookupAddIfNeeded(name1, mapping, nt);
        int t2 = LookupAddIfNeeded(name2, mapping, nt);
        int j = t1*6;  // variable number
        int k = t2*6;
        double gx1 = x1*X[j  ] + y1*X[j+1] + X[j+2];
        double gx2 = x2*X[k  ] + y2*X[k+1] + X[k+2];
        double gy1 = x1*X[j+3] + y1*X[j+4] + X[j+5];
        double gy2 = x2*X[k+3] + y2*X[k+4] + X[k+5];
        fprintf(FOUT, "MPOINTS %d %f %f %d %f %f\n", z[t1], gx1, gy1, z[t2], gx2, gy2);
        double dx = gx1 - gx2;
        double dy = gy1 - gy2;
        //printf("%s %f %f\n", name1, gx1, gy1);
        //printf("%s %f %f\n", name2, gx2, gy2);
        //printf("dx %f, dy %f\n", dx, dy);
        double err = dx*dx + dy*dy;
        sum += err;
        Error e; e.which = nn; e.amt = err;
        Errors.push_back(e);
        biggest = max(biggest, err);
        nn++;
        }
}

const char *flag = (sqrt(sum/nn) > 20.0 || sqrt(biggest) > 75.0) ? "<-------------------" : "";
printf("RMS error %8.2f, max error %8.2f %s\n", sqrt(sum/nn), sqrt(biggest), flag);
sort(Errors.begin(), Errors.end());
for(int i=max(0,int(Errors.size())-10); i<Errors.size(); i++)
    printf("%d: Constraint %d, err %f\n", i, Errors[i].which, sqrt(Errors[i].amt) );

// Now print out the 10 biggest errors
double big = Errors.size() > 10 ? Errors[Errors.size()-10].amt : 0.0;
rewind(FP);
nn = 0;

for(;;) {
	if( LS.Get( FP ) <= 0 )
		break;
    if( strncmp(LS.line,"POINT",5) == 0 ) {
        double x1, y1, x2, y2;
        if( sscanf(LS.line+5, "%s %lf %lf %s %lf %lf", name1, &x1, &y1, name2, &x2, &y2) != 6 )
	    break;
        if(InIgnoreNames(ignore_names, BadPairs, name1, name2))
	    continue;
        int t1 = LookupAddIfNeeded(name1, mapping, nt);
        int t2 = LookupAddIfNeeded(name2, mapping, nt);
        int j = t1*6;  // variable number
        int k = t2*6;
        double gx1 = x1*X[j  ] + y1*X[j+1] + X[j+2];
        double gx2 = x2*X[k  ] + y2*X[k+1] + X[k+2];
        double gy1 = x1*X[j+3] + y1*X[j+4] + X[j+5];
        double gy2 = x2*X[k+3] + y2*X[k+4] + X[k+5];
        double dx = gx1 - gx2;
        double dy = gy1 - gy2;
        double err = dx*dx + dy*dy;
        if( err >= big ) {
	    printf("'%s' '%s' %f\n", name1, name2, sqrt(err));
	    }
        nn++;
        }
}

fclose(FP);

return 0;
}
