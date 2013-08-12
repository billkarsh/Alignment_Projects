

#include	"Disk.h"
#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"TAffine.h"
#include	"Timer.h"

#include	<string.h>


// This program (and those that follow) assume all images are the same size;
int gw = -1;
int gh;

class TiffInfo {
  public:
    char *name;
    uint32 w,h;
    uint8 *raster;
    int last_used;
    };

static bool Transpose = false;      // Use translation only?








const int CSIZE = 100;
vector<TiffInfo> TiffCache;
int Tick = 0;  // for finding the oldest

uint8 *ReadTiffWithCache(char *name, uint32 &w, uint32 &h, FILE *flog, bool Transpose)
{
Tick++;
for(int i=0; i<TiffCache.size(); i++) {
    if( strcmp(TiffCache[i].name, name) == 0 ) {
	w = TiffCache[i].w;
	h = TiffCache[i].h;
        TiffCache[i].last_used = Tick;
	return TiffCache[i].raster;
	}
    }
// Not in the cache, so read it
TiffInfo ti;
ti.raster = Raster8FromTif(name, ti.w, ti.h, flog, Transpose);
ti.name = strdup(name);
ti.last_used = Tick;

// If the cache is not full, add it
int oldest = -1;
if( TiffCache.size() <= CSIZE ) {
    oldest = TiffCache.size();
    TiffCache.push_back(ti);
    }
else { // find the last used one, and replace it
    int old = Tick;
    for(int i=0; i<TiffCache.size(); i++) {
	if( TiffCache[i].last_used < old ) {
	    old = TiffCache[i].last_used;
            oldest = i;
	    }
        }
    free(TiffCache[oldest].name);
    RasterFree(TiffCache[oldest].raster);
    TiffCache[oldest] = ti;
    }
w = TiffCache[oldest].w;
h = TiffCache[oldest].h;
return TiffCache[oldest].raster;
}





// noa[0] = index string
// noa[1] = a image
// noa[2] = b image
// noa[3] = if specified, fold mask for a
// noa[4] = if specified, fold mask for b
// noa[5] = transform file (txt)
// noa[6] = transform file (tif)
//
static void ReadAFileSet(
	vector<char*>	&noa,
	FILE			*flog,
	FILE			*fcorr,
	bool			NoFolds )
{
	if( !DskExists( noa[5] ) || !DskExists( noa[6] ) )
		return;

time_t t0 = time(NULL);
char atime[32];
strcpy(atime, ctime(&t0));
//atime[24] = '\0';  // remove the newline
fprintf(flog,"start: %s", atime);
fflush(flog);
printf("------------------------------------------------------------------------------\n");
clock_t started = StartTiming();

//read the map file
uint8 *rmap;
uint32 w, h;
rmap = ReadTiffWithCache(noa[6], w, h, flog, Transpose);
if( gw == -1 ) { // never been set
    gw = w;
    gh = h;
    fprintf(fcorr,"IMAGESIZE %d %d\n", w, h);
    }

// Read the fold masks
uint8 *fold_mask_a, *fold_mask_b;
uint32 w1, w2, h1, h2;
int npixels = w*h;
if( strcmp(noa[3],"none") != 0 && !NoFolds )
    fold_mask_a = ReadTiffWithCache(noa[3], w1, h1, flog, Transpose);
else {
    fold_mask_a = (uint8*)RasterAlloc( npixels * sizeof(uint8) );
    for(int i=0; i<npixels; i++) fold_mask_a[i] = 1;
    w1 = w; h1 = h;
    }
if( strcmp(noa[4],"none") != 0 && !NoFolds )
    fold_mask_b = ReadTiffWithCache(noa[4], w2, h2, flog, Transpose);
else {
    fold_mask_b = (uint8*)RasterAlloc( npixels * sizeof(uint8) );
    for(int i=0; i<npixels; i++) fold_mask_b[i] = 1;
    w2 = w; h2 = h;
    }
if( w1 != gw || w2 != gw || h1 != gh || h2 != gh ) {
    printf("Fold masks are different size?\n");
    fprintf(flog,"Fold masks are different size?\n");
    fclose(flog);
    exit( 42 );
    }



// Read the array of transforms
FILE	*ftxt = FileOpenOrDie( noa[5], "r", flog );

vector<double> array_of_transforms;
double v;
while (fscanf(ftxt, "%lf", &v) == 1)
    array_of_transforms.push_back(v);
fclose(ftxt);
int Ntrans = array_of_transforms.size()/6;


// convert the double array back to an array of tforms.
TAffine* tfs = new TAffine[Ntrans];

for(int i=0; i<Ntrans; i++) {
    printf("Transform %3d: ", i);
    tfs[i].CopyIn( &array_of_transforms[i*6] );
    tfs[i].FromMatlab();
    tfs[i].TPrint();
    printf("\n");
    }


// Find a center for each region.   We have a map for image 'above'.  Make one for image 'below'
// Also, fill in all the pixels from b that have a map
printf("Going to allocate the counts\n");
vector<int>  pcount(Ntrans, 0);
printf("allocated the counts\n");
vector<Point> sums(Ntrans, Point(0.0,0.0));
printf("allocated the sums\n");
fflush(stdout);
MeanStd m;
for(int x=0; x<w; x++) {
    for(int y=0; y<h; y++) {
        int n = rmap[x+w*y];
        if( n-10 >= Ntrans ) {
	    printf("odd - n=%d, Ntrans=%d\n", n, Ntrans);
            exit(-1);  // illegal transform number
            }
	if( n >= 10 ) {
	    sums[n-10].x += x;
            sums[n-10].y += y;
            pcount[n-10]++;
            }
        }
    }


#if 0
for(int i=0; i<Ntrans; i++) {
    sums[i].x = sums[i].x/pcount[i];
    sums[i].y = sums[i].y/pcount[i];
    int ix = (int)sums[i].x; int iy = (int)sums[i].y;
    int patcha = -1;
    if( ix < w && iy < h )
	patcha = fold_mask_a[ix+w*iy];
    if( patcha <= 0 ) {
        printf("   WARNING: patch center not in patch at (%f %f) in patch %d of image 'above'?\n", sums[i].x, sums[i].y, patcha);
        printf("To reproduce, try this:\n");
        printf("'%s' '%s' '%s' '%s'\n", noa[1], noa[2], noa[3], noa[4]);
        }
    printf("region %d, center (%f %f), patch %d in 'above', npix=%d\n", i, sums[i].x, sums[i].y, patcha, pcount[i]);

    Point Center = sums[i];
    tfs[i].Transform( sums[i] );
    ix = (int)sums[i].x; iy = (int)sums[i].y;
    int patchb = -1;
    if( ix < w && iy < h )
	patchb = fold_mask_b[ix+w*iy];
    if( patchb <= 0 ) {
        printf("   WARNING: Point maps to (%f %f) in patch %d of image 'below'?\n", sums[i].x, sums[i].y, patchb);
        printf("To reproduce, try this:\n");
        printf("'%s' '%s' '%s' '%s'\n", noa[1], noa[2], noa[3], noa[4]);
        }
    printf("   maps to (%f %f) in patch %d of image 'below'.\n", sums[i].x, sums[i].y, patchb);

    if( patcha > 0 && patchb > 0 ) {

        fprintf(fcorr,"POINT %s::%d %f %f", noa[1], patcha, Center.x, Center.y);
        fprintf(fcorr," %s::%d %f %f\n", noa[2], patchb, sums[i].x, sums[i].y);
	}
    }
#endif

for( int i = 0; i < Ntrans; ++i ) {

	Point	p( sums[i].x / pcount[i], sums[i].y / pcount[i] );
	int		mv, ix, iy;

	ix = int(p.x);
	iy = int(p.y);

	if( ix >= 0 && ix < w && iy >= 0 && iy < h )
		mv = fold_mask_a[ix + w*iy];
	else
		mv = 0;

	fprintf( fcorr,
	"POINT %s::%d %f %f", noa[1], mv, p.x, p.y );

	tfs[i].Transform( p );

	ix = int(p.x);
	iy = int(p.y);

	if( ix >= 0 && ix < w && iy >= 0 && iy < h )
		mv = fold_mask_b[ix + w*iy];
	else
		mv = 0;

	fprintf( fcorr,
	" %s::%d %f %f\n", noa[2], mv, p.x, p.y );
}

// Free all the stuff we don't need any more.  Not important to program results, but makes leak checking easier.

delete [] tfs;
StopTiming(flog, "ReadAFileSet", started);

}


// reads in an arbitrary number of scripts, and creates a file 'corr' which contains point correlations.
int main(int argc, char* argv[])
{
bool NoFolds = false;
char *WhereToWrite = NULL;

vector<char *>noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    // process arguments here
    if( argv[i][0] != '-' )
	noa.push_back(argv[i]);
    if( strncmp(argv[i],"-f",2) == 0 )
	WhereToWrite = strdup(argv[i]+2);
    if( strcmp(argv[i],"-tr") == 0 )
	Transpose = true;
    if( strcmp(argv[i],"-nf") == 0 )
	NoFolds = true;
    }
string log_file = "align.log";

if( noa.size() < 1 ) {
    printf("reads all files specified as arguments, looking for\n");
    printf(".*deformable* <below-image> <above-image> <below-fold-map> <above-fold-map> <transform-txt-file> <map file>\n");
    exit( 42 );
    }

FILE	*flog = FileOpenOrDie( log_file.c_str(), "a" );

fprintf(flog, "args: ");
for(int i=1; i<argc; i++)
    fprintf(flog,"'%s' ", argv[i]);
fprintf(flog,"\n");

// open file for the corresponding points
FILE	*fcorr =  FileOpenOrDie(
		(WhereToWrite ? WhereToWrite : "corr_pts.txt"), "w", flog );

for( int i = 0; i < noa.size(); ++i ) {

    FILE		*fp = FileOpenOrDie( noa[i], "r", flog );
    CLineScan	LS;
    int			n;

    for( n = LS.Get( fp ); n > 0; n = LS.Get( fp ) ) {

		//printf( "line is '%s'\n", LS.line );

		char *p = strtok( LS.line, " \"\n" );

		if( p != NULL &&
			(strstr( p, "deformable" ) || strstr( p, "ptest" )) ) {

			// line starts with a token containing 'deformable' or 'ptest'.  THis is an alignment command
			vector<char *>	args;

			for(
				p = strtok( NULL, " \"\n" );
				p != NULL;
				p = strtok( NULL, " \"\n" ) ) {

				//printf( "Token is '%s'\n", p );

				if( p[0] != '-' )
					args.push_back( p );
			}

			ReadAFileSet( args, flog, fcorr, NoFolds );
		}
	}

    fclose( fp );
}

fprintf(flog, "Normal completion for align run.\n");
fclose(flog);
printf("Normal completion for align run.\n");
return 0;
}


