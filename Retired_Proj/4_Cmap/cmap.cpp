

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"TAffine.h"

#include	<string.h>

#include	<stack>
using namespace std;


static bool Transpose = false;      // Use translation only?
static bool WithinSection = false; // overlap within a section
static bool AllFiles = false;      // use files for everything
static int OVERLAP = 75;  //expect about this much overlap, in pixels
static int RADIUS = 256;  // correction for alignment should be found within this distance
static double THRESHOLD = 0.25;  // lowest correlation considered a match






bool SomeArea(int sx, int sy, void *a)
    {return sx > 1 && sy > 1;}

void ReadTransforms(const char *name, vector<TAffine> &ts)
{
FILE	*ft = FileOpenOrDie( name, "r" );

for(;;) {
    double a,b,c,d,e,f;
    if( fscanf(ft, "%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f) != 6 )
    break;
    TAffine t(a,c,e,b,d,f);  // file is in matlab order
    ts.push_back(t);
    }
printf("read %d transforms from file '%s'\n", ts.size(), name);
//for(int i=0; i<ts.size(); i++)
   //ts[i].TPrint();
fclose(ft);
}
//Reads two map files amd transform lists.
int main(int argc, char* argv[])
{

FILE	*flog = FileOpenOrDie( "cmap.log", "a" );

if( argc < 3 ) {
    printf("Usage: cmap < a directory> <b directory>\n");
    exit( 42 );
    }
bool WriteBack = false;   // overwrite the original?

time_t t0 = time(NULL);
char atime[32];
strcpy(atime, ctime(&t0));
atime[24] = '\0';  // remove the newline
fprintf(flog,"two: %s ", atime);
bool TestingNormCorr = false;
for(int i=1; i<argc; i++) {
    // process arguments here
    fprintf(flog,"%s ", argv[i]);
    if( strcmp(argv[i],"-o") == 0 && i+1 < argc )
    OVERLAP = atoi(argv[i+1]);
    if( strcmp(argv[i],"-t") == 0 && i+1 < argc )
    THRESHOLD = atof(argv[i+1]);
    if( strcmp(argv[i],"-r") == 0 && i+1 < argc )
    RADIUS = atoi(argv[i+1]);
    if( strcmp(argv[i],"-f") == 0 )
    WriteBack = true;
    if( strcmp(argv[i],"-tr") == 0 )
    Transpose = true;
    if( strcmp(argv[i],"-ws") == 0 )
    WithinSection = true;
    if( strcmp(argv[i],"-af") == 0 )
    AllFiles = true;
    if( strcmp(argv[i],"-tnc") == 0 )
    TestingNormCorr = true;
    }
if( AllFiles && argc < 7 ) {
    printf("Not enough args for '-af' option.\n");
    exit( 42 );
    }
fflush(flog);

// Read the two mapping files
uint8 *amap;
uint8 *bmap;
char file_name[256];
sprintf(file_name, "%s/map.tif", argv[1]);
uint32 w1, w2, h1, h2;
amap = Raster8FromTif(file_name, w1, h1, flog, Transpose);
sprintf(file_name, "%s/map.tif", argv[2]);
bmap = Raster8FromTif(file_name, w2, h2, flog, Transpose);
if( w2 != w1 || h2 != h1 ) {
    printf("Maps are different size than input\n");
    fprintf(flog,"Maps are different size than input\n");
    exit( 42 );
    }

// Now read the two sets of transforms
vector<TAffine> at, bt;
sprintf(file_name, "%s/t.txt", argv[1]);
ReadTransforms(file_name, at);
sprintf(file_name, "%s/t.txt", argv[2]);
ReadTransforms(file_name, bt);

// create an empty color image
vector<uint32>	raster( w1 * h1, 0xFF000000 );	// set alpha chan

// Now compare the files
int nmb = 0;  // not mapped in both
int ao = 0;   // mapped in a only
int bo = 0;   // mapped in b only
int np = 0;   // number printed
int plimit = 0;
MeanStd stats;
double peak = 0.0; // peak error
for(int i=0; i<w1*h1; i++) {

    int y = i / w1;
    int x = i - w1 * y;

    if( amap[i] == 0 && bmap[i] == 0 ) {
    nmb++;
    raster[i] |= 0x808080;  // medium gray
        }
    else if( amap[i] != 0 && bmap[i] == 0 ) {
    ao++;
    raster[i] |= 0xA0A0A0;  // lighter gray
    }
    else if( amap[i] == 0 && bmap[i] != 0 ) {
    bo++;
    raster[i] |= 0x606060;  // darker gray
    }
    else { // both are mapped
    Point pa(x,y), pb(x,y);
    at[amap[i]-10].Transform( pa );
    bt[bmap[i]-10].Transform( pb );
    double d2 = (pa.x-pb.x)*(pa.x-pb.x) + (pa.y-pb.y)*(pa.y-pb.y);
        double d = sqrt(d2);
        if( np++ < plimit )
        printf("%d,%d -> %d(%9.2f %9.2f) : %d(%9.2f %9.2f) = %f\n", x, y, amap[i], pa.x, pa.y, bmap[i], pb.x, pb.y, d);
        stats.Element(d);
        peak = max(peak, d);
    if( d <= 1.0 ) {
        int color =  127 + int(d*10)*12;
            raster[i] |= (color << 8);  // green
        }
        else if( d < 3.0 ) {
        int color =  127 + int((d-1)*10)*6;
            raster[i] |= (color << 16);  // blue
        }
        else {
            int color = 127 + int(d-3)*6;
            if (color > 255) color = 255;
        raster[i] |= color;
        }
    }
    }

// Write it out
Raster32ToTifRGBA( "emap.tif", &raster[0], w1, h1 );

printf(" %d not mapped, %d a only, %d b only, %d both\n", nmb, ao, bo, stats.HowMany() );
double mean = 0.0, std = 0.0;
if( stats.HowMany() != 0 ) {
    stats.Stats(mean, std);
    printf("Mean error %f, std dev of error %f, peak error %f\n", mean, std, peak);
    }
int limit = w1*h1/100;   // should cover the same to within a %
if( ao > limit || bo > limit || mean > 3.0 || peak > 20.0 )
    return 1;
return 0;  // return OK
}
