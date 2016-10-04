

#include	"mrc.h"

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"

#include	<string.h>

#include	<set>
using namespace std;






void FixName(char *name)
{
if (name[0] == '/') return;  // do nothing if rooted.
char *copy = strdup(name);
copy[8] = '\0';
char temp[2048];
sprintf(temp,"/groups/em/leginon/Data/%s/rawdata/%s", copy, name);
strcpy(name, temp);
free(copy);
}



int main(int argc, char **argv)
{
if( argc < 2 ) {
    printf("Usage: qcheck <list of images with offsets>\n");
    exit( 42 );
    }
vector<char *> noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    if( argv[i][0] != '-' )
    noa.push_back(argv[i]);
    else {
    printf("Unknown option %s\n", argv[3]);
    return 42;
    }
    }
for(int narg=0; narg<noa.size(); narg++) {
    printf("------ Working on %s -----------------\n", noa[narg]);
    double xscale = 1.0;
    double yscale = 1.0;
    if( strstr(noa[narg],"45Tilt") != NULL ) { // we have tilts
    //if( strstr(noa[narg],"90Rot") != NULL )
        //yscale = 1/sqrt(2.0);
    //else
        xscale = 1/sqrt(2.0);
    }
    printf( "Scales are x %f, y %f\n", xscale, yscale );

    FILE	*fp = FileOpenOrDie( noa[narg], "r" );

    // read the file once, to get the image sizes, and layer IDs
    set<int> layers;
    uint32 w=0,h=0;
    for(int n=0; ;n++) {
    char fname[2048];
    double x,y;
    int z;
    if( fscanf(fp,"%s %lf %lf %d", fname, &x, &y, &z) != 4 )
        break;
    if( w == 0 ) {
        vector<uint16*>	vras;
        FixName(fname);
        ReadRawMRCFile( vras, fname, w, h, stdout, false );
        FreeMRC( vras );
        }
    layers.insert(z);
    }
    // get the root part of the input file name.  If input is blah.txt,
    // will generate blah.N.png, where N is each layer number
    char *root = strdup(noa[narg]);
    char *ext = strstr(root, ".txt");
    if( ext == NULL ) {
    printf("Expected input name to contain '.txt'\n");
    return 42;
    }
    *ext = '\0';

    // OK, now go through all layers
    set<int>::iterator it;
    for(it = layers.begin(); it != layers.end(); it++) {
    printf("Starting layer %d\n", *it);
    rewind(fp);
    double xmin =  1.0E30, ymin =  1.0E30;
    double xmax = -1.0E30, ymax = -1.0E30;
    for(int n=0; ;n++) {
        char fname[2048];
        double x,y;
        int z;
        if( fscanf(fp,"%s %lf %lf %d", fname, &x, &y, &z) != 4 )
        break;
            x *= xscale;
            y *= yscale;
        if( z == *it ) {
        xmax = max(xmax, x+w); ymax = max(ymax, y+h);
        xmin = min(xmin, x  ); ymin = min(ymin, y);
        }
        }
    printf("Limits x=[%f %f], y = [%f %f]\n", xmin, xmax, ymin, ymax);
    int ow = int(xmax - xmin + 1.5);
    int oh = int(ymax - ymin + 1.5);
    printf("Full image %d by %d\n", ow, oh);
    int rw = ((ow-1) >> 4) + 1;
    int rh = ((oh-1) >> 4) + 1;
    printf("Reduced image %d by %d\n", rw, rh);
    int xoff = int(xmin);
    int yoff = int(ymin);

    // re-read the file, and write the image
    vector<uint32>frame(rw*rh,0);  // use ints, since we may add many 16 bit values
    vector<int> counts(rw*rh, 0);  // number of times each pixel is written to
    rewind(fp);
    for(int n=0; ;n++) {
        char fname[2048];
        double x,y;
        int z;
        if( fscanf(fp,"%s %lf %lf %d", fname, &x, &y, &z) != 4 )
        break;
        if (z != *it) continue;
        int x0 = int(x*xscale);
        int y0 = int(y*yscale);
        FixName(fname);
        vector<uint16*>	vras;
        uint16*			raster;
        ReadRawMRCFile( vras, fname, w, h, stdout, false );
        raster = vras[0];
        for(int y=0; y<h; y++) {
        int ry = (y0 + y - yoff) >> 4;
        int row_origin = rw * ry;  // origin of line in 'frame'
        int raster_y = w*y;      // origin of line in 'raster'
        for(int x=0; x<w; x++) {
            int rx = (x0 + x - xoff) >> 4;
            frame[row_origin + rx] += raster[raster_y+x];
            counts[row_origin + rx]++;
            }
        }
        FreeMRC( vras );
    //if (n > 100)  // debugging only
        //break;
        }
    int np = rw*rh;
    for(int i=0; i<np; i++) {
        if( counts[i] != 0 )
        frame[i] = frame[i]/counts[i];
        }
    MeanStd m;
    for(int i=0; i<np; i++) {
        if( frame[i] != 0 )
        m.Element(frame[i]);
        }
    double mean, std;
    m.Stats(mean, std);
    printf("Mean is %f, std = %f\n", mean, std);

    vector<uint8> out(rw*rh);

    for(int i=0; i<rw*rh; i++) {
        int pix = 127 + int((frame[i] - mean)/std * 60.0);
        if( pix < 0 )
            pix = 0;
        else if( pix > 255 )
            pix = 255;
        out[i] = pix;
        }

    char oname[2048];
    sprintf(oname,"%s.%d.png", root, *it);
        FILE *ft = fopen(oname,"w");
        if( ft != NULL )
        fclose(ft);
        else
         sprintf(oname,"/tmp/%s.%d.png", root, *it);
    Raster8ToPng8(oname, &out[0], rw, rh);  // this is the fold mask for drawing
    }
    fclose(fp);
    }
return 0;
}


