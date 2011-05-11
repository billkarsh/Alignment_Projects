

#include	"mrc.h"

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"






#if 0	// not used here, but interesting ideas on fold masking

// Here we convert an input image to a map where 0 = on fold, 1=region 1, 2=region 2, etc.
// Normally the pipeline does this, but we do it here for standalone testing.
void ImageToFoldMap(Picture &pic, uint8* FoldMask, bool remove_low_contrast = false,
 bool one_region = false, double FoldMaskThresholdOverride = 0.0)  //pipeline does this in normal operation
{

// Here, define the connected region(s) on the above layer.
int w = pic.w;
int h = pic.h;
uint8 *raster = pic.raster;
int npixels = w * h;

// if they asked for one region, do that
if( one_region ) {
    printf("Generating one uniform region\n");
    for(int i=0; i<npixels; i++)
	FoldMask[i] = 1;
    return;
    }

// First, find the mean and standard deviation of the 'real' non-saturated pixels.
// Ignore any that are too close to saturated light or dark.
int SAT=3;  // how close to the boundary do you need to be to be considered 'saturated'
MeanStd m;
for(int i=0; i<npixels; i++) {
    uint8 pix = raster[i];
    if( pix >= SAT && pix <= 255-SAT )
	m.Element(raster[i]);
    }
printf("%d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);

if(m.HowMany()/double(npixels) < 0.9) {
    printf("Saturated image!  Retrying with SAT=1\n");
    SAT = 1;
    m.Reset();
    for(int i=0; i<npixels; i++) {
	uint8 pix = raster[i];
	if( pix >= SAT && pix <= 255-SAT )
	    m.Element(raster[i]);
	}
    printf("%d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);
    }

double mean, std;
m.Stats(mean, std);  // find existing statistics
printf("Of the above image points, mean= %f and std dev = %f\n", mean, std);

// If pure white is outside the range of practical values, then it too
// should be ignored, as it is not useful info.  But if 255 is within the
// 'practical' range, we don't want to cut these pixels out
if( mean + 2.5*std < 255 ) {
    printf("Removing white pixels\n");
    for(int i=0; i<npixels; i++) {
	uint8 pix = raster[i];
	if (pix > 255-SAT)           // too bright is just as bad as too dark - not useful info
	    raster[i] = 0;
	}
    }

// Now find the connected regions.  We want to ignore very black pixels as folds, and those
// near them (distance = D) as unreliable, and might connect regions which should be disconnected.
// However, if the image as a whole is too black, this will result in lots of disconnected
// regions, as there will be many black pixels and many near them.  So reset the parameters in
// these cases.
int D=10;
int nbig = 0;
vector<ConnRegion> cr;
double thresh = 4.0;          // we would like 4 std dev, if possible
if( mean - thresh*std < 1 ) {  // but if not, pick a value that is feasible
    thresh = mean/std*0.95;   // set to 95% of the way to 0
    printf("Forced to reduce threshold to %f std dev.\n", thresh);
    if( thresh < 2.0 ) {       // we'll get too many black, and fragment the area
	thresh = (mean - 0.5)/std;  // set threshold to a pixel value of 0.5 (on scale of 255)
	printf("Desperate measures.  Changing threshold to %f, so only v=0 pixels are out.\n", thresh);
        printf("Also disabling black pixel expansion\n");
	D = 0;
	}
    }
if( FoldMaskThresholdOverride != 0.0 ) {
    thresh = FoldMaskThresholdOverride;
    printf("Explicit over-ride of threshold; set to %f\n", thresh);
    }

// pixels to vector
vector<double>	v(npixels);
for(int i=0; i<npixels; i++) {
    int y = i / w;
    int x = i - w * y;   // pixels next to background pixels should also be black
    int pix = raster[i];
    v[i] = (pix-mean)/std;
    }

if( remove_low_contrast ) {
    // If there are chunks with no contrast, get rid of them.
    int nx = (w-128)/32+1;
    double dx = (w-128)/double(nx);  // idea is to tile with overlapping 128x128 squares
    int ny = (w-128)/32+1;
    double dy = (w-128)/double(ny);
    vector<int> zap_me;        // save list of pixels to be zapped, since cannot zap on the fly without
			       // munging the overlapping squares
    for(double y=0; y< (h-127); y += dy) {
	int ymin = min(int(y), h-128);  // just in case
	int ymax = ymin+127;
	for(double x=0; x< (w-127); x += dx) {
	    int xmin = min(int(x), w-128);
	    int xmax = xmin+127;
	    vector<double>local(128*128);
	    //printf("xmin, ymin = (%d %d)\n", xmin, ymin);
	    for(int ix=xmin; ix <= xmax; ix++)
		for(int iy=ymin; iy <= ymax; iy++)
		    local[(ix-xmin) + 128*(iy-ymin)] = v[ix + w*iy];
	    if( IsLowContrast(local, std) ) {
		for(int ix=xmin; ix <= xmax; ix++)
		    for(int iy=ymin; iy <= ymax; iy++)
			zap_me.push_back(ix + w*iy);
		}
	    }
	}
    for(int i=0; i<zap_me.size(); i++)
	v[zap_me[i]] = -thresh - 100.0;  // sure to be bad
    }

// for those that should not be included, remove all pixels within distance D as well.  This will
// remove tiny filaments connected big regions together.
vector<int> remove;  // indices of points to be removed.
for(int i=0; i<npixels; i++)
    if( v[i] < -thresh )
	remove.push_back(i);
for(int ii=0; ii<remove.size(); ii++) {
    int i = remove[ii];
    int y = i / w;
    int x = i - w * y;
    int x0 = max(x-D,0);
    int x1 = min(x+D, w-1);
    int y0 = max(y-D,0);
    int y1 = min(y+D,h-1);
    for(int xx = x0; xx <= x1; xx++) {
	for(int yy = y0; yy <= y1; yy++) {
	    v[xx+w*yy] = -thresh - 1.0;   // set to a value that is more than enough black
            }
        }
    }

// Now find the connected regions

	for( int i = 0; i < npixels; ++i ) {

		if( v[i] > -thresh ) {

			ConnRegion	c;
			int			npts;

			npts = Propagate( c.pts, v, w, h, i,
					-thresh, -thresh - 1.0 );

			printf(
			"ImageToFoldMap: ConnRegion with %d pixels.\n", npts );

			if( npts > 90000 ) {	// want 100k, but we shrank it
				cr.push_back( c );
				++nbig;
			}
		}
	}

 // Final accounting

  SetBoundsAndColors( cr, FoldMask, w, h );
}

#endif








// ------------------ Read a DM3 file
void Ind(int j)
{
for(int i=0; i<j*2; i++)
    printf(" ");
}

// reads one value of the specified type from a DM3 file.
// In the case of LONG, returns a value
int ReadOneDM3(FILE *fp, int what, int indent)
{
int bytes = 0;
Ind(indent);
switch (what) {
    case 2: printf("SHORT:"); bytes = 2; break;
    case 3: printf("LONG :"); bytes = 4; break;
    case 4: printf("USHRT:"); bytes = 2; break;
    case 5: printf("ULONG:"); bytes = 4; break;
    case 6: printf("FLOAT:"); bytes = 4; break;
    case 7: printf("DBLE :"); bytes = 8; break;
    case 8: printf("BOOL :"); bytes = 1; break;
    case 10:printf("OCTET:"); bytes = 1; break;
    default: printf("Oops - got %d\n", what); exit( 42 );
    }
vector<uint8>	val(bytes);
size_t items = fread(&(val[0]), sizeof(uint8), bytes, fp);
for(int i=0; i<bytes; i++)
    printf(" %d", val[i]);
int result = 0;
if( what == 6 ) {
    float foo;
    uint8 *cp = (uint8*)&foo;
    cp[0] = val[0]; cp[1] = val[1]; cp[2] = val[2]; cp[3] = val[3];
    printf(" (%16.8f)", foo);
    }
else if( what == 7 ) {
    double foo;
    uint8 *cp = (uint8*)&foo;
    for(int i=0; i<8; i++)
	cp[i] = val[i];
    printf(" (%16.8f)", foo);
    }
else if( what == 3 ) {
    uint8 *cp = (uint8*)&result;
    cp[0] = val[0]; cp[1] = val[1]; cp[2] = val[2]; cp[3] = val[3];
    printf(" (%d)", result);
    }
printf("\n");
return result;
}

void ReadArrayOfStructs(FILE *fp, vector<int> &defs, int indent)
{
int nf = defs[3];
Ind(indent); printf("Array of structs; struct has %d fields\n", nf);
vector<int> in_record( nf );
for(int i=0; i<nf; i++)
    in_record[i] = defs[5+i*2];
int nitems = defs[4+2*nf];
Ind(indent); printf("%d structures in array\n", nitems);
for(int i=0; i<nitems; i++) {
    Ind(indent); printf("Reading entry %d\n", i);
    for(int j=0; j<in_record.size(); j++)
	ReadOneDM3(fp, in_record[j], indent);
    }
}

// Read N elements of type 'what'. Special cased since this is the main data object
void ReadSimpleArray(FILE *fp, int what, int n, vector<int> &v, int indent)
{
if( what == 3 ) { // array of longs
    v.resize(n);  // make sure it's big enough
    size_t items = fread(&(v[0]), sizeof(int), n, fp);
    for(int i=0; i<min(5,n); i++) {
	Ind(indent); printf("32 bit[%d} = %d\n", i, v[i]); }
    }
else if( what == 2 ) {  // array of shorts
    vector<uint16> dat(n);
    size_t items = fread(&(dat[0]), sizeof(int), n, fp);
    for(int i=0; i<min(5,n); i++) {
	Ind(indent); printf("16 bit[%d} = %d\n", i, dat[i]); }
    }
else { // anything else, read one at a time
    for(int l=0; l<n; l++)
	ReadOneDM3(fp, what, indent);
    }
}

void ReadDM3TagGroup(FILE *fp, vector<int> &v, uint32 &w, uint32 &h, int indent = 0)
{
// read 6 bytes to define the group
uint8	header[6];
size_t items = fread(header, sizeof(uint8), 6, fp);
//for(int i=0; i<256; i++)  // swap the words
     //header[i] = (header[i] << 16) | ((header[i] >> 16) & 0xFFFF);
Ind(indent);
printf("First 6 bytes of group: %d %d %d %d %d %d\n", header[0], header[1], header[2], header[3], header[4], header[5]);
int ntags = (header[2] << 24) + (header[3] << 16) + (header[4] << 8) + header[5];
Ind(indent); printf("Reading %d tags\n", ntags);
for(int i=0; i<ntags; i++) {
    Ind(indent);printf("Reading tag %d\n", i);
    uint8	TagHeader[6];
    size_t items = fread(TagHeader, sizeof(uint8), 3, fp);
    Ind(indent); printf("First 3 bytes of Tag: %d %d %d\n", TagHeader[0], TagHeader[1], TagHeader[2]);
    int label_size = (TagHeader[1] << 8) + TagHeader[2];
    //printf("String has size %d\n", label_size);
    vector<char> label(label_size);
    fread((void*) &(label[0]), sizeof(char), label_size, fp);
    label.push_back(0);  // add a terminator
    //for(int j=0; j<label_size; j++)
	//printf("Char '%c', int %d\n", label[j], label[j]);
    Ind(indent); printf("string is '%s'\n", &(label[0]) );
    if( TagHeader[0] == 20 )
	ReadDM3TagGroup(fp, v, w, h, indent+1);  // it's a hierarchical group
    else if( TagHeader[0] == 21 ) {
        uint8	TagType[8];
        size_t items = fread(TagType, sizeof(uint8), 8, fp);
        Ind(indent);printf("First 4 bytes: %d %d %d %d\n", TagType[0], TagType[1], TagType[2], TagType[3]);
        Ind(indent);printf(" next 4 bytes: %d %d %d %d\n", TagType[4], TagType[5], TagType[6], TagType[7]);
        int ndefs = (TagType[4] << 24) + (TagType[5] << 16) + (TagType[6] << 8) + TagType[7];
        Ind(indent);printf("Ndefs = %d\n", ndefs);
        // read as bytes and re-assemble into words
        vector<uint8>	temp(4*ndefs);
        items = fread(&(temp[0]), sizeof(uint8), ndefs*4, fp);
        vector<int> defs(ndefs);
        for(int kk=0; kk<ndefs*4; kk += 4) {
            int k = kk/4;
            defs[k] = (temp[kk] << 24) + (temp[kk+1] << 16) + (temp[kk+2] << 8) + temp[kk+3];
	    Ind(indent);printf("Def word[%d] = %d\n", k, defs[k]);
            }
        if( defs[0] == 15 ) {
	     Ind(indent);printf("This is a struct with %d components\n", defs[2]);
             for(int l=0; l<defs[2]; l++) {
		Ind(indent);printf("Component %d is a %d\n", l, defs[4+2*l]);
                ReadOneDM3(fp, defs[4+2*l], indent);
		}
	     }
        else if( defs[0] == 20 && defs[1] == 15 )
	     ReadArrayOfStructs(fp, defs, indent);
        else if( defs[0] == 20 ) {
	     Ind(indent);printf("Array of type %d, %d entries\n", defs[1], defs[2]);
             ReadSimpleArray(fp, defs[1], defs[2], v, indent);
	     }
        else if( ndefs == 1 ) {
	     Ind(indent);printf("Read a simple type\n");
             int got = ReadOneDM3(fp, defs[0], indent);  // only returned for ints
             if (strcmp(&label[0], "width") == 0)        // but w and h are integer valued
		w = got;
             if( strcmp(&label[0], "height") == 0 )
		h = got;
	     }
        else {
	     Ind(indent);printf("Not known %d\n");
	     exit( 42 );
             }
	}
    else {
	Ind(indent);printf("expected only 20 or 21 here, got %d\n", TagHeader[0]);
	exit( 42 );
	}
    }
}
uint16* ReadADM3File(const char *name, uint32 &w, uint32 &h, FILE *flog, bool Transpose)
{
printf("OK, will open dm3 file '%s'\n", name);

FILE	*fp = FileOpenOrDie( name, "r", flog );

// read the first 12 bytes - the header
uint8	header[12];
size_t items = fread(header, sizeof(uint8), 12, fp);
//for(int i=0; i<256; i++)  // swap the words
     //header[i] = (header[i] << 16) | ((header[i] >> 16) & 0xFFFF);
printf("First 4 bytes: %d %d %d %d\n", header[0], header[1], header[2], header[3]);
printf(" next 4 bytes: %d %d %d %d\n", header[4], header[5], header[6], header[7]);
vector<int> pixels;
ReadDM3TagGroup(fp, pixels, w, h);
int np = w*h;
printf("pixel array is %d pixels, w*h=%d pixels\n", pixels.size(), np );
if( np != pixels.size() ) {
    printf("Expected these two to be the same\n");
    exit( 42 );
    }
// find the mean and std dev
MeanStd m;
for(int i=0; i<w*h; i++)
    m.Element(pixels[i]);
double mean, lstd;
m.Stats(mean, lstd);
printf("Mean= %f, std dev = %f\n", mean, lstd);
uint16* raster = (uint16 *)malloc(np*sizeof(uint16));
if( raster == NULL ) {
    printf("Read DM3 malloc failed\n");
    exit( 42 );
    }
// convert to a mean of 32678 and a std of 10000
for(int i=0; i<np; i++) {
    double pix = 32768 + (pixels[i] - mean)/lstd*10000.0;
    if (pix < 0) pix = 0.0;
    if (pix > 65535) pix = 65535;
    raster[i] = uint16(pix);
    }
return raster;
}
// ------------------ End of MRC file reading




int main(int argc, char **argv)
{
if( argc < 3 ) {
    printf("Usage: importMRC <file.mrc> <file.png> [options]\n");
    exit( 42 );
    }
double EnergyThreshold = 10.0;
int Ngauss = 2;
int trim = 0;
vector<char *> noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    if( argv[i][0] != '-' )
	noa.push_back(argv[i]);
    else if( strncmp(argv[i], "-energy=",8) == 0 ) {
	EnergyThreshold = atof(argv[i]+8);
        printf("Energy Threshold now %f\n", EnergyThreshold);
	}
    else if( strncmp(argv[i], "-gauss=",7) == 0 ) {
	Ngauss = atoi(argv[i]+7);
        printf("Now using %d gaussians for image normalization\n", Ngauss);
	}
    else if( strncmp(argv[i],"-trim=",6) == 0 ) {
	trim = atoi(argv[i]+6);
        printf("Trimming %d pixels off each edge.\n", trim);
	}
    else {
	printf("Unknown option %s\n", argv[3]);
	return 42;
	}
    }
bool Transpose = false;
// if it ends with .mrc, use the MRC file reader, otherwise bomb
char *mrcp = strstr(noa[0],".mrc");
char *dm3p = strstr(noa[0],".dm3");
uint16 *raster;
uint32 w,h;
if( mrcp != NULL && *(mrcp+4) == 0 ) {
    uint8 *raster = ReadAnMRCFile(noa[0], w, h, stdout, Transpose, Ngauss);
    Raster8ToPng8(noa[1], raster, w, h);
    return 0;
    }
else if( dm3p != NULL && *(dm3p+4) == 0 )
    raster = ReadADM3File(noa[0], w, h, stdout, Transpose);
else {
    printf("Expected .mrc or .dm3 file, got %s\n", noa[0]);
    return 42;
    }

// trim the edges off the picture
int neww = w-2*trim;
int newh = h-2*trim;
int newp = neww*newh;
vector<double> vals(newp);
for(int x=0; x<neww; x++)
    for(int y=0; y<newh; y++)
	vals[x+neww*y] = raster[trim + x + w*(trim+y)];

// convert to a mean of 0 and a std dev of 1
Normalize(vals);

// Convert to 8 bit
vector<uint8> rast(newp);
for(int i=0; i<newp; i++) {
    int px = int(vals[i]*35 + 127);
    if (px < 0) px = 0;
    if (px > 255) px = 255;
    rast[i] = px;
    }

Raster8ToPng8(noa[1], &(rast[0]), neww, newh);

return 0;
}


