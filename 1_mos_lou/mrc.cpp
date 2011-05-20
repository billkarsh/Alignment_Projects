// ------------------ Read an .mrc file
vector<uint16*> ReadRawMRCFile(const char *name, uint32 &w, uint32 &h, FILE *flog, bool Transpose)
{
printf("OK, will open mrc file '%s'\n", name);
FILE *fp = fopen(name,"r");
if (fp == NULL) {
    printf("Cannot open '%s' for read\n", name);
    fprintf(flog,"Cannot open '%s' for read\n", name);
    exit(42);
    }
// read the first 1024 bytes - the header
int header[256];
size_t items = fread(header, sizeof(int), 256, fp);
//for(int i=0; i<256; i++)  // swap the words
     //header[i] = (header[i] << 16) | ((header[i] >> 16) & 0xFFFF);
printf("Header: %d %d %d %d", header[0], header[1], header[2], header[3]);
printf(" %d %d %d %d\n", header[4], header[5], header[6], header[7]);
w = header[0];
h = header[1];
int np = w*h;
vector<uint16*> result;
printf("MRC file has %d images.\n", header[2]);
if (header[3] == 6) {
    for(int i=0; i<header[2]; i++) {
	// read that number of shorts
	uint16* raster = (uint16 *)malloc(np*sizeof(uint16));
	if (raster == NULL) {
	    printf("Read MRC malloc failed\n");
	    exit(42);
	    }
	items = fread(raster, sizeof(uint16), np, fp);
	if (items != np) {
	    printf("Read MRC data read failed\n");
	    exit(42);
	    }
	result.push_back(raster);
	}
    }
else if (header[3] == 2) {  // floats
    for(int i=0; i<header[2]; i++) {
	// read that number of floats
	float* raster = (float *)malloc(np*sizeof(float));
	if (raster == NULL) {
	    printf("Read MRC malloc failed\n");
	    exit(42);
	    }
	items = fread(raster, sizeof(float), np, fp);
	if (items != np) {
	    printf("Read MRC data read failed\n");
	    exit(42);
	    }
        // convert to uint16
	uint16* r2 = (uint16 *)malloc(np*sizeof(uint16));
	if (r2 == NULL) {
	    printf("Read MRC malloc failed\n");
	    exit(42);
	    }
        for(int k=0; k<np; k++)
	    r2[k] = int(raster[k]);
        free(raster);
	result.push_back(r2);
	}
   }
else if (header[3] == 0) {  // bytes.  Assuming unsigned for here
    for(int i=0; i<header[2]; i++) {
	// read that number of floats
	unsigned char* raster = (unsigned char *)malloc(np*sizeof(float));
	if (raster == NULL) {
	    printf("Read MRC malloc failed\n");
	    exit(42);
	    }
	items = fread(raster, sizeof(unsigned char), np, fp);
	if (items != np) {
	    printf("Read MRC data read failed\n");
	    exit(42);
	    }
        // convert to uint16
	uint16* r2 = (uint16 *)malloc(np*sizeof(uint16));
	if (r2 == NULL) {
	    printf("Read MRC malloc failed\n");
	    exit(42);
	    }
        for(int k=0; k<np; k++)
	    r2[k] = raster[k];
        free(raster);
	result.push_back(r2);
	}
   }
else {
   printf("Reading MRC file; expected mode 0(bytes), 2(32 bit floats), or 6 (16 bit ints), but got %d\n", header[3]);
   exit(42);
   }
fclose(fp);
return result;
}

// fit N gaussians to the data, then find limits.
void NGaussHist(vector<double> x, vector<int> &histo, const int Ngauss, double &rmin, double &rmax, double init_width=1000.0)
{
VecDoub y,s;  // create the y array (same data but double), and the s array 
              // (points uniformly weighted; all have same std dev)
double sum = 0.0;  // used for checking for plausibility later
for(int i=0; i<x.size(); i++) {
    y.push_back(histo[i]);
    s.push_back(1.0);
    sum += histo[i];
    }

// Remove over-represented bins at the extremes.  Sum/y.size() is what you would expect with a uniform
// distribution.  Since the real distribution is (or should be) central, there should be many less than
// this counts.  But sometimes there are large black or white areas, and they screw up the fits.  So
// look for these and remove them.
for(int i=0; i<5 && y[0] > sum/y.size(); i++) { // remove up to 5 bins on the bottom, if they are over-represented
    printf("Dark entry %d over-represented - got %f, expected %f . Romoved for fitting.\n", i, y[0], sum/y.size() );
    sum -= y[0];
    x.erase(x.begin());  // OK, since we pass by value
    y.erase(y.begin());
    s.erase(s.begin());
    }
for(int i=0; i<5 && y.back() > sum/y.size(); i++) { // remove up to 5 bins on the top, if they are over-represented
    printf("Light entry %d over-represented - got %f, expected %f . Romoved for fitting.\n", i, y.back(), sum/y.size() );
    sum -= y.back();
    x.erase(x.end()-1);
    y.erase(y.end()-1);
    s.erase(s.end()-1);
    }

// Sometimes just on value sticks out like a sore thumb (a grey patch on
// the image, for example).
vector<double> sh = y;
sort(sh.begin(), sh.end());
if (sh.size() >= 2) {
    printf("Two largest histo values %.1f and %.1f\n", sh.back(), *(sh.end()-2) );
    }

// We will make up to three passes of a ransac like algorithm.  If after fitting
// one gaussian, the worst point is still the same as it was before, then that
// point is an outlier, and should be removed and the process re-started.
bool try_again = true;
for(int pass=0; pass < 3 && try_again; pass++) {
    try_again = false;
    // Now try a 1,2,..,N gaussian fit
    VecDoub residual = y;
    VecDoub a;
    vector<int> peaks(Ngauss+1);
    for(int N=1; N<=Ngauss; N++) {
        printf("--- try %d gaussians, max %d ---\n", N, Ngauss);
        a.resize(N*3);
        // find the biggest peak in the residual, and add that
	double biggest = 0.0;
	int    bigi = 0;  // just to avoid warnings
	for(int i=0; i<residual.size(); i++) {
	    if (abs(residual[i]) > biggest) {
		biggest = abs(residual[i]);
		bigi = i;
		}
	    }
	printf("Biggest remaining residual has peak of %f at x=%f\n", biggest, x[bigi]);
        peaks[N] = bigi;
        if (N >= 2 && peaks[N] == peaks[N-1]) {
	    printf("*** Odd - peak not removed by fitting?  Assume outlier and re-fit\n");
            if (bigi > 0 && bigi < x.size()-1)
		y[bigi] = (y[bigi-1] + y[bigi+1])/2.0; // set to average, if possible
            else
                y[bigi] = 0;                           // otherwise set to 0
            try_again = true;
            }
	a[(N-1)*3  ] = residual[bigi] * 0.99;
	a[(N-1)*3+1] = x[bigi];
	a[(N-1)*3+2] = init_width;  // initial width usually 1000 for EM images
	//for(int i=0; i<a.size(); i++)
	    //printf("a[%d] = %f\n", i, a[i]);
	Fitmrq f(x,y,s,  a, fgauss);
	try {
	    f.fit();
            }
        catch (int) {
            printf("Caught exception in fit\n");
	    for(int i=0; i<a.size(); i++)
	        printf("a[%d] = %f\n", i, a[i]);
            if (N == 2) { // the two gaussian fit failed.  Need a backup; first see if
                          // one gaussian fit was reasonable. See if peak is within the
                          // data bounds and amplitude is OK
                if (x[0] < a[1] && a[1] < x[x.size()-1] && a[0] > 0 && a[0] < sum) {
		    rmin = a[1] - 4/sqrt(2.0)*abs(a[2]);
		    rmax = a[1] + 4/sqrt(2.0)*abs(a[2]);
		    rmin = max(rmin, x[0]);  // rmin should not be smaller than the first x
		    rmax = min(rmax, x[x.size()-1]);  // and rmax not off the top end
		    printf("Reverting to 1 gaussian fit, +- 4 sigma: %f %f\n", rmin, rmax);
		    }
                else {
                    rmin = x[0];
                    rmax = x[x.size()-1];
                    printf("One gaussian fit looked bad, too. Use full range: %f %f\n", rmin, rmax);
                    }
		return;
		}
	    }
	rmin =  1.0e30;
	rmax = -1.0e30;
	printf("after fit:\n");
	for(int i=0; i<a.size(); i += 3) {
	    double top = f.a[i+1] + 4/sqrt(2)*abs(f.a[i+2]);  // four sigma (sqrt(2) since width un-normalized
	    rmax = max(rmax, top);
	    double bot = f.a[i+1] - 4/sqrt(2)*abs(f.a[i+2]);  // four sigma (sqrt(2) since width un-normalized
	    rmin = min(rmin, bot);
	    printf("%16.1f * exp( ((x-%7.1f)/%6.1f)^2 ) range [%7.1f %7.1f]\n", i, f.a[i], f.a[i+1], f.a[i+2], bot, top);
	    }
	char fname[10];
	sprintf(fname,"pl%d", N);  // pl1, pl2, etc. for plot files
	FILE *fp = fopen(fname,"w");
	MeanStd m;
	for(int i=0; i<x.size(); i++) {
	    VecDoub DyDa(N*3);
	    double yy;
	    fgauss(x[i], f.a, yy, DyDa);
	    if (fp != NULL)  // if we cannot write, that's OK
		fprintf(fp,"%f %f %f\n", x[i], yy, y[i]);
	    residual[i] = y[i] - yy;
	    m.Element(residual[i]);
	    }
	double mean, std;
	m.Stats(mean, std);
	printf("Residuals: mean %f, RMS about mean %f\n", mean, std);
	a = f.a;
	if(fp != NULL) fclose(fp);
	}
    }
rmin = max(rmin, x[0]);  // rmin should not be smaller than the first x
rmax = min(rmax, x[x.size()-1]);  // and rmax not off the top end
}

// Note - if image is trimmed, values of w and h are modified to reflect trimmed image
uint8* NormalizeMRCImage(uint16* raster, uint32 &w, uint32 &h, int Ngauss)
{

int trim = 20;  // should make this a parameter
char *s = getenv("MRC_TRIM");
if (s != NULL) {
    trim = atoi(s);
    }
printf("Trimming %d pixels from each edge of MRC picture.\n", trim);
// trim the edges off the picture
// Also find mean and std deviation, both overall and per square
const int SQ = 8;
MeanStd Overall;
// create a 2D array of MeanStds.  Is there a cleverer way to do this?
vector<vector<MeanStd> > Checkerboard(SQ);
for(int i=0; i<SQ; i++)
   Checkerboard[i].resize(SQ);

int neww = w-2*trim;
int newh = h-2*trim;
int newp = neww*newh;
int maxval=0, minval=0x7FFFFFFF;
vector<double> vals(newp);
for(int x=0; x<neww; x++)
    for(int y=0; y<newh; y++) {
        int pix = raster[trim + x + (trim+y)*w];
        minval = min(minval, pix);
        maxval = max(maxval, pix);
	vals[x+neww*y] = pix;
        Overall.Element(pix);
        int sx = (x*SQ)/neww;
        int sy = (y*SQ)/newh;
        Checkerboard[sx][sy].Element(pix);
	}

w = neww;     // now modify values for caller
h = newh;

printf("Range of data values is [%d, %d], to be fit with %d gaussians\n", minval, maxval, Ngauss);
// compile a histogram; x values in 'mids', counts in 'histo'
vector<int> histo(256);
for(int i=0; i<neww*newh; i++) {
    int bin = int((vals[i]-minval)*256.0/(maxval+1-minval));
    histo[bin]++;
    }
vector<double>mids(256);
int nnz = 0;  // number of non-zero values
for(int i=0; i<histo.size(); i++) {
    mids[i] = (i+0.5)*(maxval+1-minval)/256.0 + minval;
    //printf("%f %d\n", mids[i], histo[i]);
    nnz += (histo[i] > 0);
    }
bool UseGaussians = true;
printf("Number of distinct pixel values: %d\n", nnz);

if (UseGaussians) {
    // Fit a N-gaussian (default 2) curve to the histogram
    double rmin, rmax;  // range min and maximum
    NGaussHist(mids, histo, Ngauss, rmin, rmax);
    printf("Convert using range [%.1f %.1f]\n", rmin, rmax);
    // Convert to 8 bit
    uint8* rast = (uint8*) malloc(newp*sizeof(uint8));
    for(int i=0; i<newp; i++) {
        int px = int((vals[i]-rmin)/(rmax-rmin)*256.0);  // stretch to 0..255
        if (px < 0) px = 0;
        if (px > 255) px = 255;
        rast[i] = px;
        }
    return rast;
    }
// Otherwise, try our older method
// Compare std of each square to the mean.
double mean, std;
Overall.Stats(mean, std);
printf("Overall mean=%f std=%f\n", mean, std);
double highest_contrast = 0.0;
double ss[SQ][SQ];  // saved std deviations
for(int y=0; y<SQ; y++) {
    for(int x=0; x<SQ; x++) {
        double m, s;
        Checkerboard[x][y].Stats(m,s);
        ss[x][y] = s;
        printf("%10.2f", s);
        highest_contrast = max(highest_contrast, s);
        }
    printf("\n");
    }

// use only the high constrast squares to compute image normalization
MeanStd ohc;  // ohc == Only high contrast
for(int x=0; x<neww; x++)
    for(int y=0; y < newh; y++) {
	int sx = (x*SQ)/neww;
        int sy = (y*SQ)/newh;
        if (ss[sx][sy] > highest_contrast/2.0)  // at least 1/4 of the of the highest contrast square
	    ohc.Element(vals[x+neww*y]);
	}
double ohc_mean, ohc_std;
ohc.Stats(ohc_mean, ohc_std);
printf("Using only high contrast squares with %d pixels, mean=%f, std=%f\n", ohc.HowMany(), ohc_mean, ohc_std);

// convert to a mean of 0 and a std dev of 1
for(int i=0; i<vals.size(); i++)
    vals[i] = (vals[i] - ohc_mean)/ohc_std;

// Convert to 8 bit
uint8* rast = (uint8*) malloc(newp*sizeof(uint8));
for(int i=0; i<newp; i++) {
    int px = int(vals[i]*42 + 127);  // make the limits +- 3 sigma
    if (px < 0) px = 0;
    if (px > 255) px = 255;
    rast[i] = px;
    }

return rast;
}

uint8* ReadAnMRCFile(const char *name, uint32 &neww, uint32 &newh, FILE *flog, bool Transpose, int Ngauss=2)
{
printf("Using %d gaussian normalization\n", Ngauss);
uint32 w,h;
vector<uint16*>raster =  ReadRawMRCFile(name, w, h, flog, Transpose);
if (raster.size() != 1)
    printf("Only reading first image of %d in MRC file\n", raster.size() );
neww = w; //since Normalize can change these if image is trimmed
newh = h;
uint8*result =  NormalizeMRCImage(raster[0], neww, newh, Ngauss);
free(raster[0]);
return result;
}

vector<uint8*> ReadMultiImageMRCFile(const char *name, uint32 &neww, uint32 &newh, FILE *flog, bool Transpose, int Ngauss=2)
{
uint32 w,h;
vector<uint16*>raster =  ReadRawMRCFile(name, w, h, flog, Transpose);
vector<uint8*> result;
for(int i=0; i<raster.size(); i++) {
    neww = w; //since Normalize can change these if image is trimmed
    newh = h;
    printf("\n----- Normalizing image %d\n", i);
    uint8* tmp = NormalizeMRCImage(raster[i],neww, newh, Ngauss);
    result.push_back(tmp);
    free(raster[i]);
    }
return result;
}
// ------------------ End of MRC file reading

