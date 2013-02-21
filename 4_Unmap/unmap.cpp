

#include	"File.h"
#include	"ImageIO.h"
#include	"TAffine.h"


// structure for unmapping a global image


struct orig_image {
  public:
    char *fname;
    int w,h;
    int xmin, ymin, xmax, ymax; // bounding box in global image of all point that map to this
    orig_image(){fname = NULL; xmin = ymin = 1000000000; xmax = ymax = -1000000000;}
    orig_image(char *f){fname = f; xmin = ymin = 1000000000; xmax = ymax = -1000000000;}
    };

struct one_tform {
  public:
    int image_id;  // which image
    TAffine tr;      // maps from global space to individual image
  };






int main(int argc, char **argv)
{
vector<char *>noa;  // non-option arguments
for(int i=1; i<argc; i++) {
    // process arguments here
    if( argv[i][0] != '-' )
	noa.push_back(argv[i]);
    else
	printf("Ignored option '%s'\n", argv[i]);
    }
if( noa.size() < 3 ) {
    printf("Usage: unmap <image file> <map file> <file-of-transforms> [<where-to-put>] \n");
    exit( 42 );
    }
// step 1 - read the image file
uint32 w,h;
uint16* raster = Raster16FromPng(noa[0], w, h);
printf("width is %d, height %d\n", w, h);

// step 2 - read the mapping file
uint32 wm,hm;
uint16* map = Raster16FromPng(noa[1], wm, hm);
printf("width of map is %d, height %d\n", wm, hm);

// Step 3 - read the file of images and transforms
vector<orig_image> images;
vector<one_tform> tforms(1);
FILE *fp = FileOpenOrDie( noa[2], "r" );

{
	CLineScan	LS;

	for(;;) {
		if( LS.Get( fp ) <= 0 )
			break;
		if(strncmp(LS.line,"IMAGE",5) == 0) {
			int id;
			char fname[2048];
			sscanf(LS.line+5, "%d %s", &id, fname);
			printf("id %3d name %s\n", id, fname);
			if( id != images.size() ) {
			printf("Oops - bad image sequence number %d\n", id);
			return 42;
			}
			images.push_back(orig_image(strdup(strtok(fname," '\n"))));
		}
		else if(strncmp(LS.line,"TRANS",5) == 0) {
			int id, image_no;
			double a,b,c,d,e,f;
			sscanf(LS.line+5,"%d %d %lf %lf %lf %lf %lf %lf", &id, &image_no, &a, &b, &c, &d, &e, &f);
			if( id != tforms.size() ) {
			printf("Oops - bad transform sequence number %d\n", id);
			return 42;
			}
			one_tform o;
			o.image_id = image_no;
			o.tr = TAffine(a,b,c,d,e,f);
			tforms.push_back(o);
		}
		else {
			printf("UNknown line %s\n", LS.line);
		}
	}
}

fclose(fp);

// OK, find the bounding bozes for all the images
for(int y=0; y<h; y++) {
    for(int x=0; x<w; x++) {
	uint16 t = map[x + w*y];
        if( t != 0 ) {
	    int im = tforms[t].image_id;
            images[im].xmin = min(images[im].xmin, x);
            images[im].ymin = min(images[im].ymin, y);
            images[im].xmax = max(images[im].xmax, x);
            images[im].ymax = max(images[im].ymax, y);
	    }
        }
    }
//Now compute each image one at a time
for(int i=0; i<images.size(); i++) {
    int x0 = images[i].xmin;
    int x1 = images[i].xmax;
    int y0 = images[i].ymin;
    int y1 = images[i].ymax;
    printf("Image %d, x=[%6d %6d] y = [%6d %6d]\n", i, x0, x1, y0, y1);
    uint32 w, h;
    //uint8* junk = Raster8FromTif( images[i].fname, w, h );
    //RasterFree(junk);
    vector<uint8> recon_raster(w*h,0);  // create an empty raster
    printf("Original size was %d wide by %d tall\n", w, h);
    for(int y=y0; y<=y1; y++) {
        for(int x=x0; x<=x1; x++) {
	    uint16 t = map[x + wm*y];
            if( t != 0 && tforms[t].image_id == i ) {  // maps to image i
		Point pt(x,y);
		tforms[t].tr.Transform( pt );
                if( x == 8557 && y == 431 ) {  // just for debugging
                    printf("X and Y in original image: %d %d.  Pixel value is %d\n", x, y, t);
                    printf("Image id is %d. Transformation is", tforms[t].image_id);
		    tforms[t].tr.TPrint();
                    printf("Point in image: x=%f y=%f\n", pt.x, pt.y);
		    }
		// This should be within the image, but double check
                if( pt.x > -0.5 && pt.x < w-0.5 && pt.y > -0.5 && pt.y < h-0.5 ) { // it will round to a legal value
		    int ix = int(pt.x+0.5);
                    int iy = int(pt.y+0.5);
                    recon_raster[ix + w*iy] = raster[x + wm*y];
		    }
		else {
                    printf("X and Y in original image: %d %d.  Pixel value is %d\n", x, y, t);
                    printf("Image id is %d. Transformation is", tforms[t].image_id);
		    tforms[t].tr.TPrint();
                    printf("Point out of image: x=%f y=%f\n", pt.x, pt.y);
                    //return 42;
		    }
                }
	    }
	}
    char fname[256];
    sprintf(fname,"/tmp/%d.png", i);
    Raster8ToPng8(fname, &recon_raster[0], w, h);
    }
return 0;
}
