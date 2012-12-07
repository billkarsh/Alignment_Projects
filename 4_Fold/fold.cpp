

#include	"File.h"
#include	"ImageIO.h"
#include	"Maths.h"
#include	"Correlation.h"
#include	"Geometry.h"
#include	"CTForm.h"
#include	"Draw.h"






int main(int argc, char* argv[])
{
    uint32	w, h, w2, h2;
    size_t	npixels, npixels2;
    uint8*	raster = Raster8FromTif( argv[1], w, h );
    uint8*	raster2 = Raster8FromTif( argv[1], w2, h2 );

    npixels		= w * h;
    npixels2	= w2 * h2;

// Now draw the triangles onto the first file, write it out as s.tif
    uint8* r2 = (uint8*)RasterAlloc( npixels * sizeof(uint8) );
    int i;
    for(i=0; i<npixels; i++)
        r2[i] = raster[i];

    uint8* rt = (uint8*)RasterAlloc( npixels2 * sizeof(uint8) );
    for(i=0; i<npixels2; i++)
        rt[i] = raster2[i];

// Just testing
    FILE	*f = FileOpenOrDie( "fold.cmd", "r" );

    for(;;) {
        Point p1, p2, p3;
        double dx, dy;
        if (fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf",
         &(p1.x), &(p1.y), &(p2.x), &(p2.y), &(p3.x), &(p3.y), &dx, &dy) != 8)
            break;
        printf("Point %f %f %f %f %f %f %f %f\n", p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, dx, dy);
        vector<Point> tri;
        tri.push_back(p1);
        tri.push_back(p2);
        tri.push_back(p3);
        vector<Point> Plist;
        PixelListFromPolygon(Plist, tri);
        printf("Contains %d pixels\n", Plist.size() );
        vector<double> spv;  // source pixel values
        ValuesFromImageAndPoints(spv, raster, w, Plist);
        Normalize(spv);
        TForm t( 1, 0, dx, 0, 1, dy );

        // Now try to find a transform with a good match.
        for(double step=1; step > 0.05; ) {
	    vector<Point> Tpoints = Plist;   //   Start with locs in source
	    t.Transform( Tpoints );            // Transform to locations in target
	    vector<double> Tpixels;          // Get the pixel values
	    ValuesFromImageAndPoints(Tpixels, raster2, w, Tpoints);
	    Normalize(Tpixels);
            double start = CorrVectors(NULL, spv, Tpixels);
	    printf(" %f %f %f %f %f %f: correlation %f\n",
             t.t[0], t.t[1], t.t[2], t.t[3], t.t[4], t.t[5], start);
            double best_so_far = start;
            TForm tbest;  // best transform found
            int bdir = -1;  // flag to tell what the best direction is
            for(int dir=0; dir < 8; dir++) {
		double sx = step*cos(dir/8.0*2*PI);
		double sy = step*sin(dir/8.0*2*PI);
                TForm t2( t );
                t2.AddXY( sx, sy );
		vector<Point> Tpoints = Plist;   //   Start with locs in source
		t2.Transform( Tpoints );            // Transform to locations in target
		ValuesFromImageAndPoints(Tpixels, raster2, w, Tpoints);
		Normalize(Tpixels);
                double nc = CorrVectors(NULL, spv, Tpixels);
	        //printf(" case %d: correlation is %f\n", dir, nc);
                if (nc > best_so_far){
                    best_so_far = nc;
                    tbest.CopyIn( t2 );
                    bdir = dir;
                    }
                }
            Point cog;  // Center of gravity
	    cog = FindCOG(Plist);
            //printf("original GOG is %f %f\n", cog.x, cog.y);

            for(int rot = -1; rot<2; rot += 2) {  // try two rotations, too

                TForm	t2, R;
                CreateCWRot( R, rot * step, cog );
                MultiplyTrans( t2, R, t );

		vector<Point> Tpoints = Plist;   //   Start with locs in source
		t2.Transform( Tpoints );            // Transform to locations in target
		ValuesFromImageAndPoints(Tpixels, raster2, w, Tpoints);
		Normalize(Tpixels);
                double nc = CorrVectors(NULL, spv, Tpixels);
	        //printf(" rotated case %d: correlation is %f\n", rot, nc);
                if (nc > best_so_far){
                    best_so_far = nc;
                    bdir = 10;
                    tbest.CopyIn( t2 );
                    }
                }
            // Now tried 8 directions and two rotations; pick the best if any were better
            if( bdir >= 0 ) { //we found a better transform
                 t.CopyIn( tbest );
                 }
            else // nothing was better; reduce step size
                step = step/2;
            }
        // Now t is the best transform we can find.

	DrawLine(r2, w, h, p1.x, p1.y, p2.x, p2.y);
	DrawLine(r2, w, h, p2.x, p2.y, p3.x, p3.y);
	DrawLine(r2, w, h, p3.x, p3.y, p1.x, p1.y);
        t.Transform( tri );
        printf("Drew one\n");
        for(int i=0; i<3; i++) {
            int j = (i+1)%3;
	    DrawLine(rt, w2, h2, tri[i].x, tri[i].y, tri[j].x, tri[j].y);
            }
        printf("Drew two\n");
        }
    fclose(f);
   RasterFree(raster);
   RasterFree(raster2);

// Write results
	Raster8ToTif8( "s.tif", r2, w, h );
	Raster8ToTif8( "t.tif", rt, w2, h2 );

    return 0;
}
