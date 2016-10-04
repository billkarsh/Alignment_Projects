/* --------------------------------------------------------------- */
/* FindConnRegions ----------------------------------------------- */
/* --------------------------------------------------------------- */

// This is Lou's original implementation.

void FindConnRegions(vector<ConnRegion> &cr, uint8* fold_mask, int w, int h, int scale)
{
fprintf(of,"Finding connected regions: w=%d h=%d scale=%d\n", w, h, scale);
int n = w * h;
int big = 0;
for(int i=0; i<n; i++)
    big = max(big, fold_mask[i]);
cr.clear();
cr.insert(cr.begin(), big, ConnRegion());
for(int i=1; i <= big; i++) {
    cr[i-1].id = i;
    vector<vertex> vv;
    for(int j=0; j<n; j++) {
        if( fold_mask[j] == i ) {
            int iy = j / w;
            int ix = j - w * iy;
            vv.push_back( vertex( ix/scale, iy/scale ) );
        }
    }
    fprintf(of, "%d pixels in fold mask had value %d\n", vv.size(), i);
    // may contain duplicates from down-scaling.  Sort, then take unique values.
    sort(vv.begin(), vv.end());
    for(int j=0; j<vv.size(); j++) {
    if( j == 0 || !(vv[j] == vv[j-1]) ) {
        cr[i-1].pts.push_back(Point(vv[j].x, vv[j].y));
            cr[i-1].B.L = min(cr[i-1].B.L, vv[j].x);
            cr[i-1].B.R = max(cr[i-1].B.R, vv[j].x);
            cr[i-1].B.B = min(cr[i-1].B.B, vv[j].y);
            cr[i-1].B.T = max(cr[i-1].B.T, vv[j].y);
        }
    }
    fprintf(of, "CR %d, %d points, (%d %d) to (%d %d)\n", i, cr[i-1].pts.size(), cr[i-1].B.L, cr[i-1].B.B, cr[i-1].B.R, cr[i-1].B.T);
    }
// drop any that do not have enough points
for(int i=0; i<cr.size(); i++) {
    if( cr[i].pts.size() < 0.9 * MinMapArea ) {
    cr.erase(cr.begin()+i); // delete the entry
    i--;  // will need to look at this index again
    }
    }
fprintf(of, "%d CRs were big enough\n", cr.size() );
}


