/* --------------------------------------------------------------- */
/* ImproveControlPts --------------------------------------------- */
/* --------------------------------------------------------------- */

// Improve correlation by tweaking locations of control points.
//
// Return best correlation obtained.
//
// cpts			- the control points
// lambda		- points as lin. combs. of cpts
// spv			- values at those points
// image2		- target image mapped to
// w, h			- target image dims
// fsum			- log file
// describe		- string describing caller context
// threshold	- included in printed messages
//
double ImproveControlPts(
    vector<Point>					&cpts,
    const vector<vector<double> >	&lambda,
    vector<double>					&spv,
    const vector<double>			&image2,
    int								w,
    int								h,
    FILE							*fsum,
    const char						*describe,
    double							threshold )
{
    int		ncpnt	= cpts.size();

    fprintf( of, "Contains %d pixels\n", spv.size() );
    Normalize( spv );

    vector<double>	Tpixels;	// new target pixels
    vector<Point>	cderivs;	// control point derivatives

    GetBaryValues( Tpixels, cderivs, cpts, lambda, image2, w, h, spv );

    double	corr = CorrVectors( of, spv, Tpixels );

// initial plausibility check

    if( corr < DetailInitialThreshold ) {

        fprintf( of,
        "STAT:Correlation less than %f at the start of iterations,"
        " was %f\n", DetailInitialThreshold, corr );

        fprintf( fsum,
        "STAT:Correlation less than %f at the start of iterations,"
        " was %f\n", DetailInitialThreshold, corr );

        fprintf( of, "Control points are:" );

        for( int i = 0; i < ncpnt; ++i )
            fprintf( of, "(%.3f %.3f) ", cpts[i].x, cpts[i].y );

        fprintf( of, "corr=%f\n", corr );

#if 0
        for( int i = 0; i < lambda.size(); ++i ) {

            fprintf( of, "---i=%d\n", i );

            for( int j = 0; j < lambda[i].size(); ++j )
                fprintf( of, "%.4f ", lambda[i][j] );

            fprintf( of, "\n spv=%f Tpixels=%f\n", spv[i], Tpixels[i] );
        }
#endif

        return 0.0;
    }

    fprintf( of,
    "STAT: Initial %s correlation %f\n", describe, corr );

    fprintf( fsum,
    "STAT: Initial %s correlation %f\n", describe, corr );

// Try to tweak the control points for a good match

    for( double step = 10.0; step > 0.05; ) {

        fprintf( of, "\nControl points are:" );

        for( int i = 0; i < ncpnt; ++i )
            fprintf( of, "(%f %f) ", cpts[i].x, cpts[i].y );

        fprintf( of, "corr=%f, step is %f\n", corr, step );

        // compute a downhill vector of length 'step'.

        double	S = 0.0;

        for( int i = 0; i < ncpnt; ++i )
            S += cderivs[i].RSqr();

        if( !S ) {
            fprintf( of, "*** ALL DERIVS ZERO.\n" );
            break;
        }

        S = step / sqrt( S );

        // find the new (and hopefully improved) control points

        vector<Point>	newpts = cpts;

        for( int i = 0; i < ncpnt; ++i ) {
            newpts[i].x += S * cderivs[i].x;
            newpts[i].y += S * cderivs[i].y;
        }

        // is the new spot better?

        vector<Point>	newderivs;

        GetBaryValues( Tpixels, newderivs, newpts, lambda,
            image2, w, h, spv );

        double c = CorrVectors( of, spv, Tpixels );

        if( c > corr ) {

            corr	= c;
            cpts	= newpts;
            cderivs	= newderivs;
        }
        else
            step /= 2;
    }

    fprintf( of,
    "STAT: Final %s correlation %f, (threshold %f)\n",
    describe, corr, threshold );

    fprintf( fsum,
    "STAT: Final %s correlation %f, (threshold %f)\n",
    describe, corr, threshold );

    return corr;
}

