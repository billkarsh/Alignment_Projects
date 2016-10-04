

#ifdef USE_CURL

#include	"janelia.h"

#include	<curl/curl.h>
#include	<stdlib.h>
#include	<string.h>


/* --------------------------------------------------------------- */
/* Discussion ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// The Janelia FlyTEM pipeline uses JSON tile descriptors like this:
//
//{
//  "tileId" : "150311155421069028.3333.0",
//  "layout" : {
//    "sectionId" : "3333.0",
//    "temca" : "7",
//    "camera" : "0",
//    "imageRow" : 28,
//    "imageCol" : 69,
//    "stageX" : 129450.0,
//    "stageY" : 51100.0,
//    "rotation" : 180.0
//  },
//  "z" : 3333.0,
//  "minX" : 129562.0,
//  "minY" : 51127.0,
//  "maxX" : 132205.0,
//  "maxY" : 53423.0,
//  "width" : 2560.0,
//  "height" : 2160.0,
//  "minIntensity" : 0.0,
//  "maxIntensity" : 255.0,
//  "mipmapLevels" : {
//    "0" : {
//      "imageUrl" : "file:/tier2/flyTEM/eric/data/whole_fly_1/150311155421_132x216/col0069/col0069_row0028_cam0.tif",
//      "maskUrl" : "file:/tier2/flyTEM/eric/working_sets/150630_fafb_V8_delta/temp_render_masks/temca7_cam0.png"
//    }
//  },
//  "transforms" : {
//    "type" : "list",
//    "specList" : [ {
//      "type" : "ref",
//      "refId" : "150630offset_temca7_camera0"
//    }, {
//      "type" : "leaf",
//      "className" : "mpicbg.trakem2.transform.AffineModel2D",
//      "dataString" : "1.0 0.0 0.0 1.0 129450 51100"
//    } ]
//  },
//  "meshCellSize" : 64.0
//}
//
// These data map to a PipeFiles::PicSpec structure via commandline options
// -jtilea=URL -jtileb=URL
//

/* --------------------------------------------------------------- */
/* stripQuotes --------------------------------------------------- */
/* --------------------------------------------------------------- */

static string stripQuotes( string &s )
{
    if( s.find( "\"" ) == 0 )
        return s.substr( 1, s.length() - 2 );
    else
        return s;
}

/* --------------------------------------------------------------- */
/* copy_img_from_URL --------------------------------------------- */
/* --------------------------------------------------------------- */

// Given path to image file of form: "file:/path/file.ext",
// copy that file to /tmp/<uniquename>.ext.
//
// Return true if ok.
//
// http://www.labbookpages.co.uk/software/imgProc/libPNG.html
// http://stackoverflow.com/questions/12728524/save-an-image-from-jpeg-to-png-using-libcurl-and-gdkpixbuff
//
static bool copy_img_from_URL( string &output_path, string url )
{
    CURL	*easy_handle = curl_easy_init();
    bool	ok = false;

    if( easy_handle ) {

        string	suffix = url.substr( url.find_last_of(".") );

        output_path	= tmpnam( NULL ) + suffix;

        // Open file

        FILE	*fp = fopen( output_path.c_str(), "wb" );

        if( !fp ) {
            (void)perror( "Extract: The following error occurred: " );
            goto exit;
        }

        curl_easy_setopt( easy_handle, CURLOPT_URL, url.c_str() );
        curl_easy_setopt( easy_handle, CURLOPT_WRITEFUNCTION, NULL );
        curl_easy_setopt( easy_handle, CURLOPT_WRITEDATA, fp );

        // Grab image

        CURLcode	imgresult = curl_easy_perform( easy_handle );

        if( imgresult )
            fprintf( stderr, "Extract: Cannot grab the image.\n" );
        else
            ok = true;

        fclose( fp );
    }

exit:
    curl_easy_cleanup( easy_handle );
    return ok;
}

/* --------------------------------------------------------------- */
/* url_to_path --------------------------------------------------- */
/* --------------------------------------------------------------- */

// Remove 'file:' prefix from path.
//
static string url_to_path( string &url )
{
    string	prefix = "file:";

    if( url.find( prefix ) == 0 )
        return url.substr( prefix.length() );
    else
        return url;
}

/* --------------------------------------------------------------- */
/* GetArgListFromURL --------------------------------------------- */
/* --------------------------------------------------------------- */

static void Tokenize(
    vector<string>	&tokens,
    const string	&str,
    const string	&delimiters,
    bool			stripQ )
{
    string::size_type	start, stop = 0;

    for(;;) {

        start = str.find_first_not_of( delimiters, stop );

        if( start == string::npos )
            return;

        stop = str.find_first_of( delimiters, start + 1 );

        string	tok = str.substr( start, stop - start );

        if( stripQ )
            tokens.push_back( stripQuotes( tok ) );
        else
            tokens.push_back( tok );
    }
}


static size_t WriteCallback( void *src, size_t size, size_t nmemb, void *dst )
{
    size_t	bytes = size * nmemb;

    ((std::string*)dst)->append( (char*)src, bytes );

    return bytes;
}


#define	GETDOUBLE( field, name )										\
    if( !tokens2[0].compare( name ) )									\
        {field = strtod( tokens2[2].c_str(), NULL );}

#define	GETSTRING( field, name )										\
    if( !tokens2[0].compare( name ) )									\
        {field = tokens2[2];}


// URL examples:
// http://...
// file:///groups...
// file://localhost/groups...
//
// Simple parsing hack identifies colon-separated JSON pairs:
// - token[0] = name
// - token[1] = ':'
// - token[2] = value
//
bool GetTileSpecFromURL( PicSpec &P, const char *pat, char *argv )
{
    int	len = strlen( pat );

    if( !strncmp( argv, pat, len ) ) {

        std::string	readBuffer;
        CURL		*easy_handle = curl_easy_init();
        CURLcode	res;

//		curl_easy_setopt( easy_handle, CURLOPT_VERBOSE, 1L );
        curl_easy_setopt( easy_handle, CURLOPT_URL, argv + len );
        curl_easy_setopt( easy_handle, CURLOPT_WRITEFUNCTION, WriteCallback );
        curl_easy_setopt( easy_handle, CURLOPT_WRITEDATA, &readBuffer );

        res = curl_easy_perform( easy_handle );
        curl_easy_cleanup( easy_handle );

        if( res ) {
            fprintf( stderr, "jtile: Can't parse [%s].\n", argv );
            return false;
        }

        vector<string>	tokens;
        Tokenize( tokens, readBuffer, "\n", false );

        for( int i = 0, n = tokens.size(); i < n; ++i )  {

            vector<string>	tokens2;
            Tokenize( tokens2, tokens[i], " ,\n", true );

            if( tokens2.size() != 3 )
                continue;

            GETDOUBLE( P.z, "z" );
            GETSTRING( P.t2i.path, "imageUrl" );
            GETDOUBLE( P.t2i.T.t[2], "stageX" );
            GETDOUBLE( P.t2i.T.t[5], "stageY" );
            GETDOUBLE( P.t2i.col, "imageCol" );
            GETDOUBLE( P.t2i.row, "imageRow" );
            GETDOUBLE( P.t2i.cam, "camera" );
        }

//		copy_img_from_URL( P.t2i.path, P.t2i.path );
        P.t2i.path	= url_to_path( P.t2i.path );
        P.t2i.id	= -1;

        if( P.t2i.path.size() > 0 ) {

            P.id = -1;
            return true;
        }
    }

    return false;
}

#endif	// USE_CURL


