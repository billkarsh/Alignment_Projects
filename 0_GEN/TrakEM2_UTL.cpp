

#include	"File.h"
#include	"TrakEM2_UTL.h"

#include	<string.h>






/* --------------------------------------------------------------- */
/* class XML_TKEM ------------------------------------------------ */
/* --------------------------------------------------------------- */

void XML_TKEM::Open( const char *file, FILE* flog )
{
    this->flog	= flog;
    this->file	= file;

    if( !doc.LoadFile( file ) ) {
        fprintf( flog, "Can't open XML file [%s].\n", file );
        exit( 42 );
    }

    if( !doc.FirstChild( "trakem2" ) ) {
        fprintf( flog, "No <trakEM2> tag [%s].\n", file );
        exit( 42 );
    }
}


void XML_TKEM::Save( const char *name, bool copyDTD )
{
    doc.SaveFile( name );

    if( copyDTD )
        CopyDTD( file, name );
}


TiXmlNode* XML_TKEM::GetLayerset()
{
    TiXmlHandle	hdoc( &doc );
    TiXmlNode*	layerset;

    layerset = hdoc.FirstChild( "trakem2" )
                .FirstChild( "t2_layer_set" )
                .ToNode();

    if( !layerset ) {
        fprintf( flog, "No <t2_layer_set> tag [%s].\n", file );
        exit( 42 );
    }

    return layerset;
}


TiXmlElement* XML_TKEM::GetFirstLayer()
{
    TiXmlHandle		hdoc( &doc );
    TiXmlElement*	layer;

    layer = hdoc.FirstChild( "trakem2" )
                .FirstChild( "t2_layer_set" )
                .FirstChild( "t2_layer" )
                .ToElement();

    if( !layer ) {
        fprintf( flog, "No <t2_layer> tag [%s].\n", file );
        exit( 42 );
    }

    return layer;
}


TiXmlElement* XML_TKEM::GetLastLayer()
{
    TiXmlHandle		hdoc( &doc );
    TiXmlElement*	layer;

    layer = hdoc.FirstChild( "trakem2" )
                .FirstChild( "t2_layer_set" )
                .ToNode()
                ->LastChild( "t2_layer" )
                ->ToElement();

    if( !layer ) {
        fprintf( flog, "No <t2_layer> tag [%s].\n", file );
        exit( 42 );
    }

    return layer;
}


int XML_TKEM::NextOID()
{
    return ::NextOID( &doc );
}

/* -------------------------------------------------------------- */
/* IDFromPatch -------------------------------------------------- */
/* -------------------------------------------------------------- */

int IDFromPatch( TiXmlElement* ptch )
{
    int	z, id;

    sscanf( ptch->Attribute( "title" ), "%d.%d", &z, &id );

    return id;
}

/* --------------------------------------------------------------- */
/* XMLSetTFVals -------------------------------------------------- */
/* --------------------------------------------------------------- */

void XMLSetTFVals( TiXmlElement* ptch, const double *t )
{
    char	buf[256];

    sprintf( buf, "matrix(%f,%f,%f,%f,%f,%f)",
    t[0], t[3], t[1], t[4], t[2], t[5] );

    ptch->SetAttribute( "transform", buf );
}

/* --------------------------------------------------------------- */
/* NextOID ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Recursively search all children of par node for either
// 'id' or 'oid' tags.
//
// Return highest found.
//
static int NextOID_IntoNode( const TiXmlNode* par )
{
    int	highest = 0;

    const TiXmlNode*	ch = par->FirstChild();

    for( ; ch; ch = par->IterateChildren( ch ) ) {

        const char*	s;
        int			id;

        s = ((TiXmlElement*)ch)->Attribute( "oid" );

        if( s ) {

            id = atoi( s );

            if( id > highest )
                highest = id;
        }

        s = ((TiXmlElement*)ch)->Attribute( "id" );

        if( s ) {

            id = atoi( s );

            if( id > highest )
                highest = id;
        }

        id = NextOID_IntoNode( ch );

        if( id > highest )
            highest = id;
    }

    return highest;
}


// Recursively search all children of "trakem2" node for either
// 'id' or 'oid' tags.
//
// Return highest found + 1.
//
int NextOID( const TiXmlHandle hDoc )
{
    int	highest = 0;

    TiXmlNode*	par = hDoc.FirstChild( "trakem2" ).ToNode();

    if( par )
        highest = NextOID_IntoNode( par );

    return highest + 1;
}

/* --------------------------------------------------------------- */
/* SetOID -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// If par object has an 'oid' set it to nextoid and recursively
// do the same for all children of par.
//
// Return highest applied oid + 1.
//
int SetOID( TiXmlNode* par, int nextoid )
{
    const char*	s;

    s = ((TiXmlElement*)par)->Attribute( "oid" );

    if( s )
        ((TiXmlElement*)par)->SetAttribute( "oid", nextoid++ );

    TiXmlNode*	ch = par->FirstChild();

    for( ; ch; ch = par->IterateChildren( ch ) )
        nextoid = SetOID( ch, nextoid );

    return nextoid;
}

/* --------------------------------------------------------------- */
/* CopyDTD ------------------------------------------------------- */
/* --------------------------------------------------------------- */

// TrakEM2 requires input xml files to begin with DTD <!DOCTYPE>
// block, but TinyXML ignores DTD completely. Therefore, we copy
// those data from the original input file.
//
void CopyDTD( const char *dtdsrc, const char *bodysrc )
{
    CLineScan	LS;
    FILE		*fi, *fo;
    int			len;

// Open output file
    LS.bufsize	= 2048;
    LS.line		= (char*)malloc( LS.bufsize );

    len = sprintf( LS.line, "%s", dtdsrc );
    strcpy( LS.line + len - 4, "_v2.xml" );

    if( fo = fopen( LS.line, "w" ) ) {

        // Get DTD from input file
        if( fi = fopen( dtdsrc, "r" ) ) {

            while( LS.Get( fi ) > 0 ) {

                fprintf( fo, LS.line );

                if( !strncmp( LS.line, "] >", 3 ) ) {
                    fprintf( fo, "\n" );
                    break;
                }
            }

            fclose( fi );
        }

        // Get all but first header line from bodysrc
        if( fi = fopen( bodysrc, "r" ) ) {

            LS.Get( fi );

            while( LS.Get( fi ) > 0 )
                fprintf( fo, LS.line );

            fclose( fi );
        }

        fclose( fo );
    }

    remove( bodysrc );
}

/* --------------------------------------------------------------- */
/* TrakEM2WriteDTD ----------------------------------------------- */
/* --------------------------------------------------------------- */

// Required TrackEM2 document type data.
//
void TrakEM2WriteDTD( FILE *f )
{
    fprintf( f, "<!DOCTYPE trakem2_anything [\n" );
    fprintf( f, "	<!ELEMENT trakem2 (project,t2_layer_set,t2_display)>\n" );
    fprintf( f, "	<!ELEMENT project (anything)>\n" );
    fprintf( f, "	<!ATTLIST project id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project preprocessor NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project mipmaps_folder NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project storage_folder NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT anything EMPTY>\n" );
    fprintf( f, "	<!ATTLIST anything id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST anything expanded NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_layer (t2_patch,t2_label,t2_layer_set,t2_profile)>\n" );
    fprintf( f, "	<!ATTLIST t2_layer oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer thickness NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer z NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_layer_set (t2_layer,t2_pipe,t2_ball,t2_area_list,t2_calibration)>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_height NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_z NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set snapshots_quality NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set snapshots_mode NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_calibration EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelWidth NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelHeight NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelDepth NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration xOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration yOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration zOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration info NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration valueUnit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration timeUnit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration unit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_ball (t2_ball_ob)>\n" );
    fprintf( f, "	<!ATTLIST t2_ball oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_ball_ob EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob r NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_label EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_label oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_patch EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_patch oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch file_path NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch original_path NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch type NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_pipe EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe p_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_profile EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_profile oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_area_list (t2_area)>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list fill_paint NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_area (t2_path)>\n" );
    fprintf( f, "	<!ATTLIST t2_area layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_path EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_path d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_dissector (t2_dd_item)>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_dd_item EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item radius NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item tag NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item points NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_display EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_display id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display magnification NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_height NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display scroll_step NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display c_alphas NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display c_alphas_state NMTOKEN #REQUIRED>\n" );
    fprintf( f, "] >\n\n" );
}

/* --------------------------------------------------------------- */
/* TrakEM2WriteDTDEx --------------------------------------------- */
/* --------------------------------------------------------------- */

// Required TrackEM2 document type data.
//
void TrakEM2WriteDTDEx( FILE *f )
{
    fprintf( f, "<!DOCTYPE trakem2_anything [\n" );
    fprintf( f, "	<!ELEMENT trakem2 (project,t2_layer_set,t2_display)>\n" );
    fprintf( f, "	<!ELEMENT project (anything)>\n" );
    fprintf( f, "	<!ATTLIST project id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project unuid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project preprocessor NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project mipmaps_folder NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST project storage_folder NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT anything EMPTY>\n" );
    fprintf( f, "	<!ATTLIST anything id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST anything expanded NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_layer (t2_patch,t2_label,t2_layer_set,t2_profile)>\n" );
    fprintf( f, "	<!ATTLIST t2_layer oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer thickness NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer z NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_layer_set (t2_prop,t2_linked_prop,t2_annot,t2_layer,t2_pipe,t2_ball,t2_area_list,t2_calibration,t2_stack,t2_treeline)>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set layer_height NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set rot_z NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set snapshots_quality NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set color_cues NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set n_layers_color_cue NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set paint_arrows NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_layer_set paint_edge_confidence_boxes NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_calibration EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelWidth NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelHeight NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration pixelDepth NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration xOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration yOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration zOrigin NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration info NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration valueUnit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration timeUnit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_calibration unit NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_ball (t2_prop,t2_linked_prop,t2_annot,t2_ball_ob)>\n" );
    fprintf( f, "	<!ATTLIST t2_ball oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball fill NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_ball_ob EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob r NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_ball_ob layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_label (t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_label oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_label composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_patch (t2_prop,t2_linked_prop,t2_annot,ict_transform,ict_transform_list)>\n" );
    fprintf( f, "	<!ATTLIST t2_patch oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch file_path NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch original_path NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch type NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch ct NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch o_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch o_height NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_patch pps NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_pipe (t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe p_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_pipe layer_ids NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_polyline (t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_polyline d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_profile (t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_profile oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_profile d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_area_list (t2_prop,t2_linked_prop,t2_annot,t2_area)>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_area_list fill_paint NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_area (t2_path)>\n" );
    fprintf( f, "	<!ATTLIST t2_area layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_path EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_path d NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_dissector (t2_prop,t2_linked_prop,t2_annot,t2_dd_item)>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dissector composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_dd_item EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item radius NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item tag NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_dd_item points NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_stack (t2_prop,t2_linked_prop,t2_annot,(iict_transform|iict_transform_list)?)>\n" );
    fprintf( f, "	<!ATTLIST t2_stack oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack file_path CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_stack depth CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_tag EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_tag name NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_tag key NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_node (t2_area*,t2_tag*)>\n" );
    fprintf( f, "	<!ATTLIST t2_node x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_node y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_node lid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_node c NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_node r NMTOKEN #IMPLIED>\n" );
    fprintf( f, "	<!ELEMENT t2_treeline (t2_node*,t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_treeline composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_areatree (t2_node*,t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_areatree composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_connector (t2_node*,t2_prop,t2_linked_prop,t2_annot)>\n" );
    fprintf( f, "	<!ATTLIST t2_connector oid NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector transform NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector style NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector locked NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector visible NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector title NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector links NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_connector composite NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_prop EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_prop key NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_prop value NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_linked_prop EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_linked_prop target_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_linked_prop key NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_linked_prop value NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT t2_annot EMPTY>\n" );
    fprintf( f, "	<!ELEMENT t2_display EMPTY>\n" );
    fprintf( f, "	<!ATTLIST t2_display id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display layer_id NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display magnification NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_x NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_y NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_width NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display srcrect_height NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display scroll_step NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display c_alphas NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST t2_display c_alphas_state NMTOKEN #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT ict_transform EMPTY>\n" );
    fprintf( f, "	<!ATTLIST ict_transform class CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST ict_transform data CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT iict_transform EMPTY>\n" );
    fprintf( f, "	<!ATTLIST iict_transform class CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ATTLIST iict_transform data CDATA #REQUIRED>\n" );
    fprintf( f, "	<!ELEMENT ict_transform_list (ict_transform|iict_transform)*>\n" );
    fprintf( f, "	<!ELEMENT iict_transform_list (iict_transform*)>\n" );
    fprintf( f, "] >\n\n" );
}


