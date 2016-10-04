// BK - Unused fun to read TrackEM file.
// Identical versions in two.cpp and two32raster.cpp.
//


// Read in a track_em file
//
int ReadTrackEM(const char *name, FILE *flog)
{
// Read the trakEM xml file to see what's there.
vector<Picture> vp;
TiXmlDocument doc(name);
bool loadOK = doc.LoadFile();
printf("XML load gives %d\n", loadOK);
if (!loadOK)  {
    printf("Could not open XML file '%s'\n", name);
    fprintf(flog,"Could not open XML file '%s'\n", name);
    exit( 42 );
    }
TiXmlHandle hDoc(&doc);
TiXmlElement* child;
TiXmlHandle hRoot(0);

// block: should be <trakem2>
TiXmlNode*	node=0;
node =doc.FirstChild();
if( node == NULL )
    printf("No node??\n");
child=hDoc.FirstChild("trakem2").FirstChild("t2_layer_set").FirstChild("t2_layer").ToElement();
// should always have a valid root but handle gracefully if it does
if( !child ) {
    printf("No first child??\n");
    return 42;
    }
printf("child element value %s\n", child->Value() );
for( child; child; child=child->NextSiblingElement() ) {
    //printf("Got a <t2_layer>\n");
    const char *attr  = child->Attribute("z");
    //printf("z= %s\n", attr);
    // now look through all the <t2_patch> elements in each
    TiXmlElement*	c2;
    c2 = child->FirstChildElement("t2_patch");
    for( c2; c2; c2=c2->NextSiblingElement() ) {
    //printf("got a <t2_patch>\n");
        const char *tf  = c2->Attribute("transform");
        const char *fp  = c2->Attribute("file_path");
        //printf("File is '%s'\n, transform is '%s'\n", fp, tf);
        Picture p;
        p.tr.ScanTrackEM2( tf );
        p.z = int(atof(attr)+0.5);
        p.fname = fp;
        p.w = int(atof(c2->Attribute("width")));
        p.h = int(atof(c2->Attribute("height")));
        vp.push_back(p);
    }
    }
return 0;
}


