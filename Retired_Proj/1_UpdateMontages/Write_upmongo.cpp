/* --------------------------------------------------------------- */
/* Write_upmongo --------------------------------------------------- */
/* --------------------------------------------------------------- */

static void Write_upmongo()
{
	char	buf[2048], inname[256];
	FILE	*f;

	if( FileIsExt( gArgs.infile, ".xml" ) && !XMLHasMRC() ) {

		const char	*name	= FileNamePtr( gArgs.infile ),
					*dot	= FileDotPtr( name );

		sprintf( inname, "%.*s", int(dot - name), name );
	}
	else
		strcpy( inname, "PreClicks" );

	sprintf( buf, "upmongo.sht" );
	f = FileOpenOrDie( buf, "w", flog );

	fprintf( f, "#!/bin/sh\n" );
	fprintf( f, "\n" );
	fprintf( f, "# Update transforms in given 'myxml' file using lyr/montage results\n" );
	fprintf( f, "# in given temp directory. Preserves montage-montage orientation.\n" );
	fprintf( f, "# Output file named 'myxml_v2.xml'.\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# > updatemontages myxml temp\n" );
	fprintf( f, "#\n" );
	fprintf( f, "# Options:\n" );
	fprintf( f, "# -zmin=i -zmax=j\t;restricts layer range\n" );
	fprintf( f, "# -force\t\t\t;overwrite TForms with lsq solutions\n" );
	fprintf( f, "\n" );
	fprintf( f, "\n" );
	fprintf( f, "inname=%s\n",
	inname );
	fprintf( f, "\n" );
	fprintf( f, "updatemontages $inname.xml temp0 -zmin=%d -zmax=%d\n",
	gArgs.zmin, gArgs.zmax );
	fprintf( f, "\n" );
	fprintf( f, "mv $inname\"_v2.xml\" \"newmons.xml\"\n" );
	fprintf( f, "\n" );

	fclose( f );
	FileScriptPerms( buf );
}


