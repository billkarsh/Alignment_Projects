

#include	"lsq_Msg.h"

#include	"Disk.h"
#include	"File.h"
#include	"Timer.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<unistd.h>

#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






void MsgClear()
{
	system( "rm -rf MSG/*" );
}


void MsgSend( const char *msg )
{
	MsgClear();

	FILE	*f = fopen( "MSG/msg", "w" );

	if( f ) {
		fprintf( f, msg );
		fclose( f );
	}
}


bool MsgWaitAck( int nw, int timeout_ms )
{
	clock_t		t0 = StartTiming();
	vector<int>	ack( nw+1, 0 );
	int			sleep_ms = timeout_ms / 10;

	for(;;) {

		int	nack = 0;

		for( int iw = 1; iw <= nw; ++iw ) {

			if( ack[iw] )
				++nack;
			else {
				char	buf[32];
				sprintf( buf, "MSG/ack_%d", iw );

				if( DskExists( buf ) ) {
					ack[iw] = 1;
					++nack;
				}
			}
		}

		if( nack >= nw )
			return true;
		else if( DeltaSeconds( t0 ) * 1000 > timeout_ms )
			break;

		usleep( sleep_ms*1000 );	// microsec
	}

	return false;
}


bool MsgWaitRecv( char *msgbuf, int timeout_ms )
{
	clock_t	t0 = StartTiming();
	int		sleep_ms = timeout_ms / 10;

	for(;;) {

		if( DskExists( "MSG/msg" ) ) {

			FILE	*f = fopen( "MSG/msg", "r" );

			if( f ) {
				CLineScan	LS;
				LS.Get( f );
				strcpy( msgbuf, LS.line );
				fclose( f );
			}

			return true;
		}
		else if( DeltaSeconds( t0 ) * 1000 > timeout_ms )
			break;

		usleep( sleep_ms*1000 );	// microsec
	}

	return false;
}


void MsgAck( int iw )
{
	char	buf[32];
	sprintf( buf, "MSG/ack_%d", iw );

	FILE	*f = fopen( buf, "w" );

	if( f )
		fclose( f );
}


