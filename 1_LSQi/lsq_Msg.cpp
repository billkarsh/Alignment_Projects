

#include	"lsq_Msg.h"

#include	"Disk.h"
#include	"Timer.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<unistd.h>

#include	<vector>
using namespace std;


/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






void MsgSyncWorkers( int iw, int nw )
{
	if( nw <= 1 )
		return;

	clock_t	t0 = StartTiming();

	const int jobstart_s = 20*60;	// 20 min

	if( !iw ) {

		MsgSend( "waithere" );

		if( !MsgWaitAck( nw, jobstart_s ) )
			exit( 32 );

		MsgSend( "allsynced" );

		if( !MsgWaitAck( nw, jobstart_s ) )
			exit( 32 );

		MsgSend( "begin" );

		if( !MsgWaitAck( nw ) )
			exit( 32 );
	}
	else {
		if( !MsgWaitMatch( "waithere", jobstart_s ) )
			exit( 32 );

		MsgAck( iw );

		if( !MsgWaitMatch( "allsynced", jobstart_s ) )
			exit( 32 );

		MsgAck( iw );

		if( !MsgWaitMatch( "begin" ) )
			exit( 32 );

		MsgAck( iw );
	}

	StopTiming( stdout, "Sync workers", t0 );
}


void MsgClear()
{
// Remove

	system( "rm -rf MSG" );

// Verify

	const int timeout_s	= 10;	// 10 seconds

	clock_t	t0 = StartTiming();

	for(;;) {

		if( !DskExists( "MSG" ) )
			return;
		else if( DeltaSeconds( t0 ) > timeout_s )
			break;

		usleep( 10 );	// microsec
	}

	printf( "MsgClear: Timed out after %s sec.\n", timeout_s );
	exit( 32 );
}


static void MsgNew()
{
// Create

	DskCreateDir( "MSG", stdout );

// Verify

	const int timeout_s	= 10;	// 10 seconds

	clock_t	t0 = StartTiming();

	for(;;) {

		if( DskExists( "MSG" ) )
			return;
		else if( DeltaSeconds( t0 ) > timeout_s )
			break;

		usleep( 10 );	// microsec
	}

	printf( "MsgNew: Timed out after %s sec.\n", timeout_s );
	exit( 32 );
}


void MsgSend( const char *msg )
{
	MsgClear();
	MsgNew();

	FILE	*f = fopen( "MSG/msg", "w" );

	if( f ) {
		fprintf( f, msg );
		fclose( f );
	}
}


bool MsgWaitAck(
	int		nw,
	int		timeout_s )
{
	clock_t		t0 = StartTiming();
	vector<int>	ack( nw, 0 );
	int			sleep_us = (timeout_s == MSG_STD_TO_SEC ? 10 : 500000);

	for(;;) {

		int	nack = 1;

		for( int iw = 1; iw < nw; ++iw ) {

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
		else if( DeltaSeconds( t0 ) > timeout_s )
			break;

		usleep( sleep_us );	// microsec
	}

	printf( "MsgWaitAck: Timed out after %s sec.\n", timeout_s );
	return false;
}


bool MsgWaitRecv(
	char	msgbuf[128],
	int		timeout_s )
{
	clock_t	t0 = StartTiming();
	int		sleep_us = (timeout_s == MSG_STD_TO_SEC ? 10 : 500000);

	for(;;) {

		if( DskExists( "MSG/msg" ) ) {

			FILE	*f = fopen( "MSG/msg", "r" );

			if( f ) {
				fscanf( f, "%127s", msgbuf );
				fclose( f );
			}

			return true;
		}

		if( DeltaSeconds( t0 ) > timeout_s )
			break;

		usleep( sleep_us );	// microsec
	}

	printf( "MsgWaitRecv: Timed out after %s sec.\n", timeout_s );
	return false;
}


bool MsgVerify( const char msgbuf[128], const char *smatch )
{
	if( !strcmp( msgbuf, smatch ) )
		return true;

	printf( "MsgVerify: Got [%s] expected [%s].\n", msgbuf, smatch );
	return false;
}


bool MsgWaitMatch(
	const char	*smatch,
	int			timeout_s )
{
	clock_t	t0 = StartTiming();
	int		sleep_us = (timeout_s == MSG_STD_TO_SEC ? 10 : 500000);

	for(;;) {

		if( DskExists( "MSG/msg" ) ) {

			char	msg[128];
			FILE	*f = fopen( "MSG/msg", "r" );

			if( f ) {
				fscanf( f, "%127s", msg );
				fclose( f );
			}

			if( !strcmp( msg, smatch ) )
				return true;
		}

		if( DeltaSeconds( t0 ) > timeout_s )
			break;

		usleep( sleep_us );	// microsec
	}

	printf( "MsgWaitMatch: Timed out after %s sec.\n", timeout_s );
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


