

#pragma once


/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

const int MSG_STD_TO_SEC = 30;






/* --------------------------------------------------------------- */
/* Universal functions ------------------------------------------- */
/* --------------------------------------------------------------- */

void MsgSyncWorkers( int iw, int nw );

/* --------------------------------------------------------------- */
/* Master functions ---------------------------------------------- */
/* --------------------------------------------------------------- */

void MsgClear();

void MsgSend( const char *msg );

bool MsgWaitAck(
	int		nw,
	int		timeout_s = MSG_STD_TO_SEC );

/* --------------------------------------------------------------- */
/* Worker functions ---------------------------------------------- */
/* --------------------------------------------------------------- */

bool MsgWaitRecv(
	char	msgbuf[128],
	int		timeout_s = MSG_STD_TO_SEC );

bool MsgVerify( const char msgbuf[128], const char *smatch );

bool MsgWaitMatch(
	const char	*smatch,
	int			timeout_s = MSG_STD_TO_SEC );

void MsgAck( int iw );


