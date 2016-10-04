

#pragma once


#include	"GenDefs.h"

#include	<stdio.h>

#include	<vector>
using namespace std;


void FreeMRC( vector<uint16*> &vras );

int ReadRawMRCFile(
    vector<uint16*>	&vras,
    const char*		name,
    uint32			&w,
    uint32			&h,
    FILE*			flog = stdout,
    bool			transpose = false );

uint8* ReadAnMRCFile(
    const char*		name,
    uint32			&w,
    uint32			&h,
    FILE*			flog = stdout,
    bool			transpose = false,
    int				Ngauss = 2 );

int ReadMultiImageMRCFile(
    vector<uint8*>	&vras,
    const char*		name,
    uint32			&w,
    uint32			&h,
    FILE*			flog = stdout,
    bool			transpose = false,
    int				Ngauss = 2 );


