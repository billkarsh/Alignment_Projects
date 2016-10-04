

#include	"Disk.h"

#include	<errno.h>
#include	<fcntl.h>
#include	<stdlib.h>
#include	<string.h>
#include	<sys/stat.h>
#include	<unistd.h>






/* --------------------------------------------------------------- */
/* Semaphores ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// To use this function, decide upon a name string to uniquely
// identify the resource/mutex and be sure that all clients use
// the identical name. For example, "anglefile_layer2".
//
// On calling Get( name ), the return value is either false, which
// means that we had an error creating/getting a mutex, or true,
// which means that the caller has a good mutex and has exclusive
// access to the resource.
//
// In either case, the caller must call Release promptly both to
// release the mutex and to clean up after Get().
//
bool CMutex::Get( const char *name, FILE *flog )
{
// Prepend forward slash

    m_name[0] = '/';

    strcpy( m_name + 1, name + (name[0] == '/') );

// Create - if not created, else open existing
//
// Note that the initial_value param (last) sets the count
// of allowed concurrent accessors. Each call to sem_wait
// will look at the current count. If the count is greater
// than zero, wait will decrement the count and return. If
// the current count is zero, then wait will block until an
// external agent increases the count (sem_post). Hence, in
// the typical exclusive access case, set this value to one.

    mutex = sem_open( m_name, O_CREAT,
                S_IRUSR | S_IWUSR |
                S_IRGRP | S_IWGRP |
                S_IROTH | S_IWOTH, 1 );

    if( mutex == SEM_FAILED ) {
        fprintf( flog,
        "Mutex failure: sem_open errno %d [%s].\n", errno, m_name );
        return false;
    }

// Wait for lock

    int	ok = !sem_wait( mutex );

    if( !ok ) {
        fprintf( flog,
        "Mutex failure: sem_wait errno %d [%s].\n", errno, m_name );
    }

    return ok;
}


void CMutex::Release()
{
    if( mutex != SEM_FAILED ) {

        sem_post( mutex );	// release lock
        sem_close( mutex );	// release resources
    }

// Dispose - deferred by kernel until all users call sem_close

    sem_unlink( m_name );
}

/* --------------------------------------------------------------- */
/* DskExists ----------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if named file or dir exists.
//
bool DskExists( const char *path )
{
    struct stat	info;

    return !stat( path, &info );
}

/* --------------------------------------------------------------- */
/* DskBytes ------------------------------------------------------ */
/* --------------------------------------------------------------- */

// Return object size.
//
double DskBytes( const char *path )
{
    struct stat	info;

    if( stat( path, &info ) )
        info.st_size = 0;

    return info.st_size;
}

/* --------------------------------------------------------------- */
/* DskCreateDir -------------------------------------------------- */
/* --------------------------------------------------------------- */

void DskCreateDir( const char *path, FILE* flog )
{
    if( !DskExists( path ) ) {

        if( mkdir( path, 0777 ) == -1 && errno != EEXIST ) {

            fprintf( flog,
            "Error %d creating dir [%s].\n", errno, path );
            exit( 42 );
        }
    }
}

/* --------------------------------------------------------------- */
/* DskAbsPath ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Simple brute force conversion to absolute path spec.
//
// Any input path starting with slash is already absolute,
// otherwise we will append the input to the current working
// directory. The following input shorthands are accepted:
//
// - A single instance of leading ./ is trimmed off.
// - One or more ../ specify parent directories.
//
// bufsize must be large enough to handle expected result,
// 2048 is a safe guess.
//
// Returns true if no errors, false if illegal input.
//
bool DskAbsPath( char *out, int bufsize, const char *in, FILE* flog )
{
// absolute already if starts with slash

    if( in[0] == '/' ) {
        strcpy( out, in );
        return true;
    }

// otherwise need current working dir

    getcwd( out, bufsize );

    const char	*in0 = in;
    char		*end = out + strlen( out );

next_parent:
    if( in[0] == '.' ) {

        // skip leading './'
        if( in[1] == '/' ) {
            in += 2;
            goto append;
        }
        else if( in[1] == '.' ) {

            if( in[2] == '/' ) {
                in += 3;
                end = strrchr( out, '/' );
                *end = 0;
                goto next_parent;
            }
            else
                goto illegal;
        }
        else if( !in[1] ) {
            ++in;
            goto append;
        }
        else {
illegal:
            fprintf( flog, "Illegal path '%s'.\n", in0 );
            return false;
        }
    }
    else {
append:
        if( in[0] ) {
            *end++ = '/';
            strcpy( end, in );
        }
    }

    return true;
}


