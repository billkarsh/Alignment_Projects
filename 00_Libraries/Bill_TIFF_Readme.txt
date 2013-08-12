Set Up
-------
Download latest stable version (here, tiff-3.9.5.tar.gz).

Right click in Windows Explorer and select WinZip to unpack.
(In linux NX, just double-click gz file and drag folder to extract).

Read help at /html/index.html
Nav to
'Building the software distribution'
...'Building on a UNIX system.'

Useful to clarify syntax and options:
> ./configure --help

Then,
> ./configure CFLAGS=-O3 --prefix=/groups/apig/share/ClusterSupport/Libraries/TIFF --with-zlib-include-dir=/groups/apig/share/ClusterSupport/Libraries/ZLIB/include --with-zlib-lib-dir=/groups/apig/share/ClusterSupport/Libraries/ZLIB/lib

> make
> make install


