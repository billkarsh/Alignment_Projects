

#pragma once


void svdcmp(
    double	**a,
    int		m,
    int		n,
    double	*w,
    double	**v );

void svbksb(
    double	**u,
    double	w[],
    double	**v,
    int		m,
    int		n,
    double	b[],
    double	x[] );


