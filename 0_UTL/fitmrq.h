

#pragma once


typedef void (*Fitfunc)( Doub, VecDoub_I &, Doub &, VecDoub_O & );


struct Fitmrq {
    static const Int	NDONE=4, ITMAX=1000;

    Int			ndat, ma, mfit;
    VecDoub_I	&x,&y,&sig;
    Doub		tol;
    Fitfunc		funcs;
    VecBool		ia;
    VecDoub		a;
    MatDoub		covar;
    MatDoub		alpha;
    Doub		chisq;

    Fitmrq(
        VecDoub_I	&xx,
        VecDoub_I	&yy,
        VecDoub_I	&ssig,
        VecDoub_I	&aa,
        Fitfunc		funks,
        Doub		TOL = 1.e-3 )
        :
            ndat(xx.size()), ma(aa.size()),
            x(xx), y(yy), sig(ssig),
            tol(TOL), funcs(funks),
            ia(ma), a(aa),
            alpha(ma,ma), covar(ma,ma)
        {
            for( Int i = 0; i < ma; ++i )
                ia[i] = true;
        };

    void hold( Int i, Doub val )	{ia[i]=false; a[i]=val;};
    void free( Int i )				{ia[i]=true;};

    void fit();

    void mrqcof( VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta );
    void covsrt( MatDoub_IO &covar );
};


