#ifndef __GAUSSIANS_H__
#define __GAUSSIANS_H__

#include <math.h>

inline double get_S(unsigned i, unsigned j, unsigned k, unsigned l, const double alpha, const double a0)
{
    double a,A,B,C,D;
    double S=0.;
    a=a0;
    A=alpha/(i*i*i);
    B=alpha/(j*j*j);
    C=alpha/(k*k*k);
    D=alpha/(l*l*l);
    S=pow(M_PI,(3./2.))/(pow((A+B+C+D),(3./2.)));
    S=S*exp(-4.*a*a*(A+C)*(B+D)/(A+B+C+D));
    return S;
}

inline double get_T(unsigned i, unsigned j, unsigned k, unsigned l, const double alpha, const double a0)
{
    double a,A,B,C,D;
    double aux1,aux2;
    double T=0.;
    a=a0;
    A=alpha/(i*i*i);
    B=alpha/(j*j*j);
    C=alpha/(k*k*k);
    D=alpha/(l*l*l);
    aux1=(3.*(A+B)*(C+D))/(pow((A+B+C+D),(5./2.)));
    aux2=(8.*a*a*(B*C-A*D)*(B*C-A*D))/(pow((A+B+C+D),(7./2.)));
    T=pow(M_PI,(3./2.))*(aux1-aux2);
    T=T*exp(-4.*a*a*(A+C)*(B+D)/(A+B+C+D));
    return T;
}

inline double get_V(unsigned i, unsigned j, unsigned k, unsigned l, unsigned part, const double alpha, const double a0)
{
    double a,A,B,C,D;
    double aux1,aux2;
    double V=0.;
    a=a0;
    A=alpha/(i*i*i);
    B=alpha/(j*j*j);
    C=alpha/(k*k*k);
    D=alpha/(l*l*l);
    if(part==1) aux1=B+D;
    if(part==2) aux1=A+C;
    aux2=2.*(a*aux1)/(sqrt(A+B+C+D));
    V=(-1.)*pow(M_PI,(3./2.))*(1./(2.*a*aux1*sqrt(A+B+C+D)));
    V=V*erf(aux2);
    V=V*exp(-4.*a*a*(A+C)*(B+D)/(A+B+C+D));
    return V;
}

#endif