#pragma once

#ifdef __cplusplus
extern "C" {
#endif


void vec_copy(double *y, double const * x, int n)
{
    for(int i=0;i<n;i++) y[i]=x[i];
}

void vec_fill(double *y, double x, int n)
{
    for(int i=0;i<n;i++) y[i]=x;
}

/*
double vec_norm(double const *v, int n)
{
    double norm=0;
    for (int k=0;k<n;k++) norm+=v[k]*v[k];
    return sqrt(norm);
}

void vec_show(const double *x, int n, const char *name)
{
    printf("%s =\n[\n",name);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<3;j++) printf(" %+.4e",x[3*i+j]);
        printf("\n");
    }
    printf("]\n\n");
}

void headtail(double *x, int n, char const *name)
{

    int const  depth=9;
    printf("%s =\n[\n",name);
    for(int i=0;i<depth;i++)
    {
        printf(" %+.4e\n",x[i]);
    }
    printf("...\n");
    for(int i=n-depth;i<n;i++)
    {
        printf(" %+.4e\n",x[i]);
    }
    printf("]\n");

}
*/
#ifdef __cplusplus
}
#endif
