#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define dx (M_PI/10.0)
#define Y 0.4
#define x_end M_PI
#define t_end 10.0
#define D 1.0
#define n_initial 1.0
#define dt (Y*dx*dx/D)

#define x_size (int)(ceil(x_end/dx))
#define t_size (int)(ceil(t_end/dt))

void initU(double u[t_size+1][x_size+1]);
void boundaries(double u[t_size+1][x_size+1]);
void fill_1(double u_1[t_size+1][x_size+1]);
double exact(double t,double x);
double error(double exact1,double approx);
void printAll_1(double u_1[t_size+1][x_size+1]);

int main()
{
    double u_1[t_size+1][x_size+1];

    initU(u_1);

    boundaries(u_1);

    fill_1(u_1);

    printAll_1(u_1);

    printf("\n\n");

    return 0;
}

void initU(double u[t_size+1][x_size+1])
{
    for(int i=0;i<=x_size;i++)
    {
        u[0][i]=sin(i*dx*n_initial);
    }
}

void boundaries(double u[t_size+1][x_size+1])
{
    for(int n=0;n<=t_size;n++)
    {
        u[n][0]=0;
        u[n][x_size]=0;
    }

}

void fill_1(double u_1[t_size+1][x_size+1])
{
    int i,w,j,row;

    double alpha[x_size],g[x_size];

    double a=1.0+2.0*Y;
    double b=-Y;
    double c=-Y;

    alpha[1]=a;
    for(row=1;row<=t_size;row++)
    {
        g[1]=u_1[row-1][1];

        for(j=2;j<x_size;j++)
        {
            alpha[j]=a-c*b/alpha[j-1];
            g[j]=u_1[row-1][j]-b/alpha[j-1]*g[j-1];
        }

        u_1[row][x_size-1]=g[x_size-1]/alpha[x_size-1];

        for(w=x_size-2;w>0;w--)
        {
            u_1[row][w]=(g[w]-c*u_1[row][w+1])/alpha[w];
        }
    }
}

double exact(double t,double x)
{
    return exp(-n_initial*n_initial*t)*sin(n_initial*x);
}

double error(double exact1,double approx)
{
    return (fabs(approx-exact1)/fabs(approx))*100;
}

void printAll_1(double u_1[t_size+1][x_size+1])
{
    double ex=0,err=0;

    for(int n=0;n<=t_size;n++)
    {
        for(int j=0;j<=x_size;j++)
        {
            ex=exact(n*dt,j*dx);
            err=error(ex,u_1[n][j]);
            printf("u(%d)(%d): %3.20lf\t",n,j,u_1[n][j]);
            printf("Exact(%d): %3.20lf\t",j,ex);
            printf("Error(%d): %3.20lf\n",j,err);
        }
        printf("\n\n");
    }
}

