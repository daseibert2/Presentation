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
void fill(double u[t_size+1][x_size+1]);
int F(double t,double x);
double exact(double t,double x);
double error(double exact1,double approx);
void printAll(double u[t_size+1][x_size+1]);

int main()
{
    double u[t_size+1][x_size+1];

    initU(u);

    boundaries(u);

    fill(u);

    printAll(u);

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

void fill(double u[t_size+1][x_size+1])
{
    for(int m=1;m<=t_size;m++)
    {
        for(int j=1;j<x_size;j++)
        {
            u[m][j]=u[m-1][j]+Y*(u[m-1][j-1]-2*u[m-1][j]+u[m-1][j+1])+dt*F(j*dt,j*dx);
        }
    }
}

int F(double t,double x)
{
    return 0;
}

double exact(double t,double x)
{
    return exp(-n_initial*n_initial*t)*sin(n_initial*x);
}

double error(double exact1,double approx)
{
    return (fabs(approx-exact1)/fabs(approx))*100;
}

void printAll(double u[t_size+1][x_size+1])
{
    double ex=0,err=0;

    for(int n=0;n<=t_size;n++)
    {
        for(int j=0;j<=x_size;j++)
        {
            ex=exact(n*dt,j*dx);
            err=error(ex,u[n][j]);
            printf("u(%d)(%d): %3.20lf\t",n,j,u[n][j]);
            printf("Exact(%d): %3.20lf\t",j,ex);
            printf("Error(%d): %3.20lf\n",j,err);
        }
        printf("\n\n");
    }
}
