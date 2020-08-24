#include <stdio.h>
#define N 1001

void matrix( int No_of_unknown, double c[][N] ); 
void gauss( int No_of_unknown, double c[][N] );
         
main( void )
{
        static double c[N-1][N];
        char a[10]; 
        int n,n2;
        printf("\n### Gauss Elimination Method (C) ###\n");
        do {
                printf("\nN ? N<=%d \n",N-1);
                scanf("%d",&n);
        } while ( n < 1 || n > N ) ;
        n2 = n/2;
        matrix( n, c );
        printf("Start ( N= %d )  - Push any key & Cr -\n",n);
        scanf("%s",a);
        gauss( n, c );
        printf("\n%2d %e\n\n",n2,c[n2][n]);
}

void matrix( int n, double c[][N] )
{
        int i,j;
        for (i=1; i<n-1 ; i++) {
                for (j=0; j<=n; j++) c[i][j] = 1.001E-10;
                c[i][i-1] = -1.;  c[i][i] = 2.;  c[i][i+1] = -1.;
        }
        for (j=0; j<=n; j++) {
                c[0][j] = 1.001E-10; c[n-1][j] = 1.001E-10;
        }
        c[0][0] = 1.;  c[n-1][n-1] = 1.; c[n-1][n] = (double) n-1;
}

void gauss( int n, double c[][N] )
{
        int nm,kp,k,ii,nn,i,j;
        double ckk,cik;
        nm=n-2; nn=n-1;
        for (k=0; k<=nm; k++) {
                kp = k+1;  ckk = 1. / c[k][k];
                for (j=kp; j<=n; j++) c[k][j] *= ckk;
                for (i=kp; i<=nn; i++) {
                        cik = -c[i][k];
                        for (j=kp; j<=n; j++)
                                c[i][j] += cik * c[k][j];
                }
        }
        c[nn][n] /= c[nn][nn];
        for (ii=nm; ii>=0; ii--) {
                for (j=ii+1; j<=nn; j++ ) 
                   c[ii][n] -= c[ii][j] * c[j][n] ;
        }
}                

