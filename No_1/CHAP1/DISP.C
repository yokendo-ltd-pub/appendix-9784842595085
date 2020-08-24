#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g.h"

#define NP       1000 /* number of nodes */
#define NE       2000 /* number of elements */
#define NSA      21   /* number of contour lines */
#define BUF_size 256

void Input( int argc, char *argv[] );
void Element( void );
void Element_3fem( void );
void Element_tri( void );
void Contour( void );
void Vector_el( void );
void Vector( void );
void arrow( double x, double y, double z,
            double vx, double vy, double vz, double scale );
void MoveR( double dx, double dy, double cs, double sn,
            double X, double Y );
void LineR( double dx, double dy , double cs, double sn,
            double X, double Y );
void Point_setR( double dx, double dy , double cs, double sn,
            double X, double Y );
void Ellipse( double xx, double yy, double a, double b, 
              double t1, double t2, double cs, double sn,
              double X, double Y );
void Dtoe( char *s );

int np,ne,nv = 0,nvec,nsa = NSA,Noe[NE][9],Ns[NE];
double Xc[NP],Yc[NP],Zc[NP],Sd[NP],Vd[NE][10],Sa[52],Vk;
double xv = 0., yv = 0., zv = 1.;
char buf[BUF_size],*sp,*fname,*fnames,*fnamev;

int flg = 1, ri = 0, id = 4; 
 /* flg = 0: 3fem element, 1: 3bem element, 2: other element, 
          3: Contour,      4: Vector,       5: Vector on element */
 /* ri = 0: real,  1: imaginary */

main( int argc, char *argv[] )
{
   int i,ixd,iyd;
   double xd,yd;
   char text[10];
   int lp = 0;

   printf("\n <*><*><*> Disp <*><*><*>\n\n");
   Ginit( argc, argv );
   Input( argc, argv );
   printf("\n---------------------------------\n");
   printf(  " Mouse Button: Left   -> enlarge\n");
   printf(  "               Center -> move\n");
   printf(  "               Right  -> reduce\n");
   printf(  " Key         : Q or q -> quit.\n");
   printf(  "---------------------------------\n\n");

   while ( lp == 0 ) {

      XNextEvent ( Ddisplay, &Devent );
      switch ( Devent.type ) {

      case Expose:
         if ( Devent.xexpose.count == 0 )
         {  
            if( XGetGeometry( Ddisplay, Dwindow, &root, &xr, &yr, 
                    &Wwdth, &Whigt, &bw, &dp ) != 0 ) {
               IXw = (double) Wwdth;
               IYw = (double) Whigt;
               Gclr();
               Set_Vwindow( Xmin, Ymin, Xmin+Xw, Ymin+Yw );
               if( flg == 0 )      Element_3fem();
               else if( flg == 1 ) Element_tri();
               else if( flg == 2 ) Element();
               else if( flg == 3 ) Contour();
               else if( flg == 4 ) Vector();
               else                Vector_el();
               }
	  }
         break;

       case MappingNotify:
          XRefreshKeyboardMapping ( &Devent );
          break;

       case ButtonPress:
          Get_Vco( Devent.xbutton.x, Devent.xbutton.y, &xd, &yd );
          if( Devent.xbutton.button == 1 ) {
             Xw *= 0.5;
             Yw *= 0.5; }
          else if( Devent.xbutton.button == 3 ) {
             Xw *= 2.;
             Yw *= 2.; }
          Xmin = xd - 0.5 * Xw;
          Ymin = yd - 0.5 * Yw;
          Gclr();
          Set_Vwindow( Xmin, Ymin, Xmin+Xw, Ymin+Yw );
          if( flg == 0 )      Element_3fem();
          else if( flg == 1 ) Element_tri();
          else if( flg == 2 ) Element();
          else if( flg == 3 ) Contour();
          else if( flg == 4 ) Vector();
          else                Vector_el();
          break;

       case KeyPress:
          i = XLookupString ( &Devent, text, 10, &Dkey, 0 );
          if ( i == 1 && ( text[0] == 'q' || text[0] == 'Q' ) ) lp = 1;
          break;
       }
    }
   Gfinish();
   exit( 0 );
 }

void Input( int argc, char *argv[] )
{
   FILE *fp;
   int i,j,ip = 0; /* when ip=1, input data => display for debuging */ 
   double Xxmin = 1.e6, Xxmax = -1.e6, Yymin = 1.e6, Yymax = -1.e6,
          Zzmin = 1.e6, Zzmax = -1.e6, dummy;
   double ds, Smin = 1.e6, Smax = -1.e6, vabs, Vmax = -1e6;
   static char dpara[5][5] = { "2fem", "2bem", "3fem", "3bem", " " }, *dst;
   char *part,*sx,*sy,*sz;
   int igf = 0, isf = 0, ivf = 0;;

   for( i=1; i<argc; i++ )
      if( strcmp("-vew",argv[i]) == 0 )  {
         sx = argv[ i + 1 ];
         sy = argv[ i + 2 ];
         sz = argv[ i + 3 ];
         if( sx == NULL || sy == NULL || sz == NULL ) break;
         j = 0;
         if( sscanf( sx, "%lf", &xv ) == 1 ) j++;
         if( sscanf( sy, "%lf", &yv ) == 1 ) j++;
         if( sscanf( sz, "%lf", &zv ) == 1 ) j++;
         if( j != 3 ) {
            xv = 0.;   yv = 0.;   zv = 1.;
            break; }
         printf(" xvew= (%7.3lf,%7.3lf,%7.3lf )\n", xv, yv, zv );
         break; }

   for( i=1; i<argc; i++ )
      if( strcmp("-d",argv[i]) == 0 )  {
         dst = argv[ i + 1 ];
         if( dst == NULL ) break;
         for( j=0; j<4; j++) 
            if( strcmp( dst, &dpara[j][0] ) == 0 ) { 
               id = j;
               goto Find_dp; }
       }
 Find_dp:
   printf(" data = %s\n",&dpara[id][0]);

   for( i=1; i<argc; i++ )
      if( strcmp("-gf",argv[i]) == 0 )  {
         fname = argv[ i + 1 ];
         igf = 1;
         printf(" g_file= '%s'\n", fname );
         break; }

   for( i=1; i<argc; i++ )
      if( strcmp("-sf",argv[i]) == 0 )  {
         fnames = argv[ i + 1 ];
         part = argv[ i + 2 ];
         isf = 1;
         if( part != NULL && *part == 'i' ) {
            ri = 1;
            printf(" s_file= '%s' :imag\n", fnames ); }
         else {
            printf(" s_file= '%s' :real\n", fnames ); }
         break; }

   for( i=1; i<argc; i++ )
      if( strcmp("-vf",argv[i]) == 0 )  {
         fnamev = argv[ i + 1 ];
         part = argv[ i + 2 ];
         ivf = 1;
         if( part != NULL && *part == 'i' ) {
            ri = 1;
            printf(" v_file= '%s' :imag\n", fnamev ); }
         else {
            printf(" v_file= '%s' :real\n", fnamev ); }
         break; }

   if( id == 2 ) flg = 0;
   else if( igf == 0 )
           if( ivf == 1 ) flg = 4;
           else {
              printf("\n input error!\n\n");
              exit(0); }
   else if( isf == 1 )   flg = 3;
   else if( ivf == 1 )   flg = 5;
   else if( id == 3 )    flg = 1;
   else                  flg = 2;


   printf(" %1d:  ",flg );

   if( flg == 4 ) goto Vector;

   for( i=1; i<=np; i++ ) Zc[i] = 0.;

   if( ( fp = fopen( fname, "r" ) ) == NULL ) {
      printf("Can't open '%s'!\n",fname);
      exit(0); } 

   sp = fgets( buf, BUF_size, fp );
   sscanf( sp,"%d %d %d", &np, &ne, &nv);
   if( nv == 0 ) nv = 3;
   printf(" np=%d,   ne=%d,   nv=%d\n",np,ne,nv);

   for( i=1; i<=np; i++ ) {
      if( ( sp = fgets( buf, BUF_size, fp ) ) == NULL ) {
         printf("\nEOF or Err in reading Z[][]");
         exit( 0 ); }
      Dtoe( sp );
      sscanf( sp, "%le %le %le", &Xc[i],&Yc[i],&Zc[i]);
      if( ip == 1 ) printf("Z(%6d) %10.4lf %10.4lf %10.4lf\n",
                             i,Xc[i],Yc[i],Zc[i]);
    }

   for( i=1; i<=ne; i++ ) {
      if( ( sp = fgets( buf, BUF_size, fp ) ) == NULL ) {
         printf("\nEOF or Err in reading Noe[][]");
         exit( 0 ); }
      for( j=1; j<=nv; j++ ) {
         while( *sp == ' ' ) sp++;
         sscanf( sp, "%d", &Noe[i][j]);
         while( *sp != ' ' ) sp++;
       }
      if( ip == 1 ) { 
         printf("Noe(%6d) ",i);
         for( j=1; j<=nv; j++ ) printf("%6d",Noe[i][j]);
         printf("\n"); }
    }

End_input:
   fclose( fp );
   if( ip == 1 ) printf("Input end!\n\n");

   for( i=1; i<=np; i++ ) {
      if( Xc[i] > Xxmax ) Xxmax = Xc[i];
      if( Xc[i] < Xxmin ) Xxmin = Xc[i];
      if( Yc[i] > Yymax ) Yymax = Yc[i];
      if( Yc[i] < Yymin ) Yymin = Yc[i];
      if( Zc[i] > Zzmax ) Zzmax = Zc[i];
      if( Zc[i] < Zzmin ) Zzmin = Zc[i];
    }

   Set_Tmat( xv, yv, zv );
   Set_Vwindow3D( Xxmin, Yymin, Zzmin, Xxmax, Yymax, Zzmax );

   if( flg == 5 ) goto Vector;
   if( flg <= 2 ) return; 

   if( ( fp = fopen( fnames, "r" ) ) == NULL ) {
      printf("Can't open '%s'!\n",fnames);
      exit(0); } 

   if( ri == 0 ) {
      for( i=1; i<=np; i++ ) {
         if( ( sp = fgets( buf, BUF_size, fp ) ) == NULL ) {
            printf("\nEOF or Err in reading Sd[]");
            exit( 0 ); }
         Dtoe( sp );
         sscanf( sp, "%le %le", &Sd[i], &dummy );
         if( ip == 1 ) printf("Sd[%6d] %15.8le\n",i,Sd[i]);
       }
   }
   else {
      for( i=1; i<=np; i++ ) {
         if( ( sp = fgets( buf, BUF_size, fp ) ) == NULL ) {
            printf("\nEOF or Err in reading Sd[]");
            exit( 0 ); }
         Dtoe( sp );
         sscanf( sp, "%le %le", &dummy, &Sd[i] );
         if( ip == 1 ) printf("Sd[%6d] %15.8le\n",i,Sd[i]);
       } 
   }
   fclose( fp );
   for( i=1; i<=np; i++ ) {
      if( Sd[i] > Smax ) Smax = Sd[i];
      if( Sd[i] < Smin ) Smin = Sd[i];
    }
   printf(" Smin=%le,  Smax=%le\n", Smin, Smax);
   ds = ( Smax - Smin ) / 20.;
   j = (int)(Smin / ds );
   for( i=1; i<=nsa; i++ ) {
      Sa[i] = ds * ( i + j - 2 );
    }
   return;

 Vector:
   if( ( fp = fopen( fnamev, "r" ) ) == NULL ) {
      printf("Can't open '%s'!\n",fnamev);
      exit(0); } 

   nvec = 1;
   j = 1;
   for(;;) {
 Next:
      if( ( sp = fgets( buf, BUF_size, fp ) ) == NULL )
         break;
      Dtoe( sp );
      for(;;) {
         while( *sp == ' ' ) sp++;
         if( sscanf( sp, "%le", &Vd[nvec][j]) != 1 ) goto Next;
         j++;
         if( j == 10 ) {
            if( ip == 1 ) { 
               printf("Vd(%6d) ",nvec);
               for( j=1; j<=9; j++ ) {
                  printf("%10.2le",Vd[nvec][j]);
                  if( j == 3 ) printf("\n"); }
               printf("\n"); }
            j = 1;
            nvec++;
            goto Next; }
         while( *sp != ' ' && *sp != '\0' ) sp++;
      }
    }
    fclose( fp );
    nvec--;
    printf(" nvec=%d\n",nvec);
    if( ri == 1 ) 
       for( i=1; i<=nvec; i++ ) {
             for( j=1; j<=3; j++ )
                Vd[i][j+3] = Vd[i][j+6]; }

    if( flg == 5 ) goto Vscale;

    for( i=1; i<=nvec; i++ ) {
       if( Vd[i][1] > Xxmax ) Xxmax = Vd[i][1];
       if( Vd[i][1] < Xxmin ) Xxmin = Vd[i][1];
       if( Vd[i][2] > Yymax ) Yymax = Vd[i][2];
       if( Vd[i][2] < Yymin ) Yymin = Vd[i][2];
       if( Vd[i][3] > Zzmax ) Zzmax = Vd[i][3];
       if( Vd[i][3] < Zzmin ) Zzmin = Vd[i][3];
     }
/*  printf("Xmin=%10.3lf, Xmax=%10.3lf\n",Xxmin,Xxmax);
    printf("Ymin=%10.3lf, Ymax=%10.3lf\n",Yymin,Yymax);
    printf("Zmin=%10.3lf, Zmax=%10.3lf\n",Zzmin,Zzmax); */

    Set_Tmat( xv, yv, zv );
    Set_Vwindow3D( Xxmin, Yymin, Zzmin, Xxmax, Yymax, Zzmax );

 Vscale:
    for( i=1; i<=nvec; i++ ) {
       vabs = Vd[i][4] * Vd[i][4] + Vd[i][5] * Vd[i][5] + Vd[i][6] * Vd[i][6];
       if( Vmax < vabs ) Vmax = vabs; }
    Vmax = sqrt( Vmax );
    Vk = Xw / 20. / Vmax;
    printf(" Vmax=%10.3le,  Vk=%10.3le\n",Vmax,Vk);

}

void Element( void )
{
   int i,j,nj;

   for( i=1; i<=ne; i++ ) {
      nj = Noe[i][3];
      Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
      for( j=1; j<=nv; j++ ) {
         nj = Noe[i][j];
         Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
       }
    }
 }

void Element_tri( void ) /* for Triangular Elements */
{
   int i,ni,nj,nm;
   double ax,ay,az,bx,by,bz,nx,ny,nz,xg,yg,zg;

   for( i=1; i<=ne; i++ ) {

      ni = Noe[i][1];
      nj = Noe[i][2];
      nm = Noe[i][3];

      ax = Xc[nj] - Xc[ni];
      ay = Yc[nj] - Yc[ni];
      az = Zc[nj] - Zc[ni];
      bx = Xc[nm] - Xc[ni];
      by = Yc[nm] - Yc[ni];
      bz = Zc[nm] - Zc[ni];
      nx = ay * bz - az * by;
      ny = az * bx - ax * bz;
      nz = ax * by - ay * bx;
      zg = nx * xv + ny * yv + nz * zv;
      if( zg < 0. ) continue;

      Move_to3D( Xc[ni], Yc[ni], Zc[ni] );
      Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
      Line_to3D( Xc[nm], Yc[nm], Zc[nm] );
      Line_to3D( Xc[ni], Yc[ni], Zc[ni] );

      xg = ( Xc[ni] + Xc[nj] + Xc[nm] ) / 3.;
      yg = ( Yc[ni] + Yc[nj] + Yc[nm] ) / 3.;
      zg = ( Zc[ni] + Zc[nj] + Zc[nm] ) / 3.;
      Point_set3D( xg, yg, zg );
    }
 }

void Contour( void ) /* for Triangular Elements */
{

   int i,j,k,m,ni,nj,nm;
   double ax,ay,az,bx,by,bz,nx,ny,nz,x,y,z,SS,w;

   for( i=1; i<=ne; i++ ) {

      ni = Noe[i][1];
      nj = Noe[i][2];
      nm = Noe[i][3];

      ax = Xc[nj] - Xc[ni];
      ay = Yc[nj] - Yc[ni];
      az = Zc[nj] - Zc[ni];
      bx = Xc[nm] - Xc[ni];
      by = Yc[nm] - Yc[ni];
      bz = Zc[nm] - Zc[ni];
      nx = ay * bz - az * by;
      ny = az * bx - ax * bz;
      nz = ax * by - ay * bx;
      z = nx * xv + ny * yv + nz * zv;
      if( z < 0. ) continue;

      for( m=1; m<=nsa; m++ ) {
         SS = Sa[m];

         if( Sd[ni] == SS && Sd[nj] == SS ) {
            Move_to3D( Xc[ni], Yc[ni], Zc[ni] );
            Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
            continue; }
         if( Sd[nj] == SS && Sd[nm] == SS ) {
            Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
            Line_to3D( Xc[nm], Yc[nm], Zc[nm] );
            continue; }
         if( Sd[ni] == SS && Sd[nm] == SS ) {
            Move_to3D( Xc[ni], Yc[ni], Zc[ni] );
            Line_to3D( Xc[nm], Yc[nm], Zc[nm] );
            continue; }

         k = 0;
         if( ( Sd[ni] - SS ) * ( Sd[nj] - SS ) > 0. ) goto N1;
            w = ( SS - Sd[ni] ) / ( Sd[nj] - Sd[ni] );
            x = Xc[ni] + w * ( Xc[nj] - Xc[ni] );
            y = Yc[ni] + w * ( Yc[nj] - Yc[ni] );
            z = Zc[ni] + w * ( Zc[nj] - Zc[ni] );
         if( k == 0 ) {
            Move_to3D( x, y, z );
            k = 1; }
         else {
            Line_to3D( x, y, z );
            continue; }
 N1:
         if( ( Sd[nj] - SS ) * ( Sd[nm] - SS ) > 0. ) goto N2;
            w = ( SS - Sd[nj] ) / ( Sd[nm] - Sd[nj] );
            x = Xc[nj] + w * ( Xc[nm] - Xc[nj] );
            y = Yc[nj] + w * ( Yc[nm] - Yc[nj] );
            z = Zc[nj] + w * ( Zc[nm] - Zc[nj] );
         if( k == 0 ) {
            Move_to3D( x, y, z );
            k = 1; }
         else {
            Line_to3D( x, y, z );
            continue; }
 N2:
         if( ( Sd[nm] - SS ) * ( Sd[ni] - SS ) > 0. ) continue;
            w = ( SS - Sd[nm] ) / ( Sd[ni] - Sd[nm] );
            x = Xc[nm] + w * ( Xc[ni] - Xc[nm] );
            y = Yc[nm] + w * ( Yc[ni] - Yc[nm] );
            z = Zc[nm] + w * ( Zc[ni] - Zc[nm] );
         if( k == 0 ) {
            Move_to3D( x, y, z );
            k = 1; }
         else {
            Line_to3D( x, y, z );
            continue; }
       }

    }
}

void Vector_el( void )
{
   int i,ni,nj,nm;
   double ax,ay,az,bx,by,bz,nx,ny,nz,xg,yg,zg;

   for( i=1; i<=ne; i++ ) {

      ni = Noe[i][1];
      nj = Noe[i][2];
      nm = Noe[i][3];

      ax = Xc[nj] - Xc[ni];
      ay = Yc[nj] - Yc[ni];
      az = Zc[nj] - Zc[ni];
      bx = Xc[nm] - Xc[ni];
      by = Yc[nm] - Yc[ni];
      bz = Zc[nm] - Zc[ni];
      nx = ay * bz - az * by;
      ny = az * bx - ax * bz;
      nz = ax * by - ay * bx;
      zg = nx * xv + ny * yv + nz * zv;
      if( zg < 0. ) continue;

      arrow( Vd[i][1], Vd[i][2], Vd[i][3], Vd[i][4], Vd[i][5], Vd[i][6], Vk );

    }
}

void Vector( void )
{
    int i,j,k;

    for( i=1; i<=nvec; i++ ) {
       arrow( Vd[i][1], Vd[i][2], Vd[i][3], Vd[i][4], Vd[i][5], Vd[i][6], Vk );
     }
}

void arrow( double x, double y, double z,
            double vx, double vy, double vz, double scale )
{
    int i,j,k;
    double Vx,Vy,Vz,VV,Vs,ct,st,C,CC,BC,TC,WC,SCC,WS;
    double Xx,Yy;

    trans3( vx, vy, vz, &Vx, &Vy, &Vz );
    trans( x, y, z, &Xx, &Yy );

    VV = sqrt( Vx * Vx + Vy * Vy + Vz * Vz );
    SCC = Vz / VV;
    CC = sqrt( 1. - SCC * SCC );
    C = CC * SCC / fabs( SCC );
    Vs = sqrt( Vx * Vx + Vy * Vy );
    ct = Vx / Vs;
    st = Vy / Vs;
    VV *= scale;
    TC = VV * CC;
    BC = 0.6 * TC;
    WC = 0.1 * VV;
    WS = WC * SCC;

    if( C < 0. ) {
       if( CC > 0.35 ) {
          MoveR( BC, -WC, ct, st, Xx, Yy );
          LineR( TC, 0., ct, st, Xx, Yy );
          LineR( BC, WC, ct, st, Xx, Yy ); }
       Ellipse( BC, 0., WS, WC, 0., 2. * PI, ct, st, Xx, Yy );
       Move_to3D( x, y, z );
       LineR( BC, 0., ct, st, Xx, Yy ); }
    else {
       if( BC > WS ) {
          Move_to3D( x, y, z );
          LineR( BC - WS, 0., ct, st, Xx, Yy ); }
       Ellipse( BC, 0., WS, WC, 0.5 * PI, 1.5 * PI, ct, st, Xx, Yy );     
       if( CC > 0.3 ) {
          MoveR( BC, -WC, ct, st, Xx, Yy );
          LineR( TC, 0., ct, st, Xx, Yy );
          LineR( BC, WC, ct, st, Xx, Yy ); }
       else {
          Point_setR( TC, 0., ct, st, Xx, Yy );
          Ellipse( BC, 0., WS, WC, -0.5 * PI, 0.5 * PI, ct, st, Xx, Yy );
        }
     }
  }

void MoveR( double dx, double dy, double cs, double sn,
            double X, double Y )
{
   double Xx, Yy;

   Xx = dx * cs - dy * sn + X;
   Yy = dx * sn + dy * cs + Y;
   Move_to( Xx, Yy );
}

void LineR( double dx, double dy , double cs, double sn,
            double X, double Y )
{
   double Xx, Yy;

   Xx = dx * cs - dy * sn + X;
   Yy = dx * sn + dy * cs + Y;
   Line_to( Xx, Yy );
}

void Point_setR( double dx, double dy , double cs, double sn,
            double X, double Y )
{
   double Xx, Yy;

   Xx = dx * cs - dy * sn + X;
   Yy = dx * sn + dy * cs + Y;
   Point_set( Xx, Yy );
}

void Ellipse( double x0, double y0, double a, double b, 
              double t1, double t2, double cs, double sn,
              double X, double Y )
{
   static double db = 3.1415926 / 6.;
   int i,n;
   double dt, dx, dy, tt;

   n = (int)( ( t2 - t1 ) / db ) + 1;
   dt = ( t2 - t1 ) / (double)n;
   dx = a * cos( t1 );
   dy = b * sin( t1 );
   MoveR( x0 + dx, y0 + dy, cs, sn, X, Y );
   for( i=1; i<=n; i++ ) {
      tt = t1 + dt * (double)i;
      dx = a * cos( tt );
      dy = b * sin( tt );
      LineR( x0 + dx, y0 + dy, cs, sn, X, Y );
    }
}

void Element_3fem( void )
{
   int i,j,nj;

   if( nv == 8 ) goto Brick;

   for( i=1; i<=ne; i++ ) {
      nj = Noe[i][3];
      Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
      for( j=1; j<=nv; j++ ) {
         nj = Noe[i][j];
         Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
       }
      nj = Noe[i][1];
      Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
      nj = Noe[i][4];
      Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
      nj = Noe[i][2];
      Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
    }
    return;

 Brick:
   for( i=1; i<=ne; i++ ) {
      nj = Noe[i][4];
      Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
      for( j=1; j<=4; j++ ) {
         nj = Noe[i][j];
         Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
       }
      nj = Noe[i][8];
      Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
      for( j=5; j<=8; j++ ) {
         nj = Noe[i][j];
         Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
       }
      for( j=1; j<=4; j++ ) {
         nj = Noe[i][j];
         Move_to3D( Xc[nj], Yc[nj], Zc[nj] );
         nj = Noe[i][j+4];
         Line_to3D( Xc[nj], Yc[nj], Zc[nj] );
       }
    }
 }

void Dtoe( char *s )
{
   while( *s != '\0' ) {
      if( *s == 'D' || *s == 'd' ) *s = 'e';
      s++; }
 }

