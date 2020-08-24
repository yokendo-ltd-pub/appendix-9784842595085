#include <X11/Xlib.h>
#include <X11/Xutil.h>

void Ginit( int argc, char **argv );
void Gfinish( void );
void Gclr( void );
void Set_Vwindow( double Xmin, double Ymin, double Xmax, double Ymax );
void Point_set( double x, double y );
void Line_to( double x, double y );
void Move_to( double x, double y );
void Set_Tmat( double x, double y, double z );
void Set_Vwindow3D( double Xmin, double Ymin, double Zmin,
                   double Xmax, double Ymax, double Zmax );
void Line_to3D( double x, double y, double z );
void Move_to3D( double x, double y, double z );
void Point_set3D( double x, double y, double z );
void trans( double Xi, double Yi, double Zi, double *Xo, double *Yo );
void trans3( double Xi, double Yi, double Zi, 
             double *Xo, double *Yo, double *Zo );
void Get_Vco( int ixd, int iyd, double *xd, double *yd );
double T11,T12,T21,T22,T23,T31,T32,T33;
double IXw, IYw, Kx, Ky, Xmin, Ymin, Xw, Yw;
int Xcc=0, Ycc=0;
double PI;

Display *Ddisplay;
Window Dwindow;
GC Dgc;
XEvent Devent;
KeySym Dkey;
XSizeHints Dhint;
int Dscreen;
unsigned long Dforeground, Dbackground;
int xr, yr;
Window root;
unsigned int Wwdth, Whigt, bw, dp;
char *display_name = NULL;
char Dname[] = { "Data" };

void Ginit( int argc, char **argv )
{
   int i;

   for( i=1; i<argc; i++ )
      if( strcmp("-display",argv[i]) == 0 )  {
         display_name = argv[ i + 1 ]; 
         printf(" Display = %s\n\n",display_name);
         break; }

   Ddisplay = XOpenDisplay( display_name );
   Dscreen = DefaultScreen ( Ddisplay );

   Dbackground = WhitePixel ( Ddisplay,Dscreen );
   Dforeground = BlackPixel ( Ddisplay,Dscreen );

   Dhint.x = 100;   Dhint.y = 100;
   Dhint.width = 500; Dhint.height = 450;
   Dhint.flags = PPosition | PSize;

   Dwindow = XCreateSimpleWindow ( Ddisplay ,
                 DefaultRootWindow ( Ddisplay),
                 Dhint.x, Dhint.y, Dhint.width, Dhint.height,
                 5, Dforeground, Dbackground);
   XSetStandardProperties ( Ddisplay, Dwindow, Dname, Dname,
                 None, argv, argc, &Dhint );

   Dgc = XCreateGC ( Ddisplay, Dwindow, 0, 0 );
   XSetBackground ( Ddisplay, Dgc, Dbackground );
   XSetForeground ( Ddisplay, Dgc, Dforeground );

   XSelectInput ( Ddisplay, Dwindow,
        ButtonPressMask | KeyPressMask | ExposureMask );

   XMapRaised ( Ddisplay, Dwindow );

   XGetGeometry( Ddisplay, Dwindow, &root, &xr, &yr, 
                 &Wwdth, &Whigt, &bw, &dp );
   IXw = (double) Wwdth;
   IYw = (double) Whigt;

   PI = atan( 1. ) * 4.;
}

void Gclr( void )
{
   XClearWindow( Ddisplay, Dwindow );
}

void Gfinish( void )
{
   XFreeGC ( Ddisplay, Dgc );
   XDestroyWindow ( Ddisplay, Dwindow );
   XCloseDisplay ( Ddisplay );
}

void Set_Vwindow( double X1, double Y1, double X2, double Y2 )
{
   double Z;
   Xmin = X1;
   Ymin = Y1;
   Xw = X2 - X1;
   Yw = Y2 - Y1;
   Kx = IXw / Xw;
   Ky = IYw / Yw;
   if( Kx < Ky ) {
      Z = Yw * Ky / Kx;
      Ymin -= ( Z - Yw ) * 0.5;
      Yw = Z;
      Ky = Kx; }
   else {
      Z = Xw * Kx / Ky;
      Xmin -= ( Z - Xw ) * 0.5;
      Xw = Z;
      Kx = Ky; }
}

void Get_Vco( int ixd, int iyd, double *xd, double *yd )
{
   *xd = (double)ixd / Kx + Xmin;
   *yd = ( IYw - (double)iyd ) / Ky + Ymin;
}

void Point_set( double x, double y )
{
   int Xnew, Ynew;
   Xnew = (int) ( Kx * ( x - Xmin ) );
   Ynew = (int) ( IYw - Ky * ( y - Ymin ) );
   XDrawPoint( Ddisplay, Dwindow, Dgc, Xnew, Ynew );
}

void Line_to( double x, double y )
{
   int Xnew, Ynew;
   Xnew = (int) ( Kx * ( x - Xmin ) );
   Ynew = (int) ( IYw - Ky * ( y - Ymin ) );
   XDrawLine( Ddisplay, Dwindow, Dgc, Xcc, Ycc, Xnew, Ynew );
   Xcc = Xnew;
   Ycc = Ynew;
}

void Move_to( double x, double y )
{
   Xcc = (int) ( Kx * ( x - Xmin ) );
   Ycc = (int) ( IYw - Ky * ( y - Ymin ) );
}


void Line_to3D( double x, double y, double z )
{
   double X,Y;
   trans( x,y,z, &X,&Y );
   Line_to( X, Y );
}

void Move_to3D( double x, double y, double z )
{
   double X,Y;
   trans( x,y,z, &X,&Y );
   Move_to( X, Y );
}

void Point_set3D( double x, double y, double z )
{
   double X,Y;
   trans( x,y,z, &X,&Y );
   Point_set( X, Y );
}

void Set_Tmat( double x, double y, double z )
{
   double lx,ly,lz,cf,ct,sf,st,w;

   w = sqrt( x * x + y * y + z * z );
   lx = x / w;
   ly = y / w;
   lz = z / w;
   cf = lz;
   sf = sqrt( 1. - cf * cf );
   if( sf > 1.e-7 ) {
      ct = lx / sf;
      st = ly / sf; }
   else {
      ct = 0.;
      st = -1.; }
   T11 = -st;
   T12 = ct;
   T21 = -cf * ct;
   T22 = -cf * st;
   T23 = sf;
   T31 = sf * ct;
   T32 = sf * st;
   T33 = cf;
}

void trans(  double x, double y, double z, double *X, double *Y )
{
   *X = x * T11 + y * T12;
   *Y = x * T21 + y * T22 + z * T23;
}

void trans3(  double x, double y, double z, 
              double *X, double *Y, double *Z  )
{
   *X = x * T11 + y * T12;
   *Y = x * T21 + y * T22 + z * T23;
   *Z = x * T31 + y * T32 + z * T33;
}

void Set_Vwindow3D( double xmin, double ymin, double zmin,
                double xmax, double ymax, double zmax )
{
   double X,Y,Xmin = 1.e5,Xmax = -1.e5,Ymin = 1.e5,Ymax = -1.e5,xw,yw,dw;

   trans( xmin,ymin,zmin, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmax,ymin,zmin, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmin,ymax,zmin, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmax,ymax,zmin, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmin,ymin,zmax, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmax,ymin,zmax, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmin,ymax,zmax, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   trans( xmax,ymax,zmax, &X,&Y );
   if( X > Xmax ) Xmax = X;
   if( Y > Ymax ) Ymax = Y;
   if( X < Xmin ) Xmin = X;
   if( Y < Ymin ) Ymin = Y;
   xw = Xmax - Xmin;
   yw = Ymax - Ymin;
   if( xw / yw > IXw / IYw ) {
      dw = ( xw * IYw / IXw - yw ) / 2.;
      Ymin -= dw;
      Ymax += dw; }
   else {
      dw = ( yw * IXw / IYw -xw ) / 2.;
      Xmin -= dw;
      Xmax += dw; }
   Set_Vwindow( Xmin,Ymin,Xmax,Ymax );
}


