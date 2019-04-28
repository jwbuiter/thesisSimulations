#include <fstream>
#include <iomanip>
#include <limits>
#include <stdlib.h>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"

using namespace std;

const double Vsup = -60;
const double Vcol = -10;
const double Vdis = 20;
const double Vsec = -40;

const double spacing [5] = { 100e-6, 200e-6, 200e-6, 200e-6};

const double gridthickness = 25e-6;
const double wirethickness = 10e-6;
const double gridsize= 60e-6;

const double RFEAWH = gridsize*10;
const double beamsize = gridsize*5;

double griddisplace[8];

bool grid(double x, double y, double spacing, double thickness, double xoffset, double yoffset)
{
    double xpos=(x+xoffset+50*spacing + spacing/2 + thickness/2)/spacing;
    double ypos=(y+yoffset+50*spacing + spacing/2 + thickness/2)/spacing;
    bool xgrid=(((xpos-(long)(xpos))*spacing)<thickness);
    bool ygrid=(((ypos-(long)(ypos))*spacing)<thickness);
    return (xgrid || ygrid);
}


bool gridn( int n, double x, double y, double z )
{
    double gridend = 0;
    for(int i = 0; i < n; i++)
    {
        gridend += spacing[i] + gridthickness;
    }
    return (z < gridend && z > gridend - gridthickness && grid(x, y, gridsize, wirethickness, 0e-6, 0e-6));
}

bool grid1( double x, double y, double z )
{
    return gridn( 1, x + griddisplace[0], y + griddisplace[1], z);
}

bool grid2( double x, double y, double z )
{
    return gridn( 2, x + griddisplace[2], y + griddisplace[3], z);
}

bool grid3( double x, double y, double z )
{
    return gridn( 3, x + griddisplace[4], y + griddisplace[5], z);
}

void simu( int argc, char **argv )
{
    griddisplace[4]=5e-6;
    griddisplace[5]=20e-6;

    double h = (10e-6)/2;
    double depth = - gridthickness;

    for(int i = 0; i < (sizeof(spacing)/sizeof(*spacing)); i++)
    {
        depth += spacing[i] + gridthickness;
    }

    double sizereq[3] = { RFEAWH,
                          RFEAWH, 
                          depth};

    Int3D meshsize( (int)round(sizereq[0]/h)+1,
                    (int)round(sizereq[1]/h)+1,
                    (int)round(sizereq[2]/h)+1 );

    Vec3D origo( -RFEAWH/2, -RFEAWH/2, 0);
    Geometry geom( MODE_3D, meshsize, origo, h );

    Solid *s1 = new FuncSolid( grid1 );
    Transformation T;
    T.rotate_z( 0.643501109 );
    s1->set_transformation( T );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( grid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( grid3 );
    geom.set_solid( 9, s3 );
    /*Solid *s4 = new FuncSolid( grid4 );
    geom.set_solid( 10, s4 );*/

    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  5,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vcol) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vsup) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   0.0) );
    //geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vsec) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom , 1e-7);

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    solver.solve( epot, scharge);

    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   1.0) );

    geom.build_mesh();
    geom.build_surface();

    EpotField epot2( geom );

    solver.solve( epot2, scharge );

    EpotField epottotal( geom );
    epottotal += epot2;
    epottotal *= Vdis;
    epottotal += epot;

    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vcol) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vsup) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vdis) );

    geom.build_mesh();
    geom.build_surface();

    EpotField epot3( geom );

    solver.solve( epot3, scharge);

    EpotField epotcompare( geom );
    epotcompare += epottotal;
    epotcompare *= -1;
    epotcompare += epot3;

    ParticleDataBase3D pdb( geom );
    pdb.clear(); 

    geom.save( "geom.dat" );
    epotcompare.save( "epot.dat" );
    pdb.save( "pdb.dat" );       
}


int main( int argc, char **argv )
{
    try {
        //ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
        ibsimu.set_thread_count( 4 );
        simu( argc, argv );
    } catch( Error e ) {
        e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
