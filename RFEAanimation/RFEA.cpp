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
const double Vdis = -5;
const double Vsec = -40;
const double Vcol = -10;

const double spacing [5] = { 100e-6, 200e-6, 200e-6, 200e-6, 200e-6};
const bool periodic = false;

const double gridthickness = 25e-6;
const double wirethickness = 10e-6;
const double gridsize= 60e-6;

const double RFEAWH = gridsize*5;
const double beamsize = gridsize*2;

const double eMass = 1.0/1836.00;

double depth;
double griddisplace[8];

double randomDouble()
{
    return (double)rand() / (double)RAND_MAX ;
}

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

bool grid4( double x, double y, double z )
{
    return gridn( 4, x + griddisplace[6], y + griddisplace[7], z);
}

double clip(double n, double lower, double upper) 
{
    return n <= lower ? lower : n >= upper ? upper : n;
}


void simu( int argc, char **argv )
{
    // griddisplace[2] = 20e-6;
    // griddisplace[3] = 20e-6;

    // griddisplace[4] = 20e-6;
    // griddisplace[5] = -5e-6;

    double h = (10e-6)/3;
    depth = - gridthickness;

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
    // Transformation T;
    // T.rotate_z( 0.78539816339 );
    // s1->set_transformation( T );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( grid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( grid3 );
    geom.set_solid( 9, s3 );
    Solid *s4 = new FuncSolid( grid4 );
    geom.set_solid( 10, s4 );
    
    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  5,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vcol) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vsup) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   0.0) );
    // geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vsec) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom , 1e-8);

    EpotField epotconst( geom );
    MeshScalarField scharge( geom );
    solver.solve( epotconst, scharge );

    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   1.0) );
    // geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,   0.0) );

    geom.build_mesh();
    geom.build_surface();

    EpotField epotvar( geom );

    solver.solve( epotvar, scharge );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { !periodic, !periodic, !periodic, !periodic, false, false };
    pdb.set_mirror( pmirror );

    double current[2][2000];
    int i = 0;
    for (double Vdis=-50; Vdis<= 200; Vdis+=1)
    {  
        EpotField epotsum( geom );
        epotsum += epotvar;
        epotsum *= Vdis;
        epotsum += epotconst;

        EpotEfield efield( epotsum );
        field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
        efield.set_extrapolation( efldextrpl );
        MeshVectorField bfield;

        efield.recalculate();



        cout << "Vdis: " << Vdis << endl;


        pdb.clear(); 

        // Bimodal distribution
        double E1 = 20;
        double v1 = sqrt( E1*2*CHARGE_E/(MASS_U) );

        pdb.clear(); 
        pdb.add_rectangular_beam_with_velocity( 500, 0, 1.0, 1.0, 
                                            v1, v1/sqrt(E1)/2.1, 0.0, 
                                            Vec3D(0e-6,0e-6, 1e-8),
                                            Vec3D(1,0,0), 
                                            Vec3D(0,1,0),
                                            beamsize/2,
                                            beamsize/2);

        double E2 = 40;
        double v2 = sqrt( E2*2*CHARGE_E/(MASS_U) );

        pdb.add_rectangular_beam_with_velocity( 500, 0, 1.0, 1.0, 
                                            v2, v2/sqrt(E2)/3, 0.0, 
                                            Vec3D(0e-6,0e-6,1e-8),
                                            Vec3D(1,0,0), 
                                            Vec3D(0,1,0),
                                            beamsize/2,
                                            beamsize/2);


        pdb.iterate_trajectories( scharge, efield, bfield );

        GeomPlotter geomplotter( geom );
        geomplotter.set_size( 750, 750 );
        geomplotter.set_particle_database( &pdb );
        geomplotter.set_view(VIEW_ZX);
        geomplotter.plot_png( to_string(Vdis) + "plot.png" ); 

        current[0][i]=Vdis;
        current[1][i]=pdb.get_statistics().bound_collisions(6);
        cout << pdb.get_statistics().bound_collisions(6) <<endl;
                    
        i++;
    }
    for(int j=0;j<i;j++)
    {
        cout<<current[0][j] << '\t' << current[1][j] << endl;
    }
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
