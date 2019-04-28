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
const double Vsec = -40;

const double spacing [5] = { 100e-6, 300e-6, 200e-6, 200e-6, 200e-6};

const double gridthickness = 80e-6;
const double wirethickness = 40e-6;
const double gridsize= 150e-6;

const double RFEAWH = gridsize*3;
const double beamsize = gridsize;

// bool grid(double x, double y, double spacing, double thickness, double xoffset, double yoffset)
// {
//     double xpos=(x+xoffset+50*spacing + spacing/2 + thickness/2)/spacing;
//     double ypos=(y+yoffset+50*spacing + spacing/2 + thickness/2)/spacing;
//     bool xgrid=(((xpos-(long)(xpos))*spacing)<thickness);
//     bool ygrid=(((ypos-(long)(ypos))*spacing)<thickness);
//     return (xgrid || ygrid);
// }


bool planen( int n, double x, double y, double z )
{
    double planepos = 0;
    for(int i = 0; i < n; i++)
    {
        planepos += spacing[i] + gridthickness;
    }
    return (z==planepos);
}

bool plane1( double x, double y, double z )
{
    return planen( 1, x, y, z);
}

// bool grid2( double x, double y, double z )
// {
//     return gridn( 2, x, y, z);
// }

// bool grid3( double x, double y, double z )
// {
//     return gridn( 3, x, y, z);
// }

// bool grid4( double x, double y, double z )
// {
//     return gridn( 4, x, y, z);
// }

double randomDouble()
{
    return (double)rand() / (double)RAND_MAX ;
}

void simu( int argc, char **argv )
{
    double h = (10e-6)/2;
    double depth = - gridthickness;

    for(int i = 0; i < (sizeof(spacing)/sizeof(*spacing)); i++)
    {
        depth += spacing[i] + gridthickness;
    }

    depth = 100e-6;
    double sizereq[3] = { RFEAWH,
                          RFEAWH, 
                          depth};

    Int3D meshsize( (int)round(sizereq[0]/h)+1,
                    (int)round(sizereq[1]/h)+1,
                    (int)round(sizereq[2]/h)+1 );

    Vec3D origo( -RFEAWH/2, -RFEAWH/2, 0);
    Geometry geom( MODE_3D, meshsize, origo, h );


    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   1.0) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom, 1e-8 );

    EpotField epotconst( geom );
    MeshScalarField scharge( geom );
    solver.solve( epotconst, scharge );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { true, true, true, true, false, false };
    pdb.set_mirror( pmirror );

    double current[2][2000];
    int i = 0;
    for (double Vdis=-50; Vdis<=200; Vdis+=1)
    {  
        EpotField epotvar( geom );
        epotvar+= epotconst;
        epotvar *= Vdis;

        EpotField epot( geom );
        MeshScalarField scharge( geom );

        EpotEfield efield( epotvar );
        field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
        efield.set_extrapolation( efldextrpl );
        MeshVectorField bfield;

        efield.recalculate();

        cout << "Vdis: " << Vdis << endl;

        int total = 100000;
        int N1 = 50000;
        double E1=20;
        double v1 = sqrt( E1*2*CHARGE_E/(MASS_U) );

        pdb.clear(); 
        pdb.add_rectangular_beam_with_velocity( N1, 0, 1.0, 1.0, 
                                            v1, v1/sqrt(E1)/2.1, 0.0, 
                                            Vec3D(0e-6,0e-6, 1e-8),
                                            Vec3D(1,0,0), 
                                            Vec3D(0,1,0),
                                            beamsize/2,
                                            beamsize/2);

        double E2 = 40;
        double v2 = sqrt( E2*2*CHARGE_E/(MASS_U) );

        pdb.add_rectangular_beam_with_velocity( total - N1, 0, 1.0, 1.0, 
                                            v2, v2/sqrt(E2)/3, 0.0, 
                                            Vec3D(0e-6,0e-6,1e-8),
                                            Vec3D(1,0,0), 
                                            Vec3D(0,1,0),
                                            beamsize/2,
                                            beamsize/2);

        pdb.iterate_trajectories( scharge, efield, bfield );

        current[0][i]=Vdis;
        current[1][i]=pdb.get_statistics().bound_collisions(6);
        cout << pdb.get_statistics().bound_collisions(6) <<endl;

        i++;

        if (Vdis==0.0){
            geom.save( "geom.dat" );
            epotvar.save( "epot.dat" );
            pdb.save( "pdb.dat" );       
        }
    }
    for(int j=0;j<i;j++)
    {
        cout<<current[0][j] << '\t' << current[1][j] << endl;
    }
    for(int j=1;j<i;j++){
        double x = (current[0][j-1]+current[0][j])/2;
        double dx = current[0][j]-current[0][j-1];
        double dy1 = current[1][j]-current[1][j-1];

        cout << x << '\t' << (-dy1/dx) << '\t' << endl;
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
