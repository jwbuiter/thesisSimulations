#include <fstream>
#include <iomanip>
#include <limits>
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

const double gridspacing = 200e-6;
const double initialspacing = 100e-6;
const double collectorspacing = 200e-6;
const double gridthickness = 10e-6;
const double gridsize= 50e-6/3*2;

const double RFEAWH = gridsize*17*1.75;
const double beamsize = 17*gridsize;

const double energyresolution = 1.5;

bool grid(double x, double y, double spacing, double thickness, double xoffset, double yoffset)
{
    double xpos=(x+xoffset+50*spacing + spacing/2 + thickness/2)/spacing;
    double ypos=(y+yoffset+50*spacing + spacing/2 + thickness/2)/spacing;
    bool xgrid=(((xpos-(long)(xpos))*spacing)<thickness);
    bool ygrid=(((ypos-(long)(ypos))*spacing)<thickness);
    return (xgrid || ygrid);
}


bool grid1( double x, double y, double z )
{
    return(z>initialspacing && z<initialspacing+gridthickness) && grid(x, y, gridsize, gridthickness, 0e-6, 0e-6);
}

bool grid2( double x, double y, double z )
{
    return(z>initialspacing+gridspacing && z<initialspacing+gridspacing+gridthickness) && grid(x, y, gridsize, gridthickness, 0e-6, 0e-6);
}

bool grid3( double x, double y, double z )
{
    return(z>initialspacing+gridspacing*2 && z<initialspacing+gridspacing*2+gridthickness) && grid(x, y, gridsize, gridthickness, 5e-6, 15e-6);
}

double clip(double n, double lower, double upper) 
{
    return n <= lower ? lower : n >= upper ? upper : n;
}

double abs(double n)
{
    return n<0?-n:n;
}

void simu( int argc, char **argv )
{
    double h = (10e-6)/3;
    double sizereq[3] = { RFEAWH,
                          RFEAWH, 
                          initialspacing+2*gridspacing+collectorspacing};
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( -RFEAWH/2, -RFEAWH/2, 0);
    Geometry geom( MODE_3D, meshsize, origo, h );

    Solid *s1 = new FuncSolid( grid1 );
    Transformation T;
    T.rotate_z( 0.485763697 );
    s1->set_transformation( T );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( grid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( grid3 );
    geom.set_solid( 9, s3 );

    double current[3][2000];
    int i = 0;
    for (double Vdis=-0.5; Vdis<=0.5; Vdis+=1)
    {
        geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
        geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
        geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
        geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
        geom.set_boundary(  5,  Bound(BOUND_NEUMANN,     0.0) );
        geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vcol) );

        geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
        geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vsup) );
        geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vdis) );

        geom.build_mesh();
        geom.build_surface();

        EpotBiCGSTABSolver solver( geom );

        EpotField epot( geom );
        MeshScalarField scharge( geom );
        MeshScalarField scharge_ave( geom );

        EpotEfield efield( epot );
        field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                         FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                                         FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
        efield.set_extrapolation( efldextrpl );

        MeshVectorField bfield;

        ParticleDataBase3D pdb( geom );
        pdb.set_max_steps( 1000 );
        // bool pmirror[6] = { true, true, true, true, false, false };
        bool pmirror[6] = { false, false, false, false, false, false };
        pdb.set_mirror( pmirror );

        solver.solve( epot, scharge_ave );
        efield.recalculate();

        double anchor=Vdis>0?Vdis:0;
        double E=0;
        
        while(E<15)
        {
            pdb.clear(); 
            // pdb.add_cylindrical_beam_with_energy( 100000, 0, 1.0, 1.0, 
            //                                 Vdis-E, 0.0, 0.0, 
            //                                 Vec3D(0e-6,0e-6,20e-6),
            //                                 Vec3D(1,0,0), 
            //                                 Vec3D(0,1,0),
            //                                 beamsize/2);
            pdb.add_rectangular_beam_with_energy( 100000, 0.0, 1.0, 1.0, 
                                                E, 0.0, 0.0, 
                                                Vec3D(0,0,20e-6),
                                                Vec3D(1,0,0), 
                                                Vec3D(0,1,0),
                                                beamsize/2,
                                                beamsize/2);
            pdb.iterate_trajectories( scharge, efield, bfield );

            current[0][i]=E;
            current[1][i]=Vdis;
            current[2][i]=pdb.get_statistics().bound_collisions(6);
            i++;
            cout<<E<<endl;
            E+=clip(abs((E-anchor)/4),0.02,1);
        }
        
        
    }
    
    for(int j=0;j<i;j++){
        cout<<current[0][j]<<'\t'<<current[1][j]<<'\t'<<current[2][j]<<endl;
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
