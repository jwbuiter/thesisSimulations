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
const double Vdis = 15;
const double Vsec = -40;
const double Vcol = -10;

const double spacing [5] = { 100e-6, 300e-6, 200e-6, 200e-6, 200e-6};

const double gridthickness = 80e-6;
const double wirethickness = 40e-6;
const double gridsize= 150e-6;

const double RFEAWH = gridsize*3;
const double beamsize = gridsize;

double collectedCharge = 0;

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
    return gridn( 1, x, y, z);
}

bool grid2( double x, double y, double z )
{
    return gridn( 2, x, y, z);
}

bool grid3( double x, double y, double z )
{
    return gridn( 3, x, y, z);
}

bool grid4( double x, double y, double z )
{
    return gridn( 4, x, y, z);
}

double randomDouble()
{
    return (double)rand() / (double)RAND_MAX ;
}

class THCallback : public TrajectoryHandlerCallback {
public:

    THCallback() {}

    virtual ~THCallback() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) const {
        double random = randomDouble();
    }
};

class SECallback : public TrajectoryEndCallback {

    Geometry &_geom;
    double    _k;
    double    _theta;

    Vec3D getNormal(Vec3D loc){
        Int3D size = _geom.size(); 
        double w, h, d;
        double margin = 1e-6;
        w = (size[0]-1)*_geom.h();
        h = (size[1]-1)*_geom.h();
        d = (size[2]-1)*_geom.h();
        loc = loc - _geom.origo();
        if (loc[0]<margin){
            return Vec3D(1, 0, 0);
        }
        else if (loc[0]>w-margin){
            return Vec3D(-1, 0, 0);
        }
        else if (loc[1]<margin){
            return Vec3D(0, 1, 0);
        }
        else if (loc[1]>h-margin){
            return Vec3D(0, -1, 0);
        }
        else if (loc[2]<margin){
            return Vec3D(0, 0, 1);
        }
        else if (loc[2]>d-margin){
            return Vec3D(0, 0, -1);
        }
        else {
            return _geom.surface_normal( loc + _geom.origo());
        }
    }
    
public:

    SECallback( Geometry &geom ) 
        : _geom(geom){

        _k = 9.0;
        _theta = 0.5;
    }

    virtual ~SECallback() {}

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) {

        Particle3D *p3d = (Particle3D *)( particle );
        Vec3D loc = p3d->location();
        Vec3D vel = p3d->velocity();
        double E = 0.5*p3d->m()*vel.ssqr()/CHARGE_E;


        if (loc[2]>(1320e-6 - 1e-6)){
            collectedCharge += p3d->q()/CHARGE_E;
        }
       // Make secondaries based on location and energy
        if( p3d->q()>0 ) {

            // Get normal
            Vec3D normal = getNormal( loc );
            // Adjust location off the surface
            loc += 0.02*_geom.h()*normal;

            // Randomize velocity and direction
            double mass = 1.0/1836.00;
            double speed = sqrt( 20.0*randomDouble()*CHARGE_E/(mass*MASS_U) );

            // Find tangents
            Vec3D tang1 = normal.arb_perpendicular();
            Vec3D tang2 = cross( normal, tang1 );
            tang1.normalize();
            tang2.normalize();

            // Build velp in natural coordinates
            double azm_angle = 2.0*M_PI*randomDouble();
            double pol_angle = asin( sqrt(randomDouble()) );
            Vec3D velp( speed*cos(pol_angle),  speed*sin(pol_angle), 0.0 );
            Transformation t;
            t.rotate_x( azm_angle );
            velp = t.transform_vector( velp );

            // Convert to surface coordinates
            Vec3D vel = velp[0]*normal + velp[1]*tang1 + velp[2]*tang2;

            //Vec3D vel = normal*speed;

            ParticleDataBase3D *pdb3d = (ParticleDataBase3D *)( pdb );
            pdb3d->add_particle( p3d->IQ(), -1.0, mass, ParticleP3D( 0.0, 
                                                                     loc[0], vel[0], 
                                                                     loc[1], vel[1], 
                                                                     loc[2], vel[2] ) );

            if (loc[2]>(1320e-6 - 1e-6)){
                collectedCharge += p3d->q()/CHARGE_E;
            }
        }
    }
};


void simu( int argc, char **argv )
{
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
    /*Transformation T;
    T.rotate_z( 0.485763697 );
    s1->set_transformation( T );*/
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
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vdis) );
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vsec) );

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
    bool pmirror[6] = { true, true, true, true, false, false };
    pdb.set_mirror( pmirror );
    SECallback secb( geom );
    pdb.set_trajectory_end_callback( &secb );
  
    solver.solve( epot, scharge_ave );
    efield.recalculate();
    pdb.clear(); 
    // pdb.add_cylindrical_beam_with_energy( 100000, 0, 1.0, 1.0, 
    //                                     12, 0.0, 0.0, 
    //                                     Vec3D(0e-6,0e-6,20e-6),
    //                                     Vec3D(1,0,0), 
    //                                     Vec3D(0,1,0),
    //                                     beamsize/2);
    pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, 1.0, 
                                        10, 0.0, 0.1, 
                                        Vec3D(0e-6,0e-6,20e-6),
                                        Vec3D(1,0,0), 
                                        Vec3D(0,1,0),
                                        beamsize/2,
                                        beamsize/2);
    pdb.iterate_trajectories( scharge, efield, bfield );
                
    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
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
