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
const bool secondaryElectrons = false;
const bool ionization = false;
ParticleDataBase3D *pdb3dGlobal;

const double gridthickness = 25e-6;
const double wirethickness = 10e-6;
const double gridsize= 60e-6;

const double RFEAWH = gridsize*10;
const double beamsize = gridsize*7;

const double Etest = 100;

const double testrange = 10;
const double testresolution = 0.25;

//Cross section data
const double ArgonP1=15.8;
const int ArgonQ1=6;
const double ArgonP2=29.2;
const int ArgonQ2=2;
const double Argona=4.0e-14;
const double Argonb=0.62;
const double Argonc=0.4;
const double ArgonN=1e17;
const double ArgonM=40.0;

const double eMass = 1.0/1836.00;

double collectedCharge = 0;
double depth;
double griddisplace[8];

int testNumIonizations = 0;

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



class THCallback : public TrajectoryHandlerCallback {
    Geometry &_geom;
public:

    THCallback( Geometry &geom ) 
        : _geom(geom){
    }

    virtual ~THCallback() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        Particle3D *p3d = (Particle3D *)( particle );
        ParticleP3D *xcur3d = (ParticleP3D *)( xcur );
        ParticleP3D *xend3d = (ParticleP3D *)( xend );
        Vec3D vel = xcur3d->velocity();
        Vec3D pos = xcur3d->location();
        if (p3d->q()<0.0){
            double E= 0.5*p3d->m()*vel.ssqr()/CHARGE_E;
            double U = E/ArgonP1;
            double sum = 0;
            sum += max(ArgonQ1*log(E/ArgonP1)/(E*ArgonP1),0.0);
            sum += max(ArgonQ1*log(E/ArgonP1)/(E*ArgonP1),0.0);
            double crossSection = Argona*(1-Argonb*exp(-Argonc*(U-1)))*sum;
            double meanFreePath = 1/(crossSection*ArgonN);

            double distance = sqrt((xcur3d->location()-xend3d->location()).ssqr());
            double probability = (1-exp(-distance/meanFreePath));
            double random = randomDouble();
            if (random < probability){
                *xend3d = *xcur3d;
                particle->set_status( PARTICLE_COLL );

                Vec3D dir = vel;
                dir.normalize();

                double Ee1 = E-ArgonP1-0.1;
                double speed = sqrt( Ee1*2.0*CHARGE_E/(eMass*MASS_U) );

                pdb3dGlobal->add_particle( p3d->IQ(), -1.0, eMass, ParticleP3D( 0.0, 
                                                                         pos[0], speed * dir[0], 
                                                                         pos[1], speed * dir[1], 
                                                                         pos[2], speed * dir[2] ) );

                double Ee2 = 0.05;
                speed = sqrt( Ee2*2.0*CHARGE_E/(eMass*MASS_U) );

                pdb3dGlobal->add_particle( p3d->IQ(), -1.0, eMass, ParticleP3D( 0.0, 
                                                                         pos[0], speed * dir[0], 
                                                                         pos[1], speed * dir[1], 
                                                                         pos[2], speed * dir[2] ) );

                double Ei = 0.05;
                speed = sqrt( Ei*2.0*CHARGE_E/(ArgonM*MASS_U) );
                pdb3dGlobal->add_particle( p3d->IQ(), 1.0, ArgonM, ParticleP3D( 0.0, 
                                                                         pos[0], speed * dir[0], 
                                                                         pos[1], speed * dir[1], 
                                                                         pos[2], speed * dir[2] ) );

                // testNumIonizations++;
                // cout << testNumIonizations << endl;
                // cout << "Ionization occured!" << endl;
            }
        }
    }
};

class SECallback : public TrajectoryEndCallback {

    Geometry &_geom;
    double    _k;
    double    _theta;

    int getWallNumber(Vec3D loc){ // 1, 2, 3, 4, 5, 6 for the side walls, -1 for internal geometry walls
        Int3D size = _geom.size(); 
        double w, h, d;
        double margin = 1e-6;
        w = (size[0]-1)*_geom.h();
        h = (size[1]-1)*_geom.h();
        d = (size[2]-1)*_geom.h();
        loc = loc - _geom.origo();
        if (loc[0]<margin){
            return 1;
        }
        else if (loc[0]>w-margin){
            return 2;
        }
        else if (loc[1]<margin){
            return 3;
        }
        else if (loc[1]>h-margin){
            return 4;
        }
        else if (loc[2]<margin){
            return 5;
        }
        else if (loc[2]>d-margin){
            return 6;
        }
        else {
            return -1;
        }
    }

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

    void createSecondaryElectron(ParticleBase *particle, class ParticleDataBase *pdb ){
        Particle3D *p3d = (Particle3D *)( particle );
        Vec3D loc = p3d->location();
        //Vec3D vel = p3d->velocity();
        // Make secondaries based on location and energy

        // Get normal
        Vec3D normal = getNormal( loc );
        // Adjust location off the surface
        loc += 0.5*_geom.h()*normal;

        // Randomize velocity and direction
        double mass = 1.0/1836.00;
        double speed = sqrt( 20.0*randomDouble()*CHARGE_E/(mass*MASS_U) );


        Vec3D vel = normal*speed;

        ParticleDataBase3D *pdb3d = (ParticleDataBase3D *)( pdb );
        pdb3d->add_particle( p3d->IQ(), -1.0, mass, ParticleP3D( 0.0, 
                                                                 loc[0], vel[0], 
                                                                 loc[1], vel[1], 
                                                                 loc[2], vel[2] ) );

        
        collectedCharge += p3d->q()/CHARGE_E;
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
        ParticleDataBase3D *pdb3d = (ParticleDataBase3D *)( pdb );
        Vec3D loc = p3d->location();
        Vec3D vel = p3d->velocity();
        double E = 0.5*p3d->m()*vel.ssqr()/CHARGE_E;

        Int3D size = _geom.size(); 
        double w, h;
        w = (size[0]-1)*_geom.h();
        h = (size[1]-1)*_geom.h();
        switch(getWallNumber(loc)){
            case 1:
                if(periodic){
                    pdb3d->add_particle( p3d->IQ(), p3d->q(),p3d->m(), ParticleP3D( 0.0, 
                                                                     loc[0]+w, vel[0], 
                                                                     loc[1], vel[1], 
                                                                     loc[2], vel[2] ) );
                }
            break;
            case 2:
                if(periodic){
                    pdb3d->add_particle( p3d->IQ(), p3d->q(), p3d->m(), ParticleP3D( 0.0, 
                                                                     loc[0]-w, vel[0], 
                                                                     loc[1], vel[1], 
                                                                     loc[2], vel[2] ) );
                }
            break;
            case 3:
                if(periodic){
                    pdb3d->add_particle( p3d->IQ(), p3d->q(), p3d->m(), ParticleP3D( 0.0, 
                                                                     loc[0], vel[0], 
                                                                     loc[1]+h, vel[1], 
                                                                     loc[2], vel[2] ) );
                }
            break;
            case 4:
                if(periodic){
                    pdb3d->add_particle( p3d->IQ(), p3d->q(), p3d->m(), ParticleP3D( 0.0, 
                                                                     loc[0], vel[0], 
                                                                     loc[1]-h, vel[1], 
                                                                     loc[2], vel[2] ) );
                }
            break;
            case 5:
            break;
            case 6:
                collectedCharge += p3d->q()/CHARGE_E;
                if(secondaryElectrons && randomDouble()<0.2 && p3d->q()>0 )
                    createSecondaryElectron(particle, pdb);
            break;
                
            case -1:
                // if(secondaryElectrons)
                //     createSecondaryElectron(particle, pdb);
            break;
        }
    }
};


void simu( int argc, char **argv )
{
    // griddisplace[2] = 20e-6;
    // griddisplace[3] = 20e-6;

    // griddisplace[4] = -20e-6;
    // griddisplace[5] = -20e-6;

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
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vsec) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom , 1e-8);

    EpotField epotconst( geom );
    MeshScalarField scharge( geom );
    solver.solve( epotconst, scharge );

    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,   1.0) );
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,   0.0) );

    geom.build_mesh();
    geom.build_surface();

    EpotField epotvar( geom );

    solver.solve( epotvar, scharge );

    ParticleDataBase3D pdb( geom );
    pdb3dGlobal=&pdb;
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { !periodic, !periodic, !periodic, !periodic, false, false };
    pdb.set_mirror( pmirror );
    SECallback secb( geom );
    pdb.set_trajectory_end_callback( &secb );
    THCallback thcb( geom );
    pdb.set_trajectory_handler_callback( &thcb );

    double current[3][20000];
    int i = 0;
    for (double Vdis=-50; Vdis<=200; Vdis+=10)
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


        double anchor=Vdis>0?Vdis:0;
        double Etest = 0.5;

        while(Etest<400)
        {
            cout << "Vdis: " << Vdis << " Etest: " << Etest << endl;


            pdb.clear(); 

            //Mono-energetic beam
            pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
                                                Etest, 0.0, 0.5, 
                                                Vec3D(0e-6,0e-6,20e-6),
                                                Vec3D(1,0,0), 
                                                Vec3D(0,1,0),
                                                beamsize/2,
                                                beamsize/2);

            // Bimodal distribution
            // double E1 = 20;
            // double v1 = sqrt( E1*2*CHARGE_E/(MASS_U) );

            // pdb.add_rectangular_beam_with_velocity( 50000, 0, 1.0, 1.0, 
            //                                     v1, v1/sqrt(E1)/2.1, 0.0, 
            //                                     Vec3D(0e-6,0e-6, 1e-8),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);

            // double E2 = 40;
            // double v2 = sqrt( E2*2*CHARGE_E/(MASS_U) );

            // pdb.add_rectangular_beam_with_velocity( 50000, 0, 1.0, 1.0, 
            //                                     v2, v2/sqrt(E2)/3, 0.0, 
            //                                     Vec3D(0e-6,0e-6,1e-8),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);

            // Two peak distribution
            // double E1 = 25;

            // pdb.add_rectangular_beam_with_energy( 50000, 0, 1.0, 1.0, 
            //                                     E1, 0.0, 0.0, 
            //                                     Vec3D(0e-6,0e-6, 1e-8),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);

            // double E2 = 50;

            // pdb.add_rectangular_beam_with_energy( 50000, 0, 1.0, 1.0, 
            //                                     E2, 0.0, 0.0, 
            //                                     Vec3D(0e-6,0e-6,1e-8),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);


            // Electrons
            // pdb.add_rectangular_beam_with_energy( 100000, 0, -1.0, eMass, 
            //                                     10, 0.0, 0, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);


            // pdb.add_rectangular_beam_with_energy( 40000, 0, 1.0, 1.0, 
            //                                     50, 10, 0.0, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);




            pdb.iterate_trajectories( scharge, efield, bfield );

            current[0][i]=Vdis;
            current[1][i]=Etest;
            current[2][i]=collectedCharge;
            cout << collectedCharge <<endl;

            collectedCharge=0;


            // Angular distribution influence
            // pdb.clear(); 

            // pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
            //                                     Etest, 0.0, 0.0, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[2][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;

            // pdb.clear(); 

            // pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
            //                                     Etest, 0.0, 0.001, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[3][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;

            // pdb.clear(); 

            // pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
            //                                     Etest, 0.0, 0.01, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[4][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;

            // pdb.clear(); 

            // pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
            //                                     Etest, 0.0, 0.1, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[5][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;

            // pdb.clear(); 

            // pdb.add_rectangular_beam_with_energy( 100000, 0, 1.0, ArgonM, 
            //                                     Etest, 0.0, 0.5, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[6][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;
                        
            // pdb.clear(); 
            // pdb.add_rectangular_beam_with_energy( 60000, 0, 1.0, 1.0, 
            //                                     30, 15, 0.1, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.add_rectangular_beam_with_energy( 40000, 0, 1.0, 1.0, 
            //                                     50, 10, 0.1, 
            //                                     Vec3D(0e-6,0e-6,20e-6),
            //                                     Vec3D(1,0,0), 
            //                                     Vec3D(0,1,0),
            //                                     beamsize/2,
            //                                     beamsize/2);
            // pdb.iterate_trajectories( scharge, efield, bfield );

            // current[2][i]=collectedCharge;
            // cout << collectedCharge <<endl;

            // collectedCharge=0;
            Etest+=clip(fabs((Etest-anchor)/4.0),0.5,10);
            // Etest+= Etest<20?testresolution:2;
            i++;
        }
        if (-50==Vdis){
            geom.save( "geom.dat" );
            epotsum.save( "epot.dat" );
            pdb.save( "pdb.dat" );       
        }
    }
    for(int j=0;j<i;j++)
    {
        cout<<current[0][j] << '\t' << current[1][j] << '\t' << current[2][j] << endl;
    }
    // for(int j=1;j<i;j++){
    //     double x = (current[0][j-1]+current[0][j])/2;
    //     double dx = current[0][j]-current[0][j-1];
    //     double dy1 = current[1][j]-current[1][j-1];
    //     // double dy2 = current[2][j]-current[2][j-1];

    //     cout << x << '\t' << (-dy1/dx) << '\t' << endl;
    // }
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
