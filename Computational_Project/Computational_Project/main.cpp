// This code is writen by : Maedeh Seyedolmohadesin
// Computational Physics Project Fall 2018, Northeastern University, Boston, MA
// This code will simulate a monoatomic Lennard jones system consists of N atoms in an (N,V,E) ensemble.
// BCC Lattice is used for the initial positions of the atoms
// Reduced unit is used here
// RDF stands for Radial Distribution Function

// **************************************************** Include Libraries ****************************************************

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>

using namespace std;

//********************************************************* Settings *********************************************************

//system settings

int N0 =7 ; // number of one type atom in each direction (BCC lattic)
int N =  2 * pow (N0 , 3); // total number of atoms
double rho=1.2; // density of the system
double sigma1=1.0; // Lennard jones parameter
double mass1=1.0;
double sigma2=1.0;
double mass2=1.0;
double K=1.0; //Boltzman const.
double epsilon=1.0; // Lennard jones parameter
double d=pow((N/(rho)),(1.00/3.00)); // length of box in each dimension

//MD settings

int nsteps=100000; // number of steps to run the simulation
double dt=0.001; // time element
int nstxout=100; // number of steps that elapse between writing coordinates and energies to output files
double force_cut=3.5 ; // after this radius the force won't be calulated
double T=1.2; // initial temperature

// RDF setting

int nstrdf = 100; // number of steps that elapse between writing the RDF outputs
int rdfpoints= 1000; // Number of data points on RDF plot

// Energy  variables

double E=0 ; // total energy of the system
double KE=0 ; // total kinetic energy of the system
double U=0 ; // total potential energy of the system

// RDF Variables

double binsize = (d / 2) / rdfpoints ;
double *bin = new double[rdfpoints] ;

// Anderson Thermostat

double temp=1.2; // Temperature
int nsttemp=100000000; // number of steps that elapse between reinitializing the velocities
// ******************************************************* class Atom **********************************************************

class Atom
{
public:
    
    double x;
    double y;        // positions
    double z;
    double vx;
    double vy;       // velocities
    double vz;
    double fx;
    double fy;       // forces
    double fz;
    double u;       // potential energy
    double ke;      // kinetic energy
    
};

//****************************************************** List of Functions ******************************************************

Atom* gen_atom ();                                        // Generates objects for class Atom

void init_pos( Atom *p );                                 // Initializes the positions

double normal_dist (double sigma);                        // generates numbers from normal dist. with variance sigma^2 and mean 0

void init_vel( Atom *p );                                 // Initializes the velocities

Atom* initialize ();                                      // Initializes the system

void calculate_forces (Atom *p);                          // Calculates forces on atoms

void MD_run (Atom* p);                                    // Runs the simulations and does the calculations

void write_trj (ofstream &file, Atom*p, int n);           // Writes the trajectories of atoms Atom*p into 'file'

void check_boundry (Atom* p, int i);                      // Apllies the periodic boundry condition

//******************************************************* main() function *******************************************************

int main()
{
    srand(time(NULL));
    
    Atom *atom_pointer; // defining pointer to store class object's addresses
    
    atom_pointer=initialize();
    
    MD_run(atom_pointer);
    
    // ******** test *******
    
    return 0;
}

//***************************************************** Atom* gen_atom () *******************************************************

Atom* gen_atom ()
{
    Atom *p=new Atom [N];
    for (int i=0; i< N; i++)
    {
        Atom atom;
        p[i]=atom;
    }
    return p;
}

//************************************************** void init_pos( Atom *p ) ***************************************************

void init_pos( Atom *p )
{
    double r=d/N0; // distance between atoms in each direction
    double n=cbrt(N/2);
    int counter=0;
    
    
    for (int k=0; k<n; k++)
    {
        for (int j=0; j<n;j++)
        {
            for (int i=0; i<n; i++)
            {
                p[counter].x=i*r;
                p[counter].y=j*r;                 // first simple cubic lattice
                p[counter].z=k*r;
                
                p[counter + N/2].x=i*r + r/2;
                p[counter + N/2].y=j*r + r/2;     // second simple cubic lattice
                p[counter + N/2].z=k*r + r/2;
                
                counter++;
                
            }
        }
    }
}

//********************************************* double normal_dist (double sigma) ***********************************************

double normal_dist (double sigma)
{
    
    double pi =3.14159;
    double a=rand() * (1.0 / RAND_MAX);
    double b=rand() * (1.0 / RAND_MAX);
    double r;
    r=sqrt(-2.0 * log(a)) * cos(2.0 * pi * b);
    return r*sigma;
}

//*************************************************** void init_vel(Atom *p) ****************************************************

void init_vel(Atom *p)
{
    double Px=0;
    double Py=0;
    double Pz=0;
    
    
    for (int i=0; i< N/2 ; i++)
    {
        p[i].vx=normal_dist(sqrt(K * T / mass1));
        p[i].vy=normal_dist(sqrt(K * T / mass1));
        p[i].vz=normal_dist(sqrt(K * T / mass1));
                                                        // Initialize velocities based on M-B dist.
        p[i+N/2].vx=normal_dist(sqrt(K * T / mass2));
        p[i+N/2].vy=normal_dist(sqrt(K * T / mass2));
        p[i+N/2].vz=normal_dist(sqrt(K * T / mass2));
    }
    
    for (int j=0; j<N/2; j++)
    {
        Px += mass1 * p[j].vx;
        Py += mass1 * p[j].vy;
        Pz += mass1 * p[j].vz;
        
        Px += mass2 * p[j + N/2].vx;
        Py += mass2 * p[j + N/2].vy;
        Pz += mass2 * p[j + N/2].vz;
        
    }
    
    for (int j=0; j<N/2; j++)
    {
        p[j].vx -= Px/(mass1 * N);
        p[j].vy -= Py/(mass1 * N);
        p[j].vz -= Pz/(mass1 * N);
                                           // remove center of mass velocity (no translational motion)
        p[j + N/2].vx -= Px/(mass2 * N);
        p[j + N/2].vy -= Py/(mass2 * N);
        p[j + N/2].vz -= Pz/(mass2 * N);
        
    }
}

//***************************************************** Atom* initialize () *****************************************************

Atom* initialize ()
{
    Atom* p= gen_atom();
    init_pos(p);
    init_vel(p);
    
    return p;
}

//*********************************************** void calculate_forces (Atom *p) ***********************************************

void calculate_forces (Atom *p)
{
    double distance;
    double xdistance;
    double ydistance;   // Distances between pair of atoms
    double zdistance;
    double f;
    double fx;
    double fy;
    double fz;
    double ratio; // sigma/r term in L-J potential
    double u; // potential energy
    
    for (int i=0; i < N; i++)
    {
        p[i].fx= 0 ;
        p[i].fy= 0 ;
        p[i].fz= 0 ;
        p[i].u=0 ;
    }
    
    for (int i=0; i<rdfpoints; i++)
    {
        bin[i] = 0 ;
    }
    
    for (int i=0; i<N; i++)
    {
        for (int j=i+1; j<N; j++)
        {
            xdistance=p[i].x-p[j].x;
            ydistance=p[i].y-p[j].y;
            zdistance=p[i].z-p[j].z;
            
            if (xdistance < (-d/2) )
            {
                xdistance=d+xdistance;
            }
            
            if (xdistance >= (d/2) )
            {
                xdistance=xdistance-d;
            }
            
            if (ydistance < (-d/2) )
            {
                ydistance=d+ydistance;             // periodic boundry condition
            }
            
            if (ydistance >= (d/2) )
            {
                ydistance=ydistance-d;
            }
            if (zdistance < (-d/2) )
            {
                zdistance=d+zdistance;
            }
            
            if (zdistance >= (d/2) )
            {
                zdistance=zdistance-d;
            }
            
            distance=sqrt((xdistance*xdistance) + (ydistance*ydistance) + (zdistance*zdistance));
            
            ratio=sigma1/distance;
            
            if (distance > force_cut)
            {
                f=0;
            }
            else
            {
                f=((12*epsilon)/(sigma1 * sigma1)) * (pow(ratio,14.0)-pow(ratio,8));
            }
            
            fx=f*xdistance;
            fy=f*ydistance;
            fz=f*zdistance;
            
            p[i].fx += fx ;
            p[i].fy += fy ;
            p[i].fz += fz ;
            
            p[j].fx -= fx ;
            p[j].fy -= fy ;
            p[j].fz -= fz ;
            
            u = epsilon * (pow (ratio, 12) - 2 * pow (ratio, 6));  // potential energy calculation
            p[i].u += u;
            p[j].u += u;
            
            U += u; // Total potential energy
            
            if (distance <= d/2)
            {
                bin [int(distance/binsize)]++; // RDF calculations
            }
        }
    }
}

//*************************************** void write_trj (ofstream &file, Atom*p, int n) ****************************************

void write_trj (ofstream &file, Atom*p, int n) // n is the step number

{
    file << N << "\n" << "step :" << n << "\n" ;
    
    for (int i=0; i<N/2; i++ )
    {
        
        file << "C" << "\t" << p[i].x << "\t" << p[i].y <<  "\t" << p[i].z << "\n";
    }
    
    for (int i=N/2; i<N; i++ )
    {
        
        file << "O" << "\t" << p[i].x << "\t" << p[i].y <<  "\t" << p[i].z << "\n";
    }
    
}

//**************************************************** void MD_run(Atom* p) *****************************************************

void MD_run(Atom* p)

{
    double vx;
    double vy;
    double vz;
    
    double Px=0;
    double Py=0; // Center of mass momentum
    double Pz=0;
    
    ofstream trj; //trajectory file
    ofstream energies; // file
    ofstream RDFfile ; // RDF file
    
    energies.open("/Users/MaedeMohadesin/University/Fall_2018/Computational_Physics/Project/outputs/energies.txt"); // adresses
    trj.open("/Users/MaedeMohadesin/University/Fall_2018/Computational_Physics/Project/outputs/trajectory.xyz"); // to save
    RDFfile.open ("/Users/MaedeMohadesin/University/Fall_2018/Computational_Physics/Project/outputs/RDF.txt"); // the files
    
    energies << "step \t potential energy \t kinetic energy \t total energy \t temperature \n"  ;
    write_trj(trj, p, 1); // first frame
    RDFfile << "r\tg(r)\n " ;
    
    for (int i=0; i<nsteps; i++)
    {
        U=0;
        KE=0;
        E=0;
        Px=0;
        Py=0;
        Pz=0;
        
        calculate_forces(p);
        
        if ( i % nsttemp == 0 && i!=0)   // Anderson Thermostat
            
        {
            for (int j=0; j<N/2 ; j++)
            {
                p[j].vx=normal_dist(sqrt(K * temp / mass1));
                p[j].vy=normal_dist(sqrt(K * temp / mass1));
                p[j].vz=normal_dist(sqrt(K * temp / mass1));
                                                                      // Reinitializing the velocities
                p[j + N/2].vx=normal_dist(sqrt(K * temp / mass2));
                p[j + N/2].vy=normal_dist(sqrt(K * temp / mass2));
                p[j + N/2].vz=normal_dist(sqrt(K * temp / mass2));
            }
            
            for (int j=0; j<N/2; j++)
            {
                Px += mass1 * p[j].vx;
                Py += mass1 * p[j].vy;
                Pz += mass1 * p[j].vz;
                                                 // Calculating  center of mass velocity
                Px += mass2 * p[j + N/2].vx;
                Py += mass2 * p[j + N/2].vy;
                Pz += mass2 * p[j + N/2].vz;
            }
            
            for (int j=0; j<N/2; j++)
            {
                p[j].vx -= Px/(mass1 * N);
                p[j].vy -= Py/(mass1 * N);
                p[j].vz -= Pz/(mass1 * N);
                
                p[j + N/2].vx -= Px/(mass2 * N);
                p[j + N/2].vy -= Py/(mass2 * N);
                p[j + N/2].vz -= Pz/(mass2 * N);     // Removing center of mass velocity and updating positions
                
                p[j].x  += dt * p[j].vx;
                p[j].y  += dt * p[j].vy;
                p[j].z  += dt * p[j].vz;
                
                p[j + N/2].x  += dt * p[j + N/2].vx;
                p[j + N/2].y  += dt * p[j + N/2].vy;
                p[j + N/2].z  += dt * p[j + N/2].vz;
                
                // Kinetic energy calculation
                
                p[j].ke = ( mass1 / 2 ) * ( p[j].vx * p[j].vx + p[j].vy * p[j].vy + p[j].vz * p[j].vz );
                
                p[j + N/2].ke = ( mass2 / 2 ) * ( p[j + N/2].vx * p[j + N/2].vx + p[j + N/2].vy * p[j + N/2].vy + p[j + N/2].vz * p[j + N/2].vz );
                
                KE += p[j].ke;
                KE += p[j + N/2].ke;
                
            }
            
        } // end of anderson thermostat
        
        // Leap-frog algorithm : updating positions and velocities
        
        else
        {
            
            // Atom type 1
            for (int j=0; j<N/2; j++)
            {
                
                vx = p[j].vx ;
                vy = p[j].vy ;
                vz = p[j].vz ;
                
                p[j].vx += dt * p[j].fx / mass1;
                p[j].vy += dt * p[j].fy / mass1;
                p[j].vz += dt * p[j].fz / mass1;
                
                
                p[j].x  += dt * p[j].vx;
                p[j].y  += dt * p[j].vy;
                p[j].z  += dt * p[j].vz;
                
                check_boundry(p, j);
                
                // kinetic energy calculation
                
                vx = (vx + p[j].vx ) / 2;
                vy = (vy + p[j].vy ) / 2;
                vz = (vz + p[j].vz ) / 2;
                
                p[j].ke = ( mass1 / 2 ) * ( vx*vx + vy*vy + vz*vz );
                
                KE += p[j].ke ;
                
            }
            
            // Atom type 2
            
            for (int j=N/2; j<N; j++)
            {
                
                vx = p[j].vx ;
                vy = p[j].vy ;
                vz = p[j].vz ;
                
                p[j].vx += dt * p[j].fx / mass2;
                p[j].vy += dt * p[j].fy / mass2;
                p[j].vz += dt * p[j].fz / mass2;
                
                p[j].x  += dt * p[j].vx;
                p[j].y  += dt * p[j].vy;
                p[j].z  += dt * p[j].vz;
                
                check_boundry(p, j);
                
                // kinetic energy calculation
                
                vx = (vx + p[j].vx ) / 2;
                vy = (vy + p[j].vy ) / 2;
                vz = (vz + p[j].vz ) / 2;
                
                p[j].ke = ( mass2 / 2 ) * ( vx*vx + vy*vy + vz*vz );
                
                KE += p[j].ke ;
                
            }
            
        } // end of leap-frog
        
        // writting trajectories and energies
        
        if (i==0)
        {
            E = KE + U ;
            energies << i+1 << "\t" << U << "\t" << KE << "\t" << E  << "\t" << 2 * KE / (3 * N * K) <<endl;
        }
        
        if ( (i+1) % nstxout == 0  )
        {
            E = KE + U ;
            
            write_trj(trj, p, i+1);
            
            energies << i+1 << "\t" << U << "\t" << KE << "\t" << E  << "\t" << 2 * KE / (3 * N * K) << endl;
        }
        
        // Writting RDF
        
        if ( (i+1) % nstrdf == 0 || i==0)
        {
            for (int j=0; j<rdfpoints; j++)
            {
                double factor;
                double r = (j+1) * binsize;
                
                factor = rho * 4 * 3.14 * r * r * binsize * (N / 2);
                bin[j] =  bin[j]  / factor ;   // normalizing bin []
                RDFfile << r << "\t" << bin[j] << "\n";
            }
            
            RDFfile << "\n";
        }
        
        cout << "n="<< i << endl;
    }
    
    trj.close();
    energies.close();
    RDFfile.close();
}

//********************************************* void check_boundry(Atom* p, int i) **********************************************

void check_boundry(Atom* p, int i)
{
    while ( p[i].x >= d)
    {
        p[i].x = p[i].x - d;
    }
    
    while ( p[i].x < 0 )
    {
        p[i].x = p[i].x + d;
    }
    
    // y
    
    while ( p[i].y >= d)
    {
        p[i].y = p[i].y - d;
    }
    
    while ( p[i].y < 0 )
    {
        p[i].y = p[i].y + d;
    }
    
    // z
    
    while ( p[i].z >= d)
    {
        p[i].z = p[i].z - d;
    }
    
    while ( p[i].z < 0 )
    {
        p[i].z = p[i].z + d;
    }
}


