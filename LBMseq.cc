#Descritpion: Simple D2Q9 lattice boltzman
#             around an abstacle
#author: Valerio Mazzone
#

#include <iostream>

const unsigned int scale=1;
const unsigned int NX = 32*scale;
const unsigned int NY = NX;
const unsigned dir = 9;
const size_t mem_size_ndir = sizeof(double)*NX*NY*ndir;
const size_t mem_size_scalar = sizeof(double)*NX*NY;
const double w0 = 4.0/9.0; // zero weight
const double ws = 1.0/9.0; // adjacent weight
const double wd = 1.0/36.0; // diagonal weight
const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
const int dirx[] = {0,1,0,-1,0,1,-1,-1,1};
const int diry[] = {0,0,1,0,-1,1,1,-1,-1};
const double nu = 1.0/1.6;
const double tau = 3.0*nu+0.5;
const double u_max = 0.04/scale;
const double rho0 = 1.0;
const unsigned int NSTEPS = 200*scale*scale;


int main(int argc, char* argv[])
{
	double *f1 = (double*) malloc(mem_size_ndir);
	double *f2 = (double*) malloc(mem_size_ndir);
	double *rho = (double*) malloc(mem_size_scalar);
	double *ux = (double*) malloc(mem_size_scalar);
	double *uy = (double*) malloc(mem_size_scalar);

	// Init equilibrium
	init_equilibrium(f1,rho,ux,uy);

	// Main simulation loop; take NSTEPS time steps.
	for(unsigned int n = 0; n < NSTEP; ++n)
	{
		// Stream from f1 storing to f2
		stream(f1,f2);
		// Calculate post streaming density and velocity
		compute_rho_u(f1,rho,ux,uy);
		// Perform collision on f2
		collide(f2,rho,ux,uy);

		// Swap pointers
		double *temp = f1;
		f1 = f2;
		f2 = temp;
	}

	// Deallocate memory
	free(f1); free(f2);
	free(rho); free(ux); free(uy);

	return 0;
}


void init_equilibrium(double *f, double*r, double *u, double *v)
{
  for(unsigned int y = 0; y < NY; ++y)
    {
      for(unsigned int x = 0; x < NX; ++x)
	{
	  double rho = r[scalar_index(x,y)];
	  double ux = u[scalar_index(x,y)];
	  double uy = v[scalar_index(x,y)];
	  
	  for(unsigned int i = 0; i < ndir; ++i)
	    {
	      double cidotu =  dirx[i]*ux + diry[i]*uy;
	      f[field_index(x,y,i)] = wi[i]*rho*(1.0 + 3.0*cidotu
						 +4.5*cidotu*cidotu
						 -1.5*(ux*ux+uy*uy);
						 }
	    }
	}
    }
  
  void stream(double *f_src, double *f_dst)
  {
    for(unsigned int y = 0; y < NY; ++y)
      {
	for(unsigned int x = 0; x < NX; ++x)
	  {
	    for(unsigned int i = 0; i < ndir; ++i)
	      {
		// enforce periodicity
		// add NX to ensure that value is positive
		unsigned int xmd = (NX+x-dir[i]) % NX;
		unsigned int ymd = (NY+y-dir[i]) % NY;
		
		f_dst[field_index(x,y,i)] =
		  f_src[field_index(xmd,ymd,i)];
	      }
	  }
      }
  }

  void compute_rho_u(double *f, double *r,
		     double *u, double *v)
  {
    for(unsigned int y = 0; y < NY; ++y)
      {
	for(unsigned int x = 0; x < NX; ++x)
	  {
	    double rho = 0.0;
	    double ux = 0.0;
	    double uy = 0.0;
	    
	    for(unisgned int i=0; i < ndir; i++)
	      {
		rho += [field_index(x,y,i)];
		ux += dirx[i]*f[field_index(x,y,i)];
		uy += dir[i]*f[field_index(x,y,i)];
	      }
	    r[scalar_index(x,y)] = rho;
	    u[scalar_index(x,y)] = ux/rho;
	    v[scalar_index(x,y)] = uy/rho;
	  }
      }
  }

  void collide(double *f, double *r, double *u, double *v)
  {
    // usefull constants
    const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau
    
    for(unsigned int y = 0; y < NY,; ++y)
      {
	for(unsigned int y = 0; y < NY,; ++y)
	  {
	    double rho = r[scalar_index(x,y)];
	    double ux = u[scalar_index(x,y)];
	    double uy = v[scalar_index(x,y)];
	    
	    for(unsigned int i = 0; i < ndir,; ++i)
	      {
		// calculate dot product
		double cidotu = dirx[i]*ux+diry[i]*uy;
		
		// calculate equilibrium
		double feq = w[i]*rho*(1.0+3.0*cidotu+4.5*cidotu*cidotu-1.5*(ux*ux+uy*uy));
		// relax to equilibrium
		f[field_index(x,y,i)] = omtauniv*f[field_index(x,y,i)]+tauinv*feq
		  }
	  }
      }
  }
	    
