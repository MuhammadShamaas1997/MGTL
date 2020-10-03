#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include "meep.hpp"
#include "ctl-math.h"
#include "ctlgeom.h"
#include "meepgeom.hpp"
#include <math.h>
#ifndef DATADIR
#define DATADIR "./"
#endif

using namespace meep;
using namespace std;


//Argand Diagram, Directional impedance, Poynting flux density color map


typedef std::complex<double> cdouble;

double epsr=10;
//SI Conversion Factors
double a0=1e-2;//0.1mm
double c0=2.99792458e8;//Speed of Light (m/s)
double f0=c0/a0;//300GHz
double t0=1/f0;//0.33e-12 (s)
double mu0=4*pi*(1e-7);// (H/m)
double eps0=8.854187817e-12;// (F/m)
double I0=1; //(A)
double E0=I0/(a0*eps0*c0);//Electric Field
double D0=I0/(a0*c0);//Electric Displacement Field
double B0=I0/(a0*eps0*c0*c0);//Magnetic Field
double H0=I0/(a0);//Magnetizing Field
double sigmaD0=(epsr*eps0*c0)/a0;//Electric Conductivity
double J0=I0/(a0*a0);//Electric Current Density
double u0=(I0*I0)/(eps0*c0*c0*a0*a0);//Energy Density
double S0=(I0*I0)/(eps0*c0*a0*a0);//Poynting Vector
double Sc0=1/(c0);//Courant Factor
double sig0=1e4;
double b0=1e-1;

double xcen=0.0, ycen=0.0, zcen=0.0;
double dxmin=0.0, dxmax=02.5, dymin=0.0, dymax=02.5, dzmin=0.0, dzmax=15.0;
double wcore=2.0*dymax;
double winding_thickness_p=01.0, insulation_thickness_p=02.5,pml_thickness=1.0;  
//double sigma_Cu=100000000;
double mu_core=1.0;
int Np=1;//must be odd
int Ns=1;//must be odd
double margin=02.0;   
double amplitude=1.0;
double divisions=2;

int numcoord=0;
int corecoord=0;
double k=2*pi*Np/(2*dzmin);
double theta=0.0;
double rp=(wcore/2.0)+(winding_thickness_p/2.0)+insulation_thickness_p;
double * xpcoord=new double [1000];
double * ypcoord=new double [1000];
double * zpcoord=new double [1000];
double * xscoord=new double [1000];
double * yscoord=new double [1000];
double * zscoord=new double [1000];
double * xccoord=new double [1000];
double * yccoord=new double [1000];
double * zccoord=new double [1000];


cdouble line_integral_x(fields & f, component C, double dx, double xmin, double xmax, double y, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltax(dx,0.0);
  for (double x=xmin; x<=xmax; x=x+dx)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltax;
  }
  return sum;
}

cdouble line_integral_y(fields & f, component C, double dy, double ymin, double ymax, double x, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltay(dy,0.0);
  for (double y=ymin; y<=ymax; y=y+dy)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltay;
  }
  return sum;
}

cdouble line_integral_z(fields & f, component C, double dz, double zmin, double zmax, double x, double y)
{
  cdouble sum(0.0,0.0);
  cdouble deltaz(dz,0.0);
  for (double z=zmin; z<=zmax; z=z+dz)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
    sum += dF*deltaz;
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
  }
  return sum;
}

cdouble compute_Im(fields & f, double z)
{
  cdouble Iyf=line_integral_y(f,Ey,0.001,ycen-dymax,ycen+dymax,xcen+dxmax,z);
    cdouble Iyb=line_integral_y(f,Ey,0.001,ycen-dymax,ycen+dymax,xcen-dxmax,z);
    cdouble Ixt=line_integral_x(f,Ex,0.001,xcen-dxmax,xcen+dxmax,ycen-dymax,z);
    cdouble Ixb=line_integral_x(f,Ex,0.001,xcen-dxmax,xcen+dxmax,ycen+dymax,z);
    cdouble Im=Ixt+Iyf-Ixb-Iyb;
    return Im;
}

cdouble compute_Ie(fields & f, double z)
{
  cdouble Iyf=line_integral_y(f,Hy,0.001,ycen-dymax,ycen+dymax,xcen+dxmax,z);
    cdouble Iyb=line_integral_y(f,Hy,0.001,ycen-dymax,ycen+dymax,xcen-dxmax,z);
    cdouble Ixt=line_integral_x(f,Hx,0.001,xcen-dxmax,xcen+dxmax,ycen-dymax,z);
    cdouble Ixb=line_integral_x(f,Hx,0.001,xcen-dxmax,xcen+dxmax,ycen+dymax,z);
    cdouble Ie=Ixt+Iyf-Ixb-Iyb;
  return Ie;
}

cdouble compute_Vm(fields & f, double z)
{
  //cdouble Vy=line_integral_y(f,Hy,0.001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Hz,0.001,zcen-dzmax,z,xcen,ycen);
  //cdouble Vm=Vy+Vz;
  return Vz;
}

cdouble compute_Ve(fields & f, double z)
{
  //cdouble Vy=line_integral_y(f,Ey,0.001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Ez,0.001,zcen-dzmax,z,xcen,ycen);
  //cdouble Ve=Vy+Vz;
  return Vz;
}

double mu(const vec &v) 
{
  //return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,mu_core,0.0);
  return 1.0;
}

double eps(const vec &v)
{
  //return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,3.5,0.0);
  return 1.0;
}

double conductivity(const vec &v) 
{ 

}

double sigma(const vec &v) 
{ 
    bool in_middle=false;

  for (int i=0; i<corecoord; i++) {
    double dxp=v.x() - xccoord[i];
    double dyp=v.y() - yccoord[i];
    double dzp=v.z() - zccoord[i];    
    if((dxp<=(wcore/2.0))&&(dyp<=(wcore/2.0))&&(dzp<=(wcore/2.0))){
    in_middle=true;}
    
    }
  double sigp=0.0;
  if(in_middle){sigp=sig0;}
  else{sigp=0.0;}
  return sigp;
}

class core_material : public material_function {
  public:
};

class winding_material : public material_function {
  public:
};


typedef struct my_material_func_data {
  double rxInner, ryInner, rOuter;
  bool with_susceptibility;
} my_material_func_data;

void my_material_func(vector3 p, void *user_data, meep_geom::medium_struct *m) {
  my_material_func_data *data = (my_material_func_data *)user_data;
  
    bool in_middle=false;

    double dx=p.x - xcen;
    double dy=p.y - ycen;
    double dz=p.z - zcen;

  for (int i=0; i<corecoord; i++) {
    double dxp=p.x - xccoord[i];
    double dyp=p.y - yccoord[i];
    double dzp=p.z - zccoord[i];    
    //double drp=sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
    if((dxp<=(wcore/2.0))&&(dyp<=(wcore/2.0))&&(dzp<=(wcore/2.0))){
    //if(drp<=(wcore/2.0))
      in_middle=true;}
    
    }

  // set permittivity and permeability
  double nn = in_middle ? sqrt(mu_core) : 1.0;
  double mm = in_middle ? sqrt(1.0) : 1.0;
  m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = 10.0;
  m->mu_diag.x = m->mu_diag.y = m->mu_diag.z = nn*nn;
  m->E_chi2_diag.x = m->E_chi2_diag.y = m->E_chi2_diag.z = 0.0;
  m->E_chi3_diag.x = m->E_chi3_diag.y = m->E_chi3_diag.z = 0.0;
  m->H_chi2_diag.x = m->H_chi2_diag.y = m->H_chi2_diag.z = 0.0;
  m->H_chi3_diag.x = m->H_chi3_diag.y = m->H_chi3_diag.z = 0.0;

  if (in_middle)
  {
/*    m->H_susceptibilities.num_items = 1;
    m->H_susceptibilities.items = new meep_geom::susceptibility[1];

    m->H_susceptibilities.items[0].sigma_offdiag.x = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.y = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.z = 0.0;
    m->H_susceptibilities.items[0].sigma_diag.x = sig0;
    m->H_susceptibilities.items[0].sigma_diag.y = sig0;
    m->H_susceptibilities.items[0].sigma_diag.z = sig0;
    m->H_susceptibilities.items[0].bias.x = 0.0;
    m->H_susceptibilities.items[0].bias.y = 0.0;
    m->H_susceptibilities.items[0].bias.z = b0;
    m->H_susceptibilities.items[0].frequency = (0.6667e-5);
    m->H_susceptibilities.items[0].gamma = (1e-6);;
    m->H_susceptibilities.items[0].alpha = 0;
    m->H_susceptibilities.items[0].noise_amp = 0.0;
    m->H_susceptibilities.items[0].drude = true;
    m->H_susceptibilities.items[0].saturated_gyrotropy = false;
    m->H_susceptibilities.items[0].is_file = false;
    */m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = (5e-3)/sigmaD0;
  }

  for (int i=0; i<numcoord; i++) {
  double dxp=p.x - xpcoord[i];
    double dyp=p.y - ypcoord[i];
    double dzp=p.z - zpcoord[i];    
    double drp=sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
    if(drp<=(winding_thickness_p/2.0)){
      m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = (3.5e7)/sigmaD0;
    }
    double dxs=p.x - xscoord[i];
    double dys=p.y - yscoord[i];
    double dzs=p.z - zscoord[i];
    double drs=sqrt(dxs*dxs+dys*dys+dzs*dzs);
    if(drs<=(winding_thickness_p/2.0)){
      m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = (3.5e7)/sigmaD0;
    }
  }

}


class anisodisp_materialH : public material_function {
public:
  virtual void sigma_row(component c, double sigrow[3], const vec &r) {
    (void)r; // unused
    if (component_direction(c) == X) {
      sigrow[0] = sig0;
      sigrow[1] = sig0;
      sigrow[2] = sig0;
    }
    else if (component_direction(c) == Y) {
      sigrow[0] = sig0;
      sigrow[1] = sig0;
      sigrow[2] = sig0;
    }
    else {
      sigrow[0] = sig0;
      sigrow[1] = sig0;
      sigrow[2] = sig0;
    }
  }
};




int main(int argc, char *argv[]) {
  
  initialize mpi(argc,argv);
  const char *mydirname = "MMTL-out";
    std::ofstream Time;
    std::ofstream Space;
    std::ofstream FieldsIn;
    std::ofstream FieldsOut;
    std::ofstream Fluxes;
    std::ofstream Skin;
    std::ofstream Permeability;
    std::ofstream Permittivity;
    std::ofstream Chi1inv;
    std::ofstream SourceFFT;

    Time.open ("TimeEvolution.txt");
    Space.open ("SpaceEvolution.txt");
    FieldsIn.open ("FieldEvolutionIn.txt");
    FieldsOut.open ("FieldEvolutionOut.txt");
    Fluxes.open ("Flux.txt");
    Skin.open ("Skin.txt");
    Permeability.open ("Permeability.txt");
    Permittivity.open ("Permittivity.txt");
    SourceFFT.open ("SourceFFT.txt");

    //trash_output_directory(mydirname);
    double xsize=50;
    double ysize=50;
    double zsize=50;
    //double xsize=6, ysize=6, zsize=6;
    
    k=2*pi*Np/(03.0);
    double z=-9.0;
    for (theta = 0.0; theta <= (2*pi); theta=theta+0.1)
    {
      xpcoord[numcoord]=-rp*sin(theta);
      ypcoord[numcoord]=rp*cos(theta);
      zpcoord[numcoord]=z;
      xscoord[numcoord]=-rp*sin(theta);
      yscoord[numcoord]=rp*cos(theta);;
      zscoord[numcoord]=z+18.0;
      numcoord++;
    }


  for(double y=-dzmax;y<=dzmax;y=y+dzmax)
  {for (double z=-dzmax;z<=dzmax;z=z+00.1){
    xccoord[corecoord]=xcen;
    yccoord[corecoord]=y;
    zccoord[corecoord]=z;
    corecoord++;
  }}

  for(double z=-dzmax;z<=dzmax;z=z+(2*dzmax))
  {for (double y=-dzmax-(wcore/2);y<=dzmax+(wcore/2);y=y+00.1){
    xccoord[corecoord]=xcen;
    yccoord[corecoord]=y;
    zccoord[corecoord]=z;
    corecoord++;
  }}

  /*for(double z=-07.5;z<=07.5;z=z+15.0)
  {for (double y=-12.5;y<=12.5;y=y+00.1){
    xccoord[corecoord]=xcen;
    yccoord[corecoord]=y;
    zccoord[corecoord]=z;
    corecoord++;
  }}*/


    grid_volume gv = vol3d(xsize, ysize, zsize, divisions);
    gv.center_origin();
  

    structure transformer = new structure(gv, eps,pml(1.0));
    transformer.set_output_directory(mydirname);
    
    meep_geom::medium_struct my_medium_struct;
    my_medium_struct.epsilon_diag.x = 1.0;
    my_medium_struct.epsilon_diag.y = 1.0;
    my_medium_struct.epsilon_diag.z = 1.0;
    my_medium_struct.mu_diag.x=1.0;
    my_medium_struct.mu_diag.y=1.0;
    my_medium_struct.mu_diag.z=1.0;

    my_medium_struct.H_chi2_diag.x=0.0;
    my_medium_struct.H_chi2_diag.y=0.0;
    my_medium_struct.H_chi2_diag.z=0.0;
    my_medium_struct.H_chi3_diag.x=0.0;
    my_medium_struct.H_chi3_diag.y=0.0;
    my_medium_struct.H_chi3_diag.z=0.0;

    my_medium_struct.E_chi2_diag.x=0.0;
    my_medium_struct.E_chi2_diag.y=0.0;
    my_medium_struct.E_chi2_diag.z=0.0;
    my_medium_struct.E_chi3_diag.x=0.0;
    my_medium_struct.E_chi3_diag.y=0.0;
    my_medium_struct.E_chi3_diag.z=0.0;

    /*my_medium_struct.H_susceptibilities.num_items = 1;
    my_medium_struct.H_susceptibilities.items = new meep_geom::susceptibility[1];

    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.x = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.y = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.z = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.x = sig0;//(1e13)*f0*f0;
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.y = sig0;//(1e13)*f0*f0;
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.z = sig0;//(1e13)*f0*f0;
    my_medium_struct.H_susceptibilities.items[0].bias.x = 0.0;
    my_medium_struct.H_susceptibilities.items[0].bias.y = 0.0;
    my_medium_struct.H_susceptibilities.items[0].bias.z = b0;
    my_medium_struct.H_susceptibilities.items[0].frequency = 5.8e-5;
    my_medium_struct.H_susceptibilities.items[0].gamma = -0.33e-2;//*f0;
    my_medium_struct.H_susceptibilities.items[0].alpha = 0;
    my_medium_struct.H_susceptibilities.items[0].noise_amp = 0.0;
    my_medium_struct.H_susceptibilities.items[0].drude = true;
    my_medium_struct.H_susceptibilities.items[0].saturated_gyrotropy = false;
    my_medium_struct.H_susceptibilities.items[0].is_file = false;
    */my_material_func_data data;
    data.with_susceptibility = true;


  meep_geom::material_type default_material =meep_geom::make_dielectric(1.0);     
  default_material->medium=my_medium_struct;
  default_material->user_func = my_material_func;
  default_material->user_data = (void *) &data;
  default_material->do_averaging = false;  

    meep_geom::material_type my_user_material =meep_geom::make_user_material(my_material_func, (void *)&data, false);
    
    geometric_object objects[5];
  vector3 center = {0.0, 0.0, 0.0};  
  vector3 center1 = {0.0, -dzmax, 0.0};
  vector3 center2 = {0.0, 0.0, 0.0};
  vector3 center3 = {0.0, dzmax, 0.0};
  vector3 center4 = {0.0, 0.0, dzmax};
  vector3 center5 = {0.0, 0.0, -dzmax};
  double radius = 3.0;
  double height = 1.0e20;
  vector3 xhat1 = {1.0, 0.0, 0.0};
  vector3 yhat1 = {0.0, 1.0, 0.0};
  vector3 zhat1 = {0.0, 0.0, 1.0};
  vector3 size1 = {wcore, wcore, 2*dzmax};
  vector3 size2 = {wcore, wcore, 2*dzmax};
  vector3 size3 = {wcore, wcore, 2*dzmax};
  vector3 size4 = {wcore, 2*dzmax+(wcore), wcore};
  vector3 size5 = {wcore, 2*dzmax+(wcore), wcore};
  objects[0] = make_block(my_user_material, center1, xhat1, yhat1, zhat1, size1);
  objects[1] = make_block(my_user_material, center2, xhat1, yhat1, zhat1, size2);
  objects[2] = make_block(my_user_material, center3, xhat1, yhat1, zhat1, size3);
  objects[3] = make_block(my_user_material, center4, xhat1, yhat1, zhat1, size4);
  objects[4] = make_block(my_user_material, center5, xhat1, yhat1, zhat1, size5);
  geometric_object_list g = {5, objects};
  


    meep_geom::absorber_list al= new meep_geom::absorber_list_type;
    meep_geom::material_type_list mtl= meep_geom::material_type_list();    
    mtl.num_items=1;
    mtl.items=new meep_geom::material_type;
    mtl.items[0]=my_user_material;
    cout<<"MAIN "<<mtl.num_items<<endl;
    bool use_anisotropic_averaging = false;
    bool ensure_periodicity = false;
    set_materials_from_geometry(&transformer, g, center, use_anisotropic_averaging,DEFAULT_SUBPIXEL_TOL, DEFAULT_SUBPIXEL_MAXEVAL,ensure_periodicity, default_material,al,mtl);

    anisodisp_materialH anisodispmatH;
    transformer.add_susceptibility(sigma, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,b0),(5.8e-5), (-0.33e-2)/(2*3.1459),0.0,GYROTROPIC_DRUDE));
    //transformer.add_susceptibility(sigma, H_stuff, lorentzian_susceptibility(5.8e-5, (-0.33e-2)/(2*3.1459)));

    fields f(& transformer);
  

for (double fp=0.0;fp<=(1e9)/f0;fp=fp+(1e4)/f0)
    { 
      double yp = 0.0;
      {
        Permeability<<(fp*f0)<<" "<<f.get_mu(vec(0.0,0.0,0.0),fp)<<endl;
        Permittivity<<(fp*f0)<<" "<<f.get_eps(vec(0.0,0.0,0.0),fp)<<endl;
      }
      
    }
    

    //f_range;1e-5,1e-2
    double fcen = (0.5e9)/f0; // ; pulse center frequency
    double df = 0.999999*((1.0e9)/f0);    // ; df
    //continuous_src_time src(cdouble(fcen,df));
    gaussian_src_time src(fcen,df);

for (double fp=0.0;fp<=(1e9)/f0;fp=fp+(1e6)/f0)
    { 
      double yp = 0.0;
      {
        SourceFFT<<fp<<" "<<src.fourier_transform(fp).real()<<" "<<src.fourier_transform(fp).imag()<<endl;
      }
      
    }


for (int i=0;i<numcoord;i++)
    {
      theta=atan(-xpcoord[i]/ypcoord[i]);
      const volume vsrc1 =volume(vec(xpcoord[i],ypcoord[i],zpcoord[i]), vec(xpcoord[i],ypcoord[i],zpcoord[i]));
      {
        f.add_volume_source(Ex, src, vsrc1, cdouble(amplitude*cos(theta),0));
        f.add_volume_source(Ey, src, vsrc1, cdouble(-amplitude*sin(theta),0));
      }
    
    }
cout<<"DFT"<<endl;
    /*volume box1( vec(rp,-rp,-10.0), vec(-rp,rp,-8.0) );
    volume box2( vec(rp,-rp,8.0), vec(-rp,rp,10.0) );

    fcen = (5e9)/f0; // ; pulse center frequency
    df = 0.9999999*(fcen/f0);    // ; df
  
    double fmin = (1)/f0, fmax = (1e9)/f0;
    int Nfreq = 1000;
    dft_flux flux1 = f.add_dft_flux_box(box1, fmin, fmax, Nfreq);
    dft_flux flux2 = f.add_dft_flux_box(box2, fmin, fmax, Nfreq);
    double init_energy = f.field_energy_in_box(box1);
*/
cout<<"Stepping"<<endl;
  f.step();
  f.step();
  f.step();
  f.step();
  f.step();
  cout<<"Okay"<<endl;
    volume vxy=volume(vec(-xsize,-ysize,0),vec(xsize,ysize,0));
    volume vxz=volume(vec(-xsize,0,-zsize),vec(xsize,0,zsize));
    volume vyz=volume(vec(0,-ysize,-zsize),vec(0,ysize,zsize));
    
int stop=0;

    for(int i=1;i<=1000000;i++)
    {
      if(!stop)
         {f.step();}
       if (stop)
       {
         i=1000001;
       }
       

    if ((i%100)==0)
    {
    /*f.output_hdf5(Hx,vyz);
    f.output_hdf5(Hy,vyz);
    f.output_hdf5(Hz,vyz);
    f.output_hdf5(Bx,vyz);
    f.output_hdf5(By,vyz);
    f.output_hdf5(Bz,vyz);
    f.output_hdf5(Ex,vyz);
    f.output_hdf5(Ey,vyz);
    f.output_hdf5(Ez,vyz);
    f.output_hdf5(Dx,vyz);
    f.output_hdf5(Dy,vyz);
    f.output_hdf5(Dz,vyz);
    f.output_hdf5(Sx,vyz);
    f.output_hdf5(Sy,vyz);
    f.output_hdf5(Sz,vyz);
    
    double *fl1 = flux1.flux();
    double *fl2 = flux2.flux();
    cout<<"Flux Harmonics"<<endl;
    for (int i = 0; i < Nfreq; ++i) {
      Fluxes<<(fmin + i * flux1.dfreq)<<" , "<<fl1[i]<<" , "<<fl2[i]<<endl;
      //freq , fluxin , fluxout
    }*/ 

    }    
       if (i<=600)
      {
        //cdouble Vm=compute_Vm(f,zcen);
        //cdouble Im=compute_Im(f,zcen);
        //cdouble Ve=compute_Ve(f,zcen);
        //cdouble Ie=compute_Ie(f,zcen);
        //Time<<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<endl;//" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
        //Im , Vm , Ie , Ve 

        for (double z=-dzmax;z<=dzmax;z=z+1/(2*divisions))  
      {
        monitor_point pin;
        f.get_point(&pin, vec(xcen,ycen,z));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);
        
        FieldsIn<<i<<" , "<<z<<" , "<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;
      }
    
      }

      if((i==(600))) 
      {
        //cout<<"End? (1/0):";
        //cin>>stop;
        stop=1;
      }
       
    }

    
    /*double *fl1 = flux1.flux();
    double *fl2 = flux2.flux();
    cout<<"Flux Harmonics"<<endl;
    for (int i = 0; i < Nfreq; ++i) {
      Fluxes<<(fmin + i * flux1.dfreq)<<" , "<<fl1[i]<<" , "<<fl2[i]<<endl;
      //freq , fluxin , fluxout
    }*/ 
    
    cout<<"Skin Effect"<<endl;
    for (double y=-05.0;y<=05.0;y=y+00.001)  
    {
      monitor_point pin;
        f.get_point(&pin, vec(xcen,y,zcen));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);
        Skin<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;
    }

    cout<<"SpaceEvolution"<<endl;
    for (double z=-dzmax;z<=dzmax;z=z+1.0)  
    {
      cdouble Im=compute_Im(f,z);
      cdouble Ie=compute_Ie(f,z);
      cdouble Vm=compute_Vm(f,z);
      cdouble Ve=compute_Ve(f,z);
      Space <<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
      //Im , Vm , Ie , Ve
    }

/*    f.step();
    f.output_hdf5(Hx,f.gv);
*/

    Time.close();
    Space.close();
    FieldsIn.close();
    FieldsOut.close();
    Fluxes.close();
    return 0;
}

