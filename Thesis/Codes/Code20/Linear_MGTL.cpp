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





typedef std::complex<double> cdouble;

//SI Conversion Factors
double d0=1e-3;//1mm
double c0=3e8;//Speed of Light (m/s)
double f0=c0/d0;//300GHz
double t0=1/f0;//0.33e-11 (s)
double mu0=4*pi*(1e-7);// (H/m)
double eps0=8.85e-7;// (F/m)



double xcen=0.0, ycen=0.0, zcen=0.0;
double dxmin=0.0, dxmax=0.5, dymin=0.0, dymax=2.5, dzmin=0.0, dzmax=0.5;
double wcore=dymax - dymin;
double winding_thickness_p=0.2, insulation_thickness_p=0.1,pml_thickness=1.0;  
//double sigma_Cu=100000000;
double mu_core=1000;
int Np=1;//must be odd
int Ns=1;//must be odd
double margin=0.1;   
double amplitude=100.0;
double divisions=10;

//Copper
//metal_range = mp.FreqRange(min=um_scale/12.398, max=um_scale/.20664)
double um_scale = 1000.0;//1000um
double eV_um_scale = 1.0/1.23984193;
double Cu_plasma_frq = 10.83*eV_um_scale;
double Cu_f0 = 0.575;
double Cu_frq0 = 1e-10;
double Cu_gam0 = 0.030*eV_um_scale;
double Cu_sig0 = Cu_f0*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq0*Cu_frq0);
double Cu_f1 = 0.061;
double Cu_frq1 = 0.291*eV_um_scale;
double Cu_gam1 = 0.378*eV_um_scale;
double Cu_sig1 = Cu_f1*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq1*Cu_frq1);
double Cu_f2 = 0.104;
double Cu_frq2 = 2.957*eV_um_scale;      
double Cu_gam2 = 1.056*eV_um_scale;
double Cu_sig2 = Cu_f2*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq2*Cu_frq2);
double Cu_f3 = 0.723;
double Cu_frq3 = 5.300*eV_um_scale;      
double Cu_gam3 = 3.213*eV_um_scale;
double Cu_sig3 = Cu_f3*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq3*Cu_frq3);
double Cu_f4 = 0.638;
double Cu_frq4 = 11.18*eV_um_scale;
double Cu_gam4 = 4.305*eV_um_scale;
double Cu_sig4 = Cu_f4*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq4*Cu_frq4);
double sigma_Cu=Cu_sig0;


//NiFe
//NiFe_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.25)
double NiFe_frq = 1/(0.0838297450980392*um_scale);
double NiFe_gam = 1/(0.259381156903766*um_scale);
double NiFe_sig = 1;



double annulus(const vec & v, double xcen, double ycen, double zcen, double dxmin, double dxmax, double dymin, double dymax, double dzmin, double dzmax, double control, double special, double std = 0.0)
{
  double mu=control;

  double dx=v.x() - xcen;
  double dy=v.y() - ycen;
  double dz=v.z() - zcen;

  if ( (abs(dx)<=dxmax) && (abs(dy)<=dymax) && (abs(dz)<=dzmax) )
  {
    {mu= gaussian_random(special,std);}
  }

  if ( (abs(dx)<=dxmin) && (abs(dy)<=dymin) && (abs(dz)<=dzmin) )
  {
    mu= control;
  }

    return mu;  
}

bool inside_box(const vec & v, double xcen, double ycen, double zcen, double dxmax, double dymax, double dzmax)
{
  bool ans=false;

  double dx=v.x() - xcen;
  double dy=v.y() - ycen;
  double dz=v.z() - zcen;


  if ( (abs(dx)<=dxmax) && (abs(dy)<=dymax) && (abs(dz)<=dzmax) )
  {
    {ans=true;}
  }

  return ans;
}

  
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

cdouble compute_Im(fields & f, double y)
{
  cdouble Izf=line_integral_z(f,Ez,0.0001,zcen+dzmin,zcen+dzmax,xcen+dxmax,y);
    cdouble Izb=line_integral_z(f,Ez,0.0001,zcen+dzmin,zcen+dzmax,xcen-dxmax,y);
    cdouble Ixt=line_integral_x(f,Ex,0.0001,xcen-dxmax,xcen+dxmax,y,zcen+dzmax);
    cdouble Ixb=line_integral_x(f,Ex,0.0001,xcen-dxmax,xcen+dxmax,y,zcen+dzmin);
    cdouble Im=Ixt-Izf-Ixb+Izb;
    return Im;
}

cdouble compute_Ie(fields & f, double y)
{
  cdouble Izf=line_integral_z(f,Hz,0.0001,zcen+dzmin,zcen+dzmax,xcen+dxmax,y);
  cdouble Izb=line_integral_z(f,Hz,0.0001,zcen+dzmin,zcen+dzmax,xcen-dxmax,y);
  cdouble Ixt=line_integral_x(f,Hx,0.0001,xcen-dxmax,xcen+dxmax,y,zcen+dzmax);
  cdouble Ixb=line_integral_x(f,Hx,0.0001,xcen-dxmax,xcen+dxmax,y,zcen+dzmin);
  cdouble Ie=Ixt-Izf-Ixb+Izb;
  return Ie;
}

cdouble compute_Vm(fields & f, double y)
{
  cdouble Vy=line_integral_y(f,Hy,0.0001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Hz,0.0001,zcen-dzmin,zcen+dzmin,xcen,y);
  cdouble Vm=Vy+Vz;
  return Vm;
}

cdouble compute_Ve(fields & f, double y)
{
  cdouble Vy=line_integral_y(f,Ey,0.0001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Ez,0.0001,zcen-dzmin,zcen+dzmin,xcen,y);
  cdouble Ve=Vy+Vz;
  return Ve;
}

double mu(const vec &v) 
{
  return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,mu_core,0.0);
}

double eps(const vec &v){
  double dxmax1=dxmax+2*insulation_thickness_p+winding_thickness_p;
  double dxmin1=dxmax;
  double dymax1=dymax+2*insulation_thickness_p+winding_thickness_p;
  double dymin1=dymin-2*insulation_thickness_p-winding_thickness_p;
  double dzmax1=dzmax+2*insulation_thickness_p+winding_thickness_p;
  double dzmin1=dzmin-2*insulation_thickness_p-winding_thickness_p;

  return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,3.5,0.0);
  //return 3.0;
}

double conductivity(const vec &v) 
{ 
   double sig=1.0;
   double xcenp=xcen; 
    double ycenp=ycen-dymin-(0.5*wcore); 
    double zcenp=zcen;

    double dxminp=dxmax+insulation_thickness_p;
    double dxmaxp=dxminp+winding_thickness_p;
    double dyminp=(0.5*wcore)+insulation_thickness_p; 
    double dymaxp=dyminp+winding_thickness_p;
    double dzminp=0; 
    double dzmaxp=0.5*winding_thickness_p;    
  
    zcenp=zcen-(0.5*double(Np-1)*insulation_thickness_p)-(0.5*double(Np-1)*winding_thickness_p);

    for (int i=0;i<Np;i++)
    {

      if (inside_box(v,xcenp,ycenp,zcenp,dxmaxp,dymaxp,dzmaxp))
    {
      return annulus(v,xcenp,ycenp,zcenp,dxminp,dxmaxp,dyminp,dymaxp,dzminp,dzmaxp,sig,sigma_Cu,10);
    }

    zcenp=zcenp+insulation_thickness_p+winding_thickness_p;     
    }


    ycenp=ycen+dymin+(0.5*wcore); 
    zcenp=zcen-(0.5*double(Ns-1)*insulation_thickness_p)-(0.5*double(Ns-1)*winding_thickness_p);

    for (int i=0;i<Ns;i++)
    {

      if (inside_box(v,xcenp,ycenp,zcenp,dxmaxp,dymaxp,dzmaxp))
    {
      return annulus(v,xcenp,ycenp,zcenp,dxminp,dxmaxp,dyminp,dymaxp,dzminp,dzmaxp,sig,sigma_Cu,10);
    }

    zcenp=zcenp+insulation_thickness_p+winding_thickness_p;     
    }

    return sig;
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

    if ( (abs(dx)<=dxmax) && (abs(dy)<=dymax) && (abs(dz)<=dzmax) )
    {
    in_middle=true;
    }

    if ( (abs(dx)<=dxmin) && (abs(dy)<=dymin) && (abs(dz)<=dzmin) )
    {
    in_middle=false;
    }

  // set permittivity and permeability
  double nn = in_middle ? sqrt(mu_core) : 1.0;
  double mm = in_middle ? sqrt(3.5) : 1.0;
  m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = mm * mm;
  //m->epsilon_offdiag.x.re = m->epsilon_offdiag.x.im = epsilon_offdiag.y.re = m->epsilon_offdiag.y.im = epsilon_offdiag.z.re = m->epsilon_offdiag.z.im = nn * nn;
  m->mu_diag.x = m->mu_diag.y = m->mu_diag.z = nn*nn;

  //m->mu_offdiag.x.re = m->mu_offdiag.x.im = mu_offdiag.y.re = m->mu_offdiag.y.im = mu_offdiag.z.re = m->mu_offdiag.z.im = nn * nn;
  m->E_chi2_diag.x = m->E_chi2_diag.y = m->E_chi2_diag.z = 0.0;
  m->E_chi3_diag.x = m->E_chi3_diag.y = m->E_chi3_diag.z = 0.0;
  m->H_chi2_diag.x = m->H_chi2_diag.y = m->H_chi2_diag.z = nn * nn;
  m->H_chi3_diag.x = m->H_chi3_diag.y = m->H_chi3_diag.z = nn * nn;
  //m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = nn * nn;
  //m->B_conductivity_diag.x = m->B_conductivity_diag.y = m->B_conductivity_diag.z = 0.0;

  // add susceptibilities
  {

    double xcenp=xcen; 
    double ycenp=ycen-dymax; 
    double zcenp=zcen;
    double dxminp=dxmax+insulation_thickness_p;
    double dxmaxp=dxminp+winding_thickness_p;
    double dyminp=0; 
    double dymaxp=0.5*winding_thickness_p;
    double dzminp=(0.5*wcore)+insulation_thickness_p; 
    double dzmaxp=dzminp+winding_thickness_p;    
    ycenp=ycen-(0.5*double(Np-1)*insulation_thickness_p)-(0.5*double(Np-1)*winding_thickness_p);


    for (int i=0;i<Np;i++)
    {

    double dx=p.x - xcenp;
    double dy=p.y - ycenp;
    double dz=p.z - zcenp;


    if ( (abs(dx)<=dxmaxp) && (abs(dy)<=dymaxp) && (abs(dz)<=dzmaxp) )
    {
    m->E_susceptibilities.num_items = 5;
    m->E_susceptibilities.items = new meep_geom::susceptibility[5];

    m->E_susceptibilities.items[0].sigma_diag.x = Cu_sig0;
    m->E_susceptibilities.items[0].sigma_diag.y = Cu_sig0;
    m->E_susceptibilities.items[0].sigma_diag.z = Cu_sig0;
    m->E_susceptibilities.items[0].frequency = Cu_frq0;
    m->E_susceptibilities.items[0].gamma = Cu_gam0;
    m->E_susceptibilities.items[0].drude = true;

    m->E_susceptibilities.items[1].sigma_diag.x = Cu_sig1;
    m->E_susceptibilities.items[1].sigma_diag.y = Cu_sig1;
    m->E_susceptibilities.items[1].sigma_diag.z = Cu_sig1;
    m->E_susceptibilities.items[1].frequency = Cu_frq1;
    m->E_susceptibilities.items[1].gamma = Cu_gam1;
    m->E_susceptibilities.items[1].drude = true;

    m->E_susceptibilities.items[2].sigma_diag.x = Cu_sig2;
    m->E_susceptibilities.items[2].sigma_diag.y = Cu_sig2;
    m->E_susceptibilities.items[2].sigma_diag.z = Cu_sig2;
    m->E_susceptibilities.items[2].frequency = Cu_frq2;
    m->E_susceptibilities.items[2].gamma = Cu_gam2;
    m->E_susceptibilities.items[2].drude = true;

    m->E_susceptibilities.items[3].sigma_diag.x = Cu_sig3;
    m->E_susceptibilities.items[3].sigma_diag.y = Cu_sig3;
    m->E_susceptibilities.items[3].sigma_diag.z = Cu_sig3;
    m->E_susceptibilities.items[3].frequency = Cu_frq3;
    m->E_susceptibilities.items[3].gamma = Cu_gam3;
    m->E_susceptibilities.items[3].drude = true;

    m->E_susceptibilities.items[4].sigma_diag.x = Cu_sig4;
    m->E_susceptibilities.items[4].sigma_diag.y = Cu_sig4;
    m->E_susceptibilities.items[4].sigma_diag.z = Cu_sig4;
    m->E_susceptibilities.items[4].frequency = Cu_frq4;
    m->E_susceptibilities.items[4].gamma = Cu_gam4;
    m->E_susceptibilities.items[4].drude = true;

    }

    ycenp=ycenp+insulation_thickness_p+winding_thickness_p;     
    }

    /*
    ycenp=ycen+dymin+(0.5*wcore); 
    zcenp=zcen-(0.5*double(Ns-1)*insulation_thickness_p)-(0.5*double(Ns-1)*winding_thickness_p);

    for (int i=0;i<Ns;i++)
    {

    double dx=p.x - xcenp;
    double dy=p.y - ycenp;
    double dz=p.z - zcenp;


    if ( (abs(dx)<=dxmaxp) && (abs(dy)<=dymaxp) && (abs(dz)<=dzmaxp) )
    {
    m->E_susceptibilities.num_items = 5;
    m->E_susceptibilities.items = new meep_geom::susceptibility[5];

    m->E_susceptibilities.items[0].sigma_diag.x = Cu_sig0;
    m->E_susceptibilities.items[0].sigma_diag.y = Cu_sig0;
    m->E_susceptibilities.items[0].sigma_diag.z = Cu_sig0;
    m->E_susceptibilities.items[0].frequency = Cu_frq0;
    m->E_susceptibilities.items[0].gamma = Cu_gam0;
    m->E_susceptibilities.items[0].drude = true;

    m->E_susceptibilities.items[1].sigma_diag.x = Cu_sig1;
    m->E_susceptibilities.items[1].sigma_diag.y = Cu_sig1;
    m->E_susceptibilities.items[1].sigma_diag.z = Cu_sig1;
    m->E_susceptibilities.items[1].frequency = Cu_frq1;
    m->E_susceptibilities.items[1].gamma = Cu_gam1;
    m->E_susceptibilities.items[1].drude = true;

    m->E_susceptibilities.items[2].sigma_diag.x = Cu_sig2;
    m->E_susceptibilities.items[2].sigma_diag.y = Cu_sig2;
    m->E_susceptibilities.items[2].sigma_diag.z = Cu_sig2;
    m->E_susceptibilities.items[2].frequency = Cu_frq2;
    m->E_susceptibilities.items[2].gamma = Cu_gam2;
    m->E_susceptibilities.items[2].drude = true;

    m->E_susceptibilities.items[3].sigma_diag.x = Cu_sig3;
    m->E_susceptibilities.items[3].sigma_diag.y = Cu_sig3;
    m->E_susceptibilities.items[3].sigma_diag.z = Cu_sig3;
    m->E_susceptibilities.items[3].frequency = Cu_frq3;
    m->E_susceptibilities.items[3].gamma = Cu_gam3;
    m->E_susceptibilities.items[3].drude = true;

    m->E_susceptibilities.items[4].sigma_diag.x = Cu_sig4;
    m->E_susceptibilities.items[4].sigma_diag.y = Cu_sig4;
    m->E_susceptibilities.items[4].sigma_diag.z = Cu_sig4;
    m->E_susceptibilities.items[4].frequency = Cu_frq4;
    m->E_susceptibilities.items[4].gamma = Cu_gam4;
    m->E_susceptibilities.items[4].drude = true;
    }

    zcenp=zcenp+insulation_thickness_p+winding_thickness_p;     
    }
	*/

    if (in_middle)
    {
    m->H_susceptibilities.num_items = 1;
    m->H_susceptibilities.items = new meep_geom::susceptibility[1];

    m->H_susceptibilities.items[0].sigma_offdiag.x = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.y = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.z = 0.0;
    m->H_susceptibilities.items[0].sigma_diag.x = NiFe_sig;
    m->H_susceptibilities.items[0].sigma_diag.y = NiFe_sig;
    m->H_susceptibilities.items[0].sigma_diag.z = NiFe_sig;
    m->H_susceptibilities.items[0].bias.x = 0.0;
    m->H_susceptibilities.items[0].bias.y = 0.0;
    m->H_susceptibilities.items[0].bias.z = 0.0;
    m->H_susceptibilities.items[0].frequency = 0.005;
    m->H_susceptibilities.items[0].gamma = NiFe_gam;
    m->H_susceptibilities.items[0].alpha = NiFe_gam;
    m->H_susceptibilities.items[0].noise_amp = 0.0;
    m->H_susceptibilities.items[0].drude = false;
    m->H_susceptibilities.items[0].saturated_gyrotropy = true;
    m->H_susceptibilities.items[0].is_file = false;
    }
  }
}


int main(int argc, char *argv[]) {
  
  initialize mpi(argc,argv);
  const char *mydirname = "MMTL-out";
    std::ofstream Time;
    std::ofstream Space;
    std::ofstream FieldsIn;
    std::ofstream FieldsOut;
    std::ofstream Fluxes;
    Time.open ("TimeEvolution.txt");
    Space.open ("SpaceEvolution.txt");
    FieldsIn.open ("FieldEvolutionIn.txt");
    FieldsOut.open ("FieldEvolutionOut.txt");
    Fluxes.open ("Flux.txt");
    //trash_output_directory(mydirname);
    //double xsize=2*(dxmax+2*insulation_thickness_p+winding_thickness_p+pml_thickness);
    //double ysize=2*(dymax+2*insulation_thickness_p+winding_thickness_p+pml_thickness);
    //double zsize=2*(dzmax+2*insulation_thickness_p+winding_thickness_p+pml_thickness);
    double xsize=6, ysize=10, zsize=6;
  
    grid_volume gv = vol3d(xsize, ysize, zsize, divisions);
    //grid_volume vol3d(double xsize, double ysize, double zsize, double a);
    gv.center_origin();
    //void center_origin(void) { shift_origin(-icenter()); }
  

    structure transformer(gv, eps, pml(pml_thickness));
    //transformer.set_epsilon(eps,true);
    //transformer.set_mu(mu,true);
    //transformer.add_susceptibility(core_mat, H_stuff, lorentzian_susceptibility(2*pi*NiFe_frq, NiFe_gam, true));
    //transformer.add_susceptibility(core_mat, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,0.0),2*pi*NiFe_frq, NiFe_gam, true));
    //gyrotropic_susceptibility(const vec &bias, double omega_0, double gamma, double alpha = 0.0,gyrotropy_model model = GYROTROPIC_LORENTZIAN);
    transformer.set_output_directory(mydirname);
    //transformer->set_chi2(mu);
    //transformer->set_chi3(mu);

    my_material_func_data data;
    data.with_susceptibility = true;
    meep_geom::material_type my_material =
    meep_geom::make_user_material(my_material_func, (void *)&data, false);
  
    geometric_object_list g = {0,0};

    vector3 center = {0, 0, 0};
    bool use_anisotropic_averaging = true;
    bool ensure_periodicity = true;
    set_materials_from_geometry(&transformer, g, center, use_anisotropic_averaging,DEFAULT_SUBPIXEL_TOL, DEFAULT_SUBPIXEL_MAXEVAL,ensure_periodicity, my_material);

    /*    
    winding_material winding_mat;
    structure winding(gv, winding_mat);
    //winding.set_epsilon(eps);
    winding.set_conductivity(Ex,conductivity);
    winding.set_conductivity(Ey,conductivity);
    winding.set_conductivity(Ez,conductivity);    
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq0, Cu_gam0, true));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq1, Cu_gam1));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq2, Cu_gam2));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq3, Cu_gam3));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq4, Cu_gam4));
    //lorentzian_susceptibility(double omega_0, double gamma, bool no_omega_0_denominator = false): omega_0(omega_0), gamma(gamma), no_omega_0_denominator(no_omega_0_denominator)
    */      


    fields f(& transformer);
    //fields(structure *, double m = 0, double beta = 0, bool zero_fields_near_cylorigin = true);
  
 
    //f_range;1e-5,1e-2
    double fcen = 0.005; // ; pulse center frequency
    double df = 0.00499;    // ; df
    continuous_src_time src(cdouble(fcen,0));
    //gaussian_src_time src(fcen,df);

    double xcenp=xcen; 
    double ycenp=ycen-dymax; 
    double zcenp=zcen;
    double dxminp=dxmax+insulation_thickness_p;
    double dxmaxp=dxminp+winding_thickness_p;
    double dyminp=0; 
    double dymaxp=0.5*winding_thickness_p;
    double dzminp=(0.5*wcore)+insulation_thickness_p; 
    double dzmaxp=dzminp+winding_thickness_p;    
    ycenp=ycen-(0.5*double(Np-1)*insulation_thickness_p)-(0.5*double(Np-1)*winding_thickness_p);

    for (int i=0;i<Np;i++)
    {
      const volume vsrc1 =volume(vec(xcenp+dxmaxp,ycenp-dymaxp,zcenp+dzmaxp), vec(xcenp+dxminp,ycenp+dymaxp,zcenp-dzmaxp));
      const volume vsrc2 =volume(vec(xcenp+dxmaxp,ycenp-dymaxp,zcenp-dzminp), vec(xcenp-dxmaxp,ycenp+dymaxp,zcenp-dzmaxp));
      const volume vsrc3 =volume(vec(xcenp-dxmaxp,ycenp-dymaxp,zcenp+dzmaxp), vec(xcenp-dxminp,ycenp+dymaxp,zcenp-dzmaxp));
      const volume vsrc4 =volume(vec(xcenp+dxmaxp,ycenp-dymaxp,zcenp+dzmaxp), vec(xcenp-dxmaxp,ycenp+dymaxp,zcenp+dzminp));
      f.add_volume_source(Ez, src, vsrc1, cdouble(-amplitude,0));
      f.add_volume_source(Ex, src, vsrc2, cdouble(-amplitude,0));
      f.add_volume_source(Ez, src, vsrc3, cdouble(amplitude,0));
      f.add_volume_source(Ex, src, vsrc4, cdouble(amplitude,0));

      ycenp=ycenp+insulation_thickness_p+winding_thickness_p;     
    }
    //void add_point_source(component c, double freq, double width, double peaktime, double cutoff, const vec &, std::complex<double> amp = 1.0, int is_continuous = 0);
    //void add_volume_source(component c, const src_time &src, const volume &, std::complex<double> amp = 1.0);

 
    xcenp=xcen; 
    ycenp=ycen-dymax; 
    zcenp=zcen;
    volume box1( vec(xcenp+dxminp,ycenp-dymaxp,zcenp+dzminp), vec(xcenp-dxminp,ycenp+dymaxp,zcenp-dzminp) );
    ycenp=ycen+dymax; 
    volume box2( vec(xcenp+dxminp,ycenp-dymaxp,zcenp+dzminp), vec(xcenp-dxminp,ycenp+dymaxp,zcenp-dzminp) );
      
    cout<<"Start FFT1"<<endl;
    double fmin = fcen-df, fmax = fcen+df;
    int Nfreq = 10000;
    //dft_flux flux1 = f.add_dft_flux_box(box1,  fmin, fmax, Nfreq);
    cout<<"Start FFT2"<<endl;
    //dft_flux flux2 = f.add_dft_flux_box(box2, fmin, fmax, Nfreq);
    cout<<"Start FFT3"<<endl;
    double init_energy = f.field_energy_in_box(box1);
    cout<<"Start FFT4"<<endl;
   

  /*vec lb(vec()), rb(vec());
  vec lt(vec()), rt(vec());
  flux_vol *left = f.add_flux_plane(lb, lt);
  flux_vol *right = f.add_flux_plane(rb, rt);
  flux_vol *bottom = f.add_flux_plane(lb, rb);
  flux_vol *top = f.add_flux_plane(lt, rt);
  long double fluxL = 0;
  */
  //integral of flux = change in energy of box
  f.step();
      cout<<"Started"<<endl;
    volume vxy=volume(vec(-xsize,-ysize,0),vec(xsize,ysize,0));
    volume vxz=volume(vec(-xsize,0,-zsize),vec(xsize,0,zsize));
    volume vyz=volume(vec(0,-ysize,-zsize),vec(0,ysize,zsize));
    
    h5file * fEx=f.open_h5file("fEx",h5file::WRITE,0,false);    
    h5file * fEy=f.open_h5file("fEy",h5file::WRITE,0,false);
    h5file * fEz=f.open_h5file("fEz",h5file::WRITE,0,false);
    h5file * fDx=f.open_h5file("fDx",h5file::WRITE,0,false);    
    h5file * fDy=f.open_h5file("fDy",h5file::WRITE,0,false);
    h5file * fDz=f.open_h5file("fDz",h5file::WRITE,0,false);
    h5file * fHx=f.open_h5file("fHx",h5file::WRITE,0,false);    
    h5file * fHy=f.open_h5file("fHy",h5file::WRITE,0,false);
    h5file * fHz=f.open_h5file("fHz",h5file::WRITE,0,false);
    h5file * fBx=f.open_h5file("fBx",h5file::WRITE,0,false);    
    h5file * fBy=f.open_h5file("fBy",h5file::WRITE,0,false);
    h5file * fBz=f.open_h5file("fBz",h5file::WRITE,0,false);
    
    f.output_hdf5(Dielectric,vyz,fEx);
    f.output_hdf5(Permeability,vyz,fEy);
    f.output_hdf5(Ez,vyz,fEz);
    f.output_hdf5(Dx,vyz,fDx);
    f.output_hdf5(Dy,vyz,fDy);
    f.output_hdf5(Dz,vyz,fDz);
    f.output_hdf5(Hx,vyz,fHx);
    f.output_hdf5(Hy,vyz,fHy);
    f.output_hdf5(Hz,vyz,fHz);
    f.output_hdf5(Bx,vyz,fBx);
    f.output_hdf5(By,vyz,fBy);
    f.output_hdf5(Bz,vyz,fBz);



    for(int i=1;i<=1000;i++)
    {
      f.step();
          //fluxL += f.dt * (left->flux() - right->flux() + bottom->flux() - top->flux());

       if ((i%10)==0)
      {

    f.output_hdf5(Hx,vyz);
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

        cdouble Vm=compute_Im(f,ycen);
        cdouble Im=compute_Im(f,ycen);
        cdouble Ve=compute_Im(f,ycen);
        cdouble Ie=compute_Im(f,ycen);

        monitor_point pin;
        f.get_point(&pin, vec(xcen,ycen-dymin-0.5*wcore,zcen));
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

        monitor_point po;
        f.get_point(&po, vec(xcen,ycen+dymin+0.5*wcore,zcen));
        cdouble E1o = po.get_component(Ex);
        cdouble E2o = po.get_component(Ey);
        cdouble E3o = po.get_component(Ez);
        cdouble D1o = po.get_component(Dx);
        cdouble D2o = po.get_component(Dy);
        cdouble D3o = po.get_component(Dz);
        cdouble H1o = po.get_component(Hx);
        cdouble H2o = po.get_component(Hy);
        cdouble H3o = po.get_component(Hz);
        cdouble B1o = po.get_component(Bx);
        cdouble B2o = po.get_component(By);
        cdouble B3o = po.get_component(Bz);        
        
        FieldsIn<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;
        FieldsOut<<H1o.real()<<" , "<<H1o.imag()<<" , "<<H2o.real()<<" , "<<H2o.imag()<<" , "<<H3o.real()<<" , "<<H3o.imag()<<" , "<<B1o.real()<<" , "<<B1o.imag()<<" , "<<B2o.real()<<" , "<<B2o.imag()<<" , "<<B3o.real()<<" , "<<B3o.imag()<<" , "<<E1o.real()<<" , "<<E1o.imag()<<" , "<<E2o.real()<<" , "<<E2o.imag()<<" , "<<E3o.real()<<" , "<<E3o.imag()<<" , "<<D1o.real()<<" , "<<D1o.imag()<<" , "<<D2o.real()<<" , "<<D2o.imag()<<" , "<<D3o.real()<<" , "<<D3o.imag()<<endl;    
        //Hx , Hy , Hz , Bx , By , Bz , Ex , Ey , Ez , Dx , Dy , Dz 
        Time<<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
        //Im , Vm , Ie , Ve 
      }
       
    }
 

    f.output_hdf5(Dielectric,vyz,fEx);
    f.output_hdf5(Permeability,vyz,fEy);
    f.output_hdf5(Ez,vyz,fEz);
    f.output_hdf5(Dx,vyz,fDx);
    f.output_hdf5(Dy,vyz,fDy);
    f.output_hdf5(Dz,vyz,fDz);
    f.output_hdf5(Hx,vyz,fHx);
    f.output_hdf5(Hy,vyz,fHy);
    f.output_hdf5(Hz,vyz,fHz);
    f.output_hdf5(Bx,vyz,fBx);
    f.output_hdf5(By,vyz,fBy);
    f.output_hdf5(Bz,vyz,fBz);


 	/*
    double *fl1 = flux1.flux();
    double *fl2 = flux2.flux();
    cout<<"Flux Harmonics"<<endl;
    for (int i = 0; i < Nfreq; ++i) {
      Fluxes<<(fmin + i * flux1.dfreq)<<" , "<<fl1[i]<<" , "<<fl2[i]<<endl;
      //freq , fluxin , fluxout
    } 
	*/

    cout<<"SpaceEvolution"<<endl;
    //for (double y=(ycen-dymin);y<=(ycen+dymin);y=y+0.001)
    for (double y=(ycen-dymax);y<=(ycen+dymax);y=y+0.1)  
    {
      cdouble Im=compute_Im(f,y);
      cdouble Ie=compute_Ie(f,y);
      cdouble Vm=compute_Vm(f,y);
      cdouble Ve=compute_Ve(f,y);
      Space <<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
      //Im , Vm , Ie , Ve
    }
     



    Time.close();
    Space.close();
    FieldsIn.close();
    FieldsOut.close();
    Fluxes.close();
    return 0;
}
