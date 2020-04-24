/* Copyright (C) 2005-2019 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include "meep.hpp"
#include "ctl-math.h"
#include <mpb.h>
#include "ctlgeom.h"
#include "meepgeom.hpp"
#include <math.h>
#ifndef DATADIR
#define DATADIR "./"
#endif

using namespace meep;
using namespace std;
typedef std::complex<double> cdouble;

double eps(const vec &) { return 1.0; }


int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  const char *mydirname = "MMTL-out";

  double a = 10.0;
  double ttot = 170.0;

  double fcen = (3e9)/(3e12); // ; pulse center frequency
  double df = 0.000001*((3e9)/(3e12));    // ; df
  //continuous_src_time src(cdouble(fcen,0));
  gaussian_src_time src(fcen,df);


  grid_volume gv = vol3d(2.0, 2.0, 2.0, a);
  gv.center_origin();
  structure s(gv, eps);
  s.set_output_directory(mydirname);

  fields f(&s);
  const volume vsrc1 =volume(vec(0.0,0.0,0.0), vec(0.0,0.0,0.0));
  f.add_volume_source(Ex, src, vsrc1, cdouble(1.0,0.0));
  
    volume vxy=volume(vec(-1.0,-1.0,0.0),vec(1.0,1.0,0.0));
 
 std::ofstream FieldsIn;
    FieldsIn.open ("FieldEvolutionIn.txt");
   

  while (f.time() < ttot) {
    f.step();

        monitor_point pin;
        f.get_point(&pin, vec(0.01,0.01,0.01));
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
        FieldsIn<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;


  }



  return 0;
}
