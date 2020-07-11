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

double eps(const vec &) { return 1000.0; }


int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  const char *mydirname = "MMTL-out";

  double a = 5.0;
  double ttot = 1000.0;

  double fcen = 0.1; // ; pulse center frequency
  double df = 0.0;//00001*(fcen);    // ; df
  //continuous_src_time src(cdouble(fcen,0));
  continuous_src_time src(cdouble(fcen,0.0));


  grid_volume gv = vol3d(4.0, 4.0, 4.0, a);
  gv.center_origin();
  structure s(gv, eps, pml(1.0));
  s.set_output_directory(mydirname);

  fields f(&s);
  const volume vsrc1 =volume(vec(0.0,0.0,0.0), vec(0.0,0.0,0.0));
  f.add_volume_source(Ex, src, vsrc1, cdouble(1000.0,0.0));
  
    volume vxy=volume(vec(-2.0,-2.0,0.0),vec(2.0,2.0,0.0));
 
 std::ofstream FieldsIn;
    FieldsIn.open ("FieldEvolutionIn.txt");
   

  while (f.time() < ttot) {
        f.step();

        /*monitor_point pin;
        f.get_point(&pin, vec(0.5,0.00,0.00));
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
        */
  }

        f.output_hdf5(Ex,vxy);
        f.output_hdf5(Ey,vxy);
        f.output_hdf5(Ez,vxy);
        f.output_hdf5(Dx,vxy);
        f.output_hdf5(Dy,vxy);
        f.output_hdf5(Dz,vxy);
        f.output_hdf5(Hx,vxy);
        f.output_hdf5(Hy,vxy);
        f.output_hdf5(Hz,vxy);
        f.output_hdf5(Bx,vxy);
        f.output_hdf5(By,vxy);
        f.output_hdf5(Bz,vxy);
        f.output_hdf5(Sx,vxy);
        f.output_hdf5(Sy,vxy);
        f.output_hdf5(Sz,vxy);

  return 0;
}
