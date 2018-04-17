/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/plaquette.h

Copyright (C) 2017

Author: David Preti <david.preti@to.infn.it>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */




#ifndef HMC_DSDETA_H
#define HMC_DSDETA_H

namespace Grid {
namespace QCD {



struct dSdetaParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(dSdetaParameters,
      int, interval,
      RealD, betaT,
      RealD, ct_SF);  

    dSdetaParameters(int interval = 1):
        interval(interval){}

    template <class ReaderClass >
      dSdetaParameters(Reader<ReaderClass>& Reader){
        read(Reader, "uSF_Measurement", *this);
  }
};

// this is only defined for a gauge theory
template <class Impl>
class dSdeta : public HmcObservable<typename Impl::Field> {
      dSdetaParameters Pars;
 public:
  // here forces the Impl to be of gauge fields
  // if not the compiler will complain
  INHERIT_GIMPL_TYPES(Impl);

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;


    dSdeta(int interval = 1):
        Pars(interval){}
    
    dSdeta(dSdetaParameters P):Pars(P){
        std::cout << GridLogDebug << "Creating uSF coupling " << std::endl;
    }

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {


   int T   = U._grid->GlobalDimensions()[3];
   int X   = U._grid->GlobalDimensions()[0];
   int Y   = U._grid->GlobalDimensions()[1];
   int Z   = U._grid->GlobalDimensions()[2];

   assert((X==Y)&&(Y==Z)); 
   RealD theta = 1./3. * M_PI * 1./X/X;
   RealD norm = 12.*(X*X)*( std::sin(theta) + std::sin(2.*theta) ); 

  if (traj%Pars.interval == 0){
    RealD u = WilsonLoops<Impl>::dSdeta(U,Pars.betaT,Pars.ct_SF);

    int def_prec = std::cout.precision();

    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "uSF: [ " << traj << " ] "<< norm/u << std::endl;
    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "k:   [ " << traj << " ] "<< norm << std::endl;

    std::cout.precision(def_prec);
  }

  }
};


 



}  // namespace QCD
}  // namespace Grid

#endif  // HMC_DSDETA_H
