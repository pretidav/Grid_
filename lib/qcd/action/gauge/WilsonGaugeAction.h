/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

namespace Grid {
namespace QCD {

////////////////////////////////////////////////////////////////////////
// Wilson Gauge Action .. should I template the Nc etc..
////////////////////////////////////////////////////////////////////////
template <class Gimpl>
class WilsonGaugeAction : public Action<typename Gimpl::GaugeField> {
 public:  
  INHERIT_GIMPL_TYPES(Gimpl);

  /////////////////////////// constructors
  explicit WilsonGaugeAction(RealD beta_):beta(beta_){};

  virtual std::string action_name() {return "WilsonGaugeAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[WilsonGaugeAction] Beta: " << beta << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  virtual RealD S(const GaugeField &U) {
    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
    RealD vol = U._grid->gSites();
    RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
    return action;
  };

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    // not optimal implementation FIXME
    // extend Ta to include Lorentz indexes

    RealD factor = 0.5 * beta / RealD(Nc);

    //GaugeLinkField Umu(U._grid);
    GaugeLinkField dSdU_mu(U._grid);
    for (int mu = 0; mu < Nd; mu++) {
      //Umu = PeekIndex<LorentzIndex>(U, mu);

      // Staple in direction mu
      //WilsonLoops<Gimpl>::Staple(dSdU_mu, U, mu);
      //dSdU_mu = Ta(Umu * dSdU_mu) * factor;

  
      WilsonLoops<Gimpl>::StapleMult(dSdU_mu, U, mu);
      dSdU_mu = Ta(dSdU_mu) * factor;

      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
    }
  }
private:
  RealD beta;  
};




// Anisotropic Wilson gauge action
template <class Gimpl>
class WilsonGaugeActionAnisotropic : public Action<typename Gimpl::GaugeField> {
 public:  
  INHERIT_GIMPL_TYPES(Gimpl);

// Anisotropic bare coupling (Nd betas).
 private:
  std::vector<double> vBeta;
  

  /////////////////////////// constructors
  explicit WilsonGaugeActionAnisotropic(std::vector<double> vBeta_):vBeta(vBeta_){};

  virtual std::string action_name() {return "WilsonGaugeActionAnisotropic";}

  virtual std::string LogParameters(){
    std::stringstream sstream;

    //explicit loop over the volume: TODO
    for (int mu=0;mu<Nd;mu++){
      sstream << GridLogMessage << "[WilsonGaugeActionAnisotropic] Beta[" << mu <<"]:" << vBeta[mu] << std::endl;
    }
    return sstream.str();
  }


//MOVE THIS INTO A TEST, AND CHECK EXPLICITLY

virtual void implementSF(const GaugeField &U){
/*
  assert(Nd==4);  //only 4 dim.

    GridBase *grid = U._grid;
    int Tmax = U._grid->GlobalDimensions()[3];

    Lattice<iScalar<vInteger>> coor(grid);
    LatticeCoordinate(coor, 3);

    GaugeLinkField tmp(grid);
    tmp = where(coor == Lmu-1, Cprime, tmp);
    return Cshift(tmp, mu, -1); // moves towards positive mu


//SF parameters
std::vector<double> phi_SF,phiprime_SF,omega;
double eta=0, nuSF=0;

omega[0] = 1;
omega[1] = -1/2 + nuSF;
omega[2] = -1/2 - nuSF;

phi_SF[0] = eta*omega[0] - M_PI/3;
phi_SF[1] = eta*omega[1;
phi_SF[2] = eta*omega[2] + M_PI/3;

phiprime_SF[0] = -phi_SF[0] - 4/3 * M_PI;
phiprime_SF[1] = -phi_SF[2] + 2/3 * M_PI;
phiprime_SF[2] = -phi_SF[1] + 2/3 * M_PI;

ColourMatrix C,Cprime;
C=zero;
Cprime=zero;


*/

//C(0,0)=phi_SF[0];
//C(1,1)=phi_SF[1];
//C(2,2)=phi_SF[2];

//Cprime(0,0)=phiprime_SF[0];
//Cprime(1,1)=phiprime_SF[1];
//Cprime(2,2)=phiprime_SF[2];

// for(int mu=0;mu<Nd;mu++){
//      pokeLorentz(U,Umu,mu);
// }


//now insert C @ t=0 and Cprime @ T-1  <----- this have to enter into the hot, tepid and cold start configuration! 
//then set U4(x,x0=T-1)=0 
//fix the gauge force to 0 to prevent the update of timeslices t=0 and t=T-1 and then we are done!!!!! 
};





  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms








  virtual RealD S(const GaugeField &U) {

    // FIXME. I have to weight with beta before the average!!!!!
    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);

    //RealD beta_OneMinusPlaq = WilsonLoops<Gimpl>::OneMinusavgPlaquetteAnisotropy(U,Lbeta);    

    //explicit loop over volume: TODO
    RealD vol = U._grid->gSites();
    int X = U._grid->GlobalDimensions()[0];
    int Y = U._grid->GlobalDimensions()[1];
    int Z = U._grid->GlobalDimensions()[2];
    
    if (!Gimpl::isPeriodicGaugeField()) vol -= X*Y*Z;  //subtraction of the Wall slice (needed for avoiding periodic). 



   RealD action = vBeta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
  // RealD action = ( beta_OneMinusPlaq ) * (Nd * (Nd - 1.0)) * vol * 0.5;

    return action;
  };





  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    // not optimal implementation FIXME
    // extend Ta to include Lorentz indexes

    RealD factor;
 
    GaugeLinkField dSdU_mu(U._grid);
    for (int mu = 0; mu < Nd; mu++) {
    factor = 0.5 * vBeta[mu] / RealD(Nc);
  
  
      WilsonLoops<Gimpl>::StapleMult(dSdU_mu, U, mu); 
      dSdU_mu = Ta(dSdU_mu) * factor;
      
      int isDirichlet = 0;
      if (isDirichlet){
      Lattice<iScalar<vInteger>> coor(dSdU_mu.grid_);
      LatticeCoordinate(coor, 3);
      int Tmax = dSdU_mu._grid->GlobalDimensions()[3]-1;
      dSdU_mu = where((coor==(Tmax-1) || coor==0), 0.*dSdU_mu, dSdU_mu);
      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
      }

    }
  }
private:
  RealD beta;  
};


}
}

#endif
