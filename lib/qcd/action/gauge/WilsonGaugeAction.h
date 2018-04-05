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







////////////////////////////////////////////////////////////////////////
// Wilson Gauge Action .. should I template the Nc etc..
////////////////////////////////////////////////////////////////////////
template <class Gimpl>
class WilsonGaugeAnisotropicAction : public Action<typename Gimpl::GaugeField> {
 public:  
  INHERIT_GIMPL_TYPES(Gimpl);

// Anisotropic bare coupling (Nd betas).
 private:
  std::vector<double> vbeta;


  /////////////////////////// constructors
  public:
   WilsonGaugeAnisotropicAction(std::vector<double> vbeta_):vbeta(vbeta_){};

  virtual std::string action_name() {return "WilsonGaugeAnisotropicAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    for (int j=0;j<2;j++){
    sstream << GridLogMessage << "[WilsonGaugeAnisotropicAction] Beta["<< j <<"]: " << vbeta[j] << std::endl;
    }
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  virtual RealD S(const GaugeField &U) {
//    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U); //change this!

//modifications
LatticeColourMatrix tmp(U._grid);
std::vector<LatticeColourMatrix> Link(Nd, U._grid);
tmp=zero;
std::vector<RealD> PlaqDir(2); 
PlaqDir[0]=0.0;
PlaqDir[1]=0.0;

for (int mu=0;mu<Nd;mu++){
  Link[mu] = PeekIndex<LorentzIndex>(U, mu);
}
 
for (int mu=1;mu<Nd-1;mu++){
  for (int nu=0;nu<mu;nu++){
    WilsonLoops<Gimpl>::dirPlaquette(tmp, Link, mu, nu);
    PlaqDir[0]+=TensorRemove(sum(trace(tmp))).real();
  }
}

tmp=zero;
for (int nu=0;nu<Nd-1;nu++){
  WilsonLoops<Gimpl>::dirPlaquette(tmp, Link, 3, nu);
  PlaqDir[1]+=TensorRemove(sum(trace(tmp))).real();
}
//end modifications
int T   = U._grid->GlobalDimensions()[3];
int X   = U._grid->GlobalDimensions()[0];
int Y   = U._grid->GlobalDimensions()[1];
int Z   = U._grid->GlobalDimensions()[2];
    RealD vol = U._grid->gSites();
    RealD action = 0;
    RealD Norm=((Nd-1) * (Nd-2) * vol * Nc / 2 ); 
    
    for (int j=0;j<2;j++){
      if (j==1) vol-=X*Y*Z; 
      PlaqDir[j]*=1./((Nd-1) * (Nd-2) * vol * Nc / 2 );   
      action += vbeta[j] * (1.0 - PlaqDir[j]) *  ((Nd-1) * (Nd-2)) * vol * 0.5 ;
      std::cout<< "PlaqDir=" << PlaqDir[j]<<std::endl;
      std::cout<< "action="  << action <<std::endl;
    }

    return action;
  };

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    // not optimal implementation FIXME
    // extend Ta to include Lorentz indexes

    RealD factor;
 
    GaugeLinkField dSdU_mu(U._grid);
    for (int mu = 0; mu < Nd; mu++) {
    factor = 0.5 * vbeta[mu] / RealD(Nc);
  
  
      WilsonLoops<Gimpl>::StapleMult(dSdU_mu, U, mu); 
      dSdU_mu = Ta(dSdU_mu) * factor;
      
      int isDirichlet = 0;
      if (isDirichlet){
      Lattice<iScalar<vInteger>> coor(dSdU_mu._grid);
      LatticeCoordinate(coor, 3);
      int Tmax = dSdU_mu._grid->GlobalDimensions()[3]-1;
      dSdU_mu = where((coor==(Tmax-1) || coor==0), 0.*dSdU_mu, dSdU_mu);
      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
      }
  }}
private:
  RealD beta;  
};

}
}

#endif
