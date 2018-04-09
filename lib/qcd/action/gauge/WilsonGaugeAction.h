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
class WilsonGaugeSFAction : public Action<typename Gimpl::GaugeField> {
 public:  
  INHERIT_GIMPL_TYPES(Gimpl);

// Anisotropic bare coupling (Nd betas).
 private:
  std::vector<double> vbeta;
  RealD cs_SF, ct_SF;


  /////////////////////////// constructors
  public:
   WilsonGaugeSFAction(std::vector<double> vbeta_, RealD cs_SF_, RealD ct_SF_): vbeta(vbeta_), cs_SF(cs_SF_), ct_SF(ct_SF_){};

  virtual std::string action_name() {return "WilsonGaugeSFAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    for (int j=0;j<2;j++){
    sstream << GridLogMessage << "[WilsonGaugeSFAction] Beta["<< j <<"]: " << vbeta[j] << std::endl;
    }
    sstream << GridLogMessage << "[WilsonGaugeSFAction] ct = " << ct_SF << std::endl;
    sstream << GridLogMessage << "[WilsonGaugeSFAction] cs = " << cs_SF << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  virtual RealD S(const GaugeField &U) {

int T   = U._grid->GlobalDimensions()[3];
int X   = U._grid->GlobalDimensions()[0];
int Y   = U._grid->GlobalDimensions()[1];
int Z   = U._grid->GlobalDimensions()[2];

LatticeColourMatrix tmp(U._grid), tmpBLK(U._grid),tmpBC(U._grid);
std::vector<LatticeColourMatrix> Link(Nd, U._grid);
tmp=zero;
tmpBLK=zero;
tmpBC=zero;

std::vector<RealD> PlaqDirBLK(2); 
PlaqDirBLK[0]=0.0;
PlaqDirBLK[1]=0.0;
std::vector<RealD> PlaqDirBC(2); 
PlaqDirBC[0]=0.0;
PlaqDirBC[1]=0.0;


for (int mu=0;mu<Nd;mu++){
  Link[mu] = PeekIndex<LorentzIndex>(U, mu);
}
 
for (int mu=1;mu<Nd-1;mu++){
  for (int nu=0;nu<mu;nu++){
    WilsonLoops<Gimpl>::dirPlaquette(tmp, Link, mu, nu);
  Lattice<iScalar<vInteger>> coorS(tmp._grid);
  LatticeCoordinate(coorS, 3);  

    tmpBC   = where( ((coorS==0) || (coorS==(T-1))), tmp, 0.*tmp); 
    tmpBLK  = where( ((coorS==0) || (coorS==(T-1))), 0.*tmp, tmp); 
    PlaqDirBLK[0]+=TensorRemove(sum(trace(tmpBLK))).real();
    PlaqDirBC[0] +=TensorRemove(sum(trace(tmpBC))).real();
  }
}

for (int nu=0;nu<Nd-1;nu++){
  WilsonLoops<Gimpl>::dirPlaquette(tmp, Link, 3, nu);
  Lattice<iScalar<vInteger>> coorT(tmp._grid);
  LatticeCoordinate(coorT, 3);  

  tmpBC   = where( ((coorT==0) || (coorT==(T-2))), tmp, 0.*tmp); //no faces moltiplicity ?
  tmpBLK  = where( ((coorT==0) || (coorT==(T-2))), 0.*tmp, tmp);
  PlaqDirBLK[1]+=TensorRemove(sum(trace(tmpBLK))).real();
  PlaqDirBC[1] +=TensorRemove(sum(trace(tmpBC))).real();
}

RealD vol = U._grid->gSites();
RealD action = 0;
RealD NORM; 
  for (int j=0;j<2;j++){
    PlaqDirBLK[j]*=1./(Nc);
    PlaqDirBC[j] *=1./(Nc); 
    if (j==0) NORM=0.5*cs_SF;
    if (j==1) NORM=ct_SF;
    action += vbeta[j] * (((Nd-1) * (Nd-2))*0.5 * (X*Y*Z*T-2*X*Y*Z - (j)*X*Y*Z)  - PlaqDirBLK[j]); //this modification in the volume account for the missing temporal plaqs @ x4=T
    action += vbeta[j] * (((Nd-1) * (Nd-2))*0.5 * (2*X*Y*Z)  - PlaqDirBC[j]) *NORM;
  }
  return action;
};

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {

    std::vector<RealD> factor(2);
    std::vector<LatticeColourMatrix> dSdU_mu(Nd, U._grid);
    std::vector<LatticeColourMatrix> U_mu(Nd, U._grid);
    LatticeColourMatrix staple(U._grid);
    LatticeColourMatrix tmp1(U._grid),tmp2(U._grid);

    for (int mu = 0; mu < Nd; mu++) {
      U_mu[mu] = PeekIndex<LorentzIndex>(U, mu);
    }

    staple = zero;    
    factor[0] = 0.5 * vbeta[0] / RealD(Nc);
    factor[1] = 0.5 * vbeta[1] / RealD(Nc);
    
    Lattice<iScalar<vInteger>> coorStaple(staple._grid);
    LatticeCoordinate(coorStaple, 3);
    int Time = staple._grid->GlobalDimensions()[3];

 for (int mu=0;mu<Nd-1;mu++){
   staple=zero;
  for (int nu=0; nu<Nd; nu++) {
      if (nu != mu) {
        tmp1 = Cshift(U_mu[nu], mu, 1);
        tmp2 = Cshift(U_mu[mu], nu, 1);
        if (nu==3){ 
          staple += tmp1* adj(U_mu[nu]*tmp2)*factor[1];
          staple = where((coorStaple==Time-2), ct_SF*staple, staple);
          std::cout << "ciao" << std::endl; 
        } else {
          staple += tmp1* adj(U_mu[nu]*tmp2)*factor[0];
         }
        tmp2 = adj(U_mu[mu]*tmp1)*U_mu[nu];
        if (nu==3){
        staple += Cshift(tmp2, nu, -1)*factor[1];
        staple = where((coorStaple==1), ct_SF*staple, staple);
        } else {
        staple += Cshift(tmp2, nu, -1)*factor[0];
        }
      }
  }
  dSdU_mu[mu] = U_mu[mu]*staple;
  dSdU_mu[mu] = Ta(dSdU_mu[mu]);
 }

  staple=zero;
  for (int nu=0;nu<Nd-1;nu++){
    tmp1 = Cshift(U_mu[nu], 3, 1);
    tmp2 = Cshift(U_mu[3], nu, 1);
    staple += tmp1* adj(U_mu[nu]*tmp2);
    tmp2 = adj(U_mu[3]*tmp1)*U_mu[nu];
    staple += Cshift(tmp2, nu, -1);
  }
  staple = where((coorStaple==0 || coorStaple==Time-2), ct_SF*staple, staple); //only T-2 slice is relevant here.
  dSdU_mu[3] = U_mu[3]*staple;
  dSdU_mu[3] = Ta(dSdU_mu[3]) * factor[1];  //all those staples carry beta[1] 
  
  Lattice<iScalar<vInteger>> coor(dSdU_mu[3]._grid);
  LatticeCoordinate(coor, 3);
  int T = dSdU_mu[3]._grid->GlobalDimensions()[3];
  for (int mu=0;mu<Nd;mu++){
    dSdU_mu[mu] = where((coor==0 || coor==T-1), 0.*dSdU_mu[mu], dSdU_mu[mu]); 
    PokeIndex<LorentzIndex>(dSdU, dSdU_mu[mu], mu);
  }
 }
//private:
//  RealD beta;  
};

}
}

#endif
