/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/utils/WilsonLoops.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef QCD_UTILS_WILSON_LOOPS_H
#define QCD_UTILS_WILSON_LOOPS_H
namespace Grid {
namespace QCD {

// Common wilson loop observables
template <class Gimpl> class WilsonLoops : public Gimpl {
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  //////////////////////////////////////////////////
  // directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    // Annoyingly, must use either scope resolution to find dependent base
    // class,
    // or this-> ; there is no "this" in a static method. This forces explicit
    // Gimpl scope
    // resolution throughout the usage in this file, and rather defeats the
    // purpose of deriving
    // from Gimpl.
    /*
    plaq = Gimpl::CovShiftBackward(
        U[mu], mu, Gimpl::CovShiftBackward(
                       U[nu], nu, Gimpl::CovShiftForward(U[mu], mu, U[nu])));
                       */
    // _
    //|< _|
    plaq = Gimpl::CovShiftForward(U[mu],mu,
           Gimpl::CovShiftForward(U[nu],nu,
           Gimpl::CovShiftBackward(U[mu],mu,
           Gimpl::CovShiftIdentityBackward(U[nu], nu))));
  }

  
  static void WilsonLoopNxM(GaugeMat &plaq, const std::vector<GaugeMat> &U,
			    const int mu, const int nu,const int Nmu, const int Mnu) {

    //  ^
    //  |mu 
    //  | 
    //  -----> nu 

    GaugeMat tmp1(U[0]._grid), tmp2(U[0]._grid);
    tmp1=Gimpl::CovShiftIdentityBackward(U[nu],nu); 

    for (int i=0;i<Mnu-1;i++){
      tmp2=Gimpl::CovShiftBackward(U[nu],nu,tmp1);
      tmp1=tmp2;
    }
    for (int j=0;j<Nmu;j++){
      tmp2=Gimpl::CovShiftBackward(U[mu],mu,tmp1);
      tmp1=tmp2;
    }
 for (int i=0;i<Mnu;i++){
      tmp2=Gimpl::CovShiftForward(U[nu],nu,tmp1);
      tmp1=tmp2;
    }
 for (int j=0;j<Nmu;j++){
      tmp2=Gimpl::CovShiftForward(U[mu],mu,tmp1);
      tmp1=tmp2;
    }
   plaq=tmp1;
  }
  
  static void dirPlaquette2x2(GaugeMat &plaq, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    

    //  P(y)=
    //  x->x->x
    //  |     |
    //  x  y  x
    //  |     |
    //  x<-x<-x
    //
    //  ^
    //  |mu 
    //  | 
    //  -----> nu 
    

    GaugeMat tmp(plaq._grid);
    tmp = Gimpl::CovShiftForward(U[mu],mu,Gimpl::CovShiftForward(U[mu],mu,
           Gimpl::CovShiftForward(U[nu],nu,Gimpl::CovShiftForward(U[nu],nu,
           Gimpl::CovShiftBackward(U[mu],mu,Gimpl::CovShiftBackward(U[mu],mu,
           Gimpl::CovShiftBackward(U[nu],nu,Gimpl::CovShiftIdentityBackward(U[nu], nu))))))));
     
     plaq = Cshift(Cshift(plaq, mu, -1), nu, -1);          
  }

  //////////////////////////////////////////////////
  // trace of directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceDirPlaquette(ComplexField &plaq,
                                const std::vector<GaugeMat> &U, const int mu,
                                const int nu) {
    GaugeMat sp(U[0]._grid);
    dirPlaquette(sp, U, mu, nu);
    plaq = trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of plaquette
  //////////////////////////////////////////////////
  static void sitePlaquette(ComplexField &Plaq,
                            const std::vector<GaugeMat> &U) {
    ComplexField sitePlaq(U[0]._grid);
    Plaq = zero;
    for (int mu = 1; mu < Nd; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirPlaquette(sitePlaq, U, mu, nu);
        Plaq = Plaq + sitePlaq;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumPlaquette(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu._grid);
    // inefficient here
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    ComplexField Plaq(Umu._grid);

    sitePlaquette(Plaq, U);
    auto Tp = sum(Plaq);
    auto p = TensorRemove(Tp);
    return p.real();
  }


  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgPlaquette(const GaugeLorentz &Umu) {
    RealD sumplaq = sumPlaquette(Umu);
    double vol = Umu._grid->gSites();
    double faces = (1.0 * Nd * (Nd - 1)) / 2.0;
    return sumplaq / vol / faces / Nc; // Nd , Nc dependent... FIXME
  }

  static Real avgPlaquetteSF(const GaugeLorentz &Umu) {
    int ndim = Umu._grid->_ndimension;
    int X   = Umu._grid->GlobalDimensions()[0];
    int Y   = Umu._grid->GlobalDimensions()[1];
    int Z   = Umu._grid->GlobalDimensions()[2];
    int T   = Umu._grid->GlobalDimensions()[3];

    LatticeColourMatrix tmp(Umu._grid), tmpBLK(Umu._grid),tmpBC(Umu._grid);
    std::vector<LatticeColourMatrix> Link(Nd, Umu._grid);
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
      Link[mu] = PeekIndex<LorentzIndex>(Umu, mu);
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

tmpBC=zero;
tmpBLK=zero;
tmp=zero;

    for (int nu=0;nu<Nd-1;nu++){
      WilsonLoops<Gimpl>::dirPlaquette(tmp, Link, 3, nu);
      Lattice<iScalar<vInteger>> coorT(tmp._grid);
      LatticeCoordinate(coorT, 3);

      tmpBC   = where( ((coorT==0) || (coorT==(T-2))), tmp, 0.*tmp);
      tmpBLK  = where( ((coorT==0) || (coorT==(T-2))), 0.*tmp, tmp);
      PlaqDirBLK[1]+=TensorRemove(sum(trace(tmpBLK))).real();
      PlaqDirBC[1] +=TensorRemove(sum(trace(tmpBC))).real();
    }

    RealD vol = Umu._grid->gSites();
    RealD sumplaq = 0.;
    RealD NORM;
    for (int j=0;j<2;j++){
           // PlaqDirBLK[j]*=1./(Nc);
           // PlaqDirBC[j] *=1./(Nc);
      if (j==0) NORM=0.5;
      if (j==1) NORM=1.0;
      //std::cout << "plaqDirBLK[" << j << "]: " << PlaqDirBLK[j] << std::endl;
      //std::cout << "plaqDirBC[" << j << "]: " << PlaqDirBC[j] << std::endl; 
      sumplaq += double(PlaqDirBLK[j]);      // ( (Nd-1)*(Nd-2)*0.5 * (X*Y*Z*T-2.*X*Y*Z - (j)*X*Y*Z) ); 
      //std::cout << "j: " << j << " -> " << sumplaq << std::endl;
      sumplaq += double(PlaqDirBC[j] * NORM);       //  ( (Nd-1)*(Nd-2)*0.5 * (2.*X*Y*Z) ) ;
      //std::cout << "j: " << j << " -> " << sumplaq << std::endl;
    }
    return sumplaq / Nc / (X*Y*Z*(T-1)*Nd*(Nd-1)*0.5);   //( (Nd-1)*(Nd-2)* (X*Y*Z*(T-1)))/Nc   ; 
    //is used to normalize the averaged plaq to be =1 for trivial cnfg.
  }
  
  //////////////////////////////////////////////////
  // average over all x,y,z the temporal loop
  //////////////////////////////////////////////////
  static ComplexD avgPolyakovLoop(const GaugeField &Umu) {  //assume Nd=4
    GaugeMat Ut(Umu._grid), P(Umu._grid);
    ComplexD out;
    int T = Umu._grid->GlobalDimensions()[3];
    int X = Umu._grid->GlobalDimensions()[0];
    int Y = Umu._grid->GlobalDimensions()[1];
    int Z = Umu._grid->GlobalDimensions()[2];

    Ut = peekLorentz(Umu,3); //Select temporal direction
    P = Ut;
    for (int t=1;t<T;t++){ 
      P = Gimpl::CovShiftForward(Ut,3,P);
    }
   RealD norm = 1.0/(Nc*X*Y*Z*T);
   out = sum(trace(P))*norm;
   return out;   
}

  //////////////////////////////////////////////////
  // average over traced single links
  //////////////////////////////////////////////////
  static RealD linkTrace(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu._grid);

    ComplexField Tr(Umu._grid);
    Tr = zero;
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
      Tr = Tr + trace(U[mu]);
    }

    auto Tp = sum(Tr);
    auto p = TensorRemove(Tp);

    double vol = Umu._grid->gSites();

    return p.real() / vol / 4.0 / 3.0;
  };


  //////////////////////////////////////////////////
  //           g2SF= k / dSdeta                   //
  //////////////////////////////////////////////////
  static RealD dSdeta(const GaugeLorentz &Umu, RealD BetaT, RealD CT) {
    std::vector<GaugeMat> U(Nd, Umu._grid);
    ColourMatrix lambda8;
    RealD out;
    out=0.0;

    int T   = Umu._grid->GlobalDimensions()[3];
    int X   = Umu._grid->GlobalDimensions()[0];
    int Y   = Umu._grid->GlobalDimensions()[1];
    int Z   = Umu._grid->GlobalDimensions()[2];

    lambda8=zero;
    lambda8()()(0,0) = 1.0;
    lambda8()()(1,1) = -0.5;
    lambda8()()(2,2) = -0.5;
    std::vector <LatticeColourMatrix> E8(Nd-1, Umu._grid), E8prime(Nd-1, Umu._grid), B(Nd, Umu._grid);  
    Lattice<iScalar<vInteger>> coor(Umu._grid);
    LatticeCoordinate(coor, 3);  

    for (int mu=0;mu<Nd;mu++){
     U[mu]=peekLorentz(Umu,mu); 
     B[mu]=ComplexD(0.0, 1.)*lambda8*U[mu];
    }

    for (int mu=0;mu<Nd-1;mu++){
      E8[mu] =      Gimpl::CovShiftForward(B[mu],mu,
                    Gimpl::CovShiftForward(U[3],3,
                    Gimpl::CovShiftBackward(U[mu],mu,
                    Gimpl::CovShiftIdentityBackward(U[3], 3))));
      E8[mu] = where(coor==0, E8[mu], 0.*E8[mu]);
    
std::cout << "E8[" << mu << "]: " << TensorRemove(sum(trace(E8[mu]))).real() << std::endl;

      E8prime[mu] = Gimpl::CovShiftForward(U[3],3,
                    Gimpl::CovShiftForward(B[mu],mu,
                    Gimpl::CovShiftBackward(U[3],3,
                    Gimpl::CovShiftIdentityBackward(U[mu], mu))));
      E8prime[mu] = where(coor==T-2, E8prime[mu], 0.*E8prime[mu]);
      

 std::cout << "E8prime[" << mu << "]: " << TensorRemove(sum(trace(E8prime[mu]))).real() << std::endl;

     out+= TensorRemove(sum(trace(E8[mu]))).real();
     out-= TensorRemove(sum(trace(E8prime[mu]))).real();
    }

    RealD norm = - BetaT*CT/(Nc*(T-1));   
    return norm*out;
  };


  //////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
                     int nu) {

    GridBase *grid = Umu._grid;

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = zero;

    if (nu != mu) {

      // mu
      // ^
      // |__>  nu

      //    __
      //      |
      //    __|
      //

      staple += Gimpl::ShiftStaple(
          Gimpl::CovShiftForward(
              U[nu], nu,
              Gimpl::CovShiftBackward(
                  U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
          mu);

      //  __
      // |
      // |__
      //
      //
      staple += Gimpl::ShiftStaple(
          Gimpl::CovShiftBackward(U[nu], nu,
                                  Gimpl::CovShiftBackward(U[mu], mu, U[nu])),
          mu);
    }
  }


// For the force term
static void StapleMult(GaugeMat &staple, const GaugeLorentz &Umu, int mu) {
    GridBase *grid = Umu._grid;
    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      // this operation is taking too much time
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = zero;
    GaugeMat tmp1(grid);
    GaugeMat tmp2(grid);

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
        // this is ~10% faster than the Staple
        tmp1 = Cshift(U[nu], mu, 1);
        tmp2 = Cshift(U[mu], nu, 1);
        staple += tmp1* adj(U[nu]*tmp2);
        tmp2 = adj(U[mu]*tmp1)*U[nu];
        staple += Cshift(tmp2, nu, -1);
      }
    }
    staple = U[mu]*staple;
}

  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple, const GaugeLorentz &Umu, int mu) {

    GridBase *grid = Umu._grid;

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = zero;

    for (int nu = 0; nu < Nd; nu++) {

      if (nu != mu) {

        // mu
        // ^
        // |__>  nu

        //    __
        //      |
        //    __|
        //
     
        staple += Gimpl::ShiftStaple(
            Gimpl::CovShiftForward(
                U[nu], nu,
                Gimpl::CovShiftBackward(
                    U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
            mu);

        //  __
        // |
        // |__
        //
        //

        staple += Gimpl::ShiftStaple(
            Gimpl::CovShiftBackward(U[nu], nu,
                                    Gimpl::CovShiftBackward(U[mu], mu, U[nu])), mu);
      }
    }
  }

  //////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu, upper part
  //////////////////////////////////////////////////
  static void StapleUpper(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
                          int nu) {
    if (nu != mu) {
      GridBase *grid = Umu._grid;

      std::vector<GaugeMat> U(Nd, grid);
      for (int d = 0; d < Nd; d++) {
        U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
      }

      // mu
      // ^
      // |__>  nu

      //    __
      //      |
      //    __|
      //

      staple = Gimpl::ShiftStaple(
          Gimpl::CovShiftForward(
              U[nu], nu,
              Gimpl::CovShiftBackward(
                  U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
          mu);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu, lower part
  ////////////////////////////////////////////////////////////////////////
  static void StapleLower(GaugeMat &staple, const GaugeLorentz &Umu, int mu, int nu) {
    if (nu != mu) {
      GridBase *grid = Umu._grid;

      std::vector<GaugeMat> U(Nd, grid);
      for (int d = 0; d < Nd; d++) {
        U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
      }

      // mu
      // ^
      // |__>  nu

      //  __
      // |
      // |__
      //
      //
      staple = Gimpl::ShiftStaple(
          Gimpl::CovShiftBackward(U[nu], nu,
                                  Gimpl::CovShiftBackward(U[mu], mu, U[nu])),
          mu);

    }
  }

  //////////////////////////////////////////////////////
  //  Field Strength
  //////////////////////////////////////////////////////
  static void FieldStrength(GaugeMat &FS, const GaugeLorentz &Umu, int mu, int nu){
      // Fmn +--<--+  Ut +--<--+
      //     |     |     |     |
      //  (x)+-->--+     +-->--+(x)  - h.c.
      //     |     |     |     |
      //     +--<--+     +--<--+

      GaugeMat Vup(Umu._grid), Vdn(Umu._grid);
      StapleUpper(Vup, Umu, mu, nu);
      StapleLower(Vdn, Umu, mu, nu);
      GaugeMat v = Vup - Vdn;
      GaugeMat u = PeekIndex<LorentzIndex>(Umu, mu);  // some redundant copies
      GaugeMat vu = v*u;
      //FS = 0.25*Ta(u*v + Cshift(vu, mu, -1));
      FS = (u*v + Cshift(vu, mu, -1));
      FS = 0.125*(FS - adj(FS));
  }

  static Real TopologicalCharge(GaugeLorentz &U){
    // 4d topological charge
    assert(Nd==4);
    // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
    GaugeMat Bx(U._grid), By(U._grid), Bz(U._grid);
    FieldStrength(Bx, U, Ydir, Zdir);
    FieldStrength(By, U, Zdir, Xdir);
    FieldStrength(Bz, U, Xdir, Ydir);

    // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
    GaugeMat Ex(U._grid), Ey(U._grid), Ez(U._grid);
    FieldStrength(Ex, U, Tdir, Xdir);
    FieldStrength(Ey, U, Tdir, Ydir);
    FieldStrength(Ez, U, Tdir, Zdir);

    double coeff = 8.0/(32.0*M_PI*M_PI);

    ComplexField qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
    auto Tq = sum(qfield);
    return TensorRemove(Tq).real();
  }


  //////////////////////////////////////////////////////
  // Similar to above for rectangle is required
  //////////////////////////////////////////////////////
  static void dirRectangle(GaugeMat &rect, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    rect = Gimpl::CovShiftForward(
               U[mu], mu, Gimpl::CovShiftForward(U[mu], mu, U[nu])) * // ->->|
           adj(Gimpl::CovShiftForward(
               U[nu], nu, Gimpl::CovShiftForward(U[mu], mu, U[mu])));
    rect = rect +
           Gimpl::CovShiftForward(
               U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[nu])) * // ->||
               adj(Gimpl::CovShiftForward(
                   U[nu], nu, Gimpl::CovShiftForward(U[nu], nu, U[mu])));
  }
  static void traceDirRectangle(ComplexField &rect,
                                const std::vector<GaugeMat> &U, const int mu,
                                const int nu) {
    GaugeMat sp(U[0]._grid);
    dirRectangle(sp, U, mu, nu);
    rect = trace(sp);
  }
  static void siteRectangle(ComplexField &Rect,
                            const std::vector<GaugeMat> &U) {
    ComplexField siteRect(U[0]._grid);
    Rect = zero;
    for (int mu = 1; mu < Nd; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirRectangle(siteRect, U, mu, nu);
        Rect = Rect + siteRect;
      }
    }
  }

  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumRectangle(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu._grid);

    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    ComplexField Rect(Umu._grid);

    siteRectangle(Rect, U);

    auto Tp = sum(Rect);
    auto p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgRectangle(const GaugeLorentz &Umu) {

    RealD sumrect = sumRectangle(Umu);

    double vol = Umu._grid->gSites();

    double faces = (1.0 * Nd * (Nd - 1)); // 2 distinct orientations summed

    return sumrect / vol / faces / Nc; // Nd , Nc dependent... FIXME
  }

  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void RectStapleDouble(GaugeMat &U2, const GaugeMat &U, int mu) {
    U2 = U * Cshift(U, mu, 1);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Hop by two optimisation strategy does not work nicely with Gparity. (could
  // do,
  // but need to track two deep where cross boundary and apply a conjugation).
  // Must differentiate this in Gimpl, and use Gimpl::isPeriodicGaugeField to do
  // so .
  ////////////////////////////////////////////////////////////////////////////
  static void RectStapleOptimised(GaugeMat &Stap, std::vector<GaugeMat> &U2,
                                  std::vector<GaugeMat> &U, int mu) {

    Stap = zero;

    GridBase *grid = U[0]._grid;

    GaugeMat Staple2x1(grid);
    GaugeMat tmp(grid);

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {

        // Up staple    ___ ___
        //             |       |
        tmp = Cshift(adj(U[nu]), nu, -1);
        tmp = adj(U2[mu]) * tmp;
        tmp = Cshift(tmp, mu, -2);

        Staple2x1 = Gimpl::CovShiftForward(U[nu], nu, tmp);

        // Down staple
        //             |___ ___|
        //
        tmp = adj(U2[mu]) * U[nu];
        Staple2x1 += Gimpl::CovShiftBackward(U[nu], nu, Cshift(tmp, mu, -2));

        //              ___ ___
        //             |    ___|
        //             |___ ___|
        //

        Stap += Cshift(Gimpl::CovShiftForward(U[mu], mu, Staple2x1), mu, 1);

        //              ___ ___
        //             |___    |
        //             |___ ___|
        //

        //  tmp= Staple2x1* Cshift(U[mu],mu,-2);
        //  Stap+= Cshift(tmp,mu,1) ;
        Stap += Cshift(Staple2x1, mu, 1) * Cshift(U[mu], mu, -1);
        ;

        //       --
        //      |  |
        //
        //      |  |

        tmp = Cshift(adj(U2[nu]), nu, -2);
        tmp = Gimpl::CovShiftBackward(U[mu], mu, tmp);
        tmp = U2[nu] * Cshift(tmp, nu, 2);
        Stap += Cshift(tmp, mu, 1);

        //      |  |
        //
        //      |  |
        //       --

        tmp = Gimpl::CovShiftBackward(U[mu], mu, U2[nu]);
        tmp = adj(U2[nu]) * tmp;
        tmp = Cshift(tmp, nu, -2);
        Stap += Cshift(tmp, mu, 1);
      }
    }
  }

  static void RectStaple(GaugeMat &Stap, const GaugeLorentz &Umu, int mu) {
    RectStapleUnoptimised(Stap, Umu, mu);
  }
  static void RectStaple(const GaugeLorentz &Umu, GaugeMat &Stap,
                         std::vector<GaugeMat> &U2, std::vector<GaugeMat> &U,
                         int mu) {
    if (Gimpl::isPeriodicGaugeField()) {
      RectStapleOptimised(Stap, U2, U, mu);
    } else {
      RectStapleUnoptimised(Stap, Umu, mu);
    }
  }

  static void RectStapleUnoptimised(GaugeMat &Stap, const GaugeLorentz &Umu,
                                    int mu) {
    GridBase *grid = Umu._grid;

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }

    Stap = zero;

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
        //           __ ___
        //          |    __ |
        //
        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftForward(
                U[mu], mu,
                Gimpl::CovShiftForward(
                    U[nu], nu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu,
                        Gimpl::CovShiftBackward(
                            U[mu], mu,
                            Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
            mu);

        //              __
        //          |__ __ |

        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftForward(
                U[mu], mu,
                Gimpl::CovShiftBackward(
                    U[nu], nu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu, Gimpl::CovShiftBackward(U[mu], mu, U[nu])))),
            mu);

        //           __
        //          |__ __ |

        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftBackward(
                U[nu], nu,
                Gimpl::CovShiftBackward(
                    U[mu], mu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[mu])))),
            mu);

        //           __ ___
        //          |__    |

        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftForward(
                U[nu], nu,
                Gimpl::CovShiftBackward(
                    U[mu], mu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu, Gimpl::CovShiftBackward(U[nu], nu, U[mu])))),
            mu);

        //       --
        //      |  |
        //
        //      |  |

        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftForward(
                U[nu], nu,
                Gimpl::CovShiftForward(
                    U[nu], nu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu,
                        Gimpl::CovShiftBackward(
                            U[nu], nu,
                            Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
            mu);

        //      |  |
        //
        //      |  |
        //       --

        Stap += Gimpl::ShiftStaple(
            Gimpl::CovShiftBackward(
                U[nu], nu,
                Gimpl::CovShiftBackward(
                    U[nu], nu,
                    Gimpl::CovShiftBackward(
                        U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[nu])))),
            mu);
      }
    }
  }
};

typedef WilsonLoops<PeriodicGimplR> ColourWilsonLoops;
typedef WilsonLoops<PeriodicGimplR> U1WilsonLoops;
typedef WilsonLoops<PeriodicGimplR> SU2WilsonLoops;
typedef WilsonLoops<PeriodicGimplR> SU3WilsonLoops;
}
}

#endif
