// $Id$
//==============================================================================
//!
//! \file AnalyticSolutions.C
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for linear elasticity problems.
//!
//==============================================================================

#include "AnalyticSolutions.h"
#include "Vec3.h"


/*!
  \class Hole

  Smooth plane strain problem.

  Reference: Zienkiewicz & Zhu, "The superconvergent patch recovery...",
  IJNME, 33, 1331-1364, 1992, page 1355 (eq. 36).
*/

SymmTensor Hole::evaluate (const Vec3& X) const
{
  double R  = hypot(X.x,X.y);
  double th = atan2(X.y,X.x);
  double C2 = cos(2.0*th);
  double C4 = cos(4.0*th);
  double S2 = sin(2.0*th);
  double S4 = sin(4.0*th);
  double R2 = R <= a ? 1.0 : a*a/(R*R);
  double R4 = R <= a ? 1.0 : R2*R2;

  SymmTensor sigma(is3D ? 3 : 2, true);
  sigma(1,1) = F0 * (1.0 - R2*(1.5*C2 + C4) + 1.5*R4*C4);
  sigma(2,2) = F0 * (    - R2*(0.5*C2 - C4) - 1.5*R4*C4);
  sigma(1,2) = F0 * (    - R2*(0.5*S2 + S4) + 1.5*R4*S4);
  sigma(3,3) = F0 * nu*(1.0 - 2.0*R2*C2);

  return sigma;
}


/*!
  \class Lshape

  Symmetric (mode 1) loading. Plane strain problem with a singular point.

  Reference: Szabo & Babuska, "Finite Element Analysis", 1991, page 188-192.

  \code
                   y ^
                     |

     (-a,a) +--------+(0,a)
            |        |
            |        |
            |        |
            |        +--------+ (a,0) - - -> x
            |   (0,0)         |
            |                 |
            |                 |
    (-a,-a) +-----------------+ (a,-a)

   \endcode
*/

Lshape::Lshape (double r, double f, double P, bool use3D)
  : STensorFunc(use3D ? 3 : 2, true), a(r), F0(f), nu(P), is3D(use3D), T(2)
{
  // Set up the local-to-global transformation tensor
  T(1,1) = T(2,2) = T(2,1) = -sqrt(0.5);
  T(1,2) = sqrt(0.5);
}


SymmTensor Lshape::evaluate (const Vec3& X) const
{
  // Some constants (see Szabo & Babuska for an elaboration)
  const double lambda = 0.544483737;
  const double q      = 0.543075579;
  const double lp1    = lambda + 1.0;
  const double lm1    = lambda - 1.0;
  const double lm3    = lambda - 3.0;
  const double tol    = a*1.0e-32;

  // Find local (polar) coordinates
  double x = X.x*T(1,1) + X.y*T(2,1);
  double y = X.x*T(1,2) + X.y*T(2,2);
  double r = hypot(x,y);
  if (r < tol) r = tol; // truncate the singularity to avoid NaN values
  double theta = atan2(y,x);

  // Set up the stress tensor in local system
  SymmTensor sigma(is3D ? 3 : 2, true);
  double c0  = F0*lambda*pow(r,lm1);
  sigma(1,1) = c0*((2.0-q*lp1)*cos(lm1*theta) - lm1*cos(lm3*theta));
  sigma(2,2) = c0*((2.0+q*lp1)*cos(lm1*theta) + lm1*cos(lm3*theta));
  sigma(1,2) = c0*(     q*lp1 *sin(lm1*theta) + lm1*sin(lm3*theta));
  sigma(3,3) = nu * (sigma(1,1)+sigma(2,2));

  // Transform to global coordinates
  return sigma.transform(T);
}


SymmTensor CanTS::evaluate (const Vec3& X) const
{
  double x = X.x/L;
  double y = (is3D ? X.z : X.y)/H - 0.5;
  double I = H*H*H / 12.0;
  size_t n = is3D ? 3 : 2;

  SymmTensor sigma(n);
  sigma(1,1) = F0*L*H/I * (x-1.0)*y;
  sigma(1,n) = F0*H*H/I * 0.5*(0.25-y*y);

  return sigma;
}


SymmTensor CanTM::evaluate (const Vec3& X) const
{
  double y = (is3D ? X.z : X.y)/H - 0.5;
  double I = H*H*H / 12.0;

  SymmTensor sigma(is3D ? 3 : 2);
  sigma(1,1) = M0*H/I * y;

  return sigma;
}


/*!
  \class CurvedBeam

  Smooth plane stress problem with a prescribed edge displacement.

  Reference: Zienkiewicz & Taylor, "The Finite Element Method",
  pages 42-45.
*/

CurvedBeam::CurvedBeam (double u0, double Ri, double Ro, double E, bool use3D)
  : STensorFunc(use3D ? 3 : 2), a(Ri), b(Ro), is3D(use3D)
{
  PN = -u0*E/(M_PI*(a*a+b*b));
}


SymmTensor CurvedBeam::evaluate (const Vec3& X) const
{
  // Find polar coordinates
  double r     = hypot(X.x,X.y);
  double theta = atan2(X.y,X.x);

  // Set up the stress tensor in polar coordinates
  SymmTensor sigma(is3D ? 3 : 2);
  double c1 = a*a*b*b/(r*r*r);
  double c2 = (a*a + b*b)/r;
  double ct = cos(theta);
  double st = sin(theta);
  sigma(1,1) = PN*(    r + c1 - c2)*st;
  sigma(2,2) = PN*(3.0*r - c1 - c2)*st;
  sigma(1,2) = PN*(   -r - c1 + c2)*ct;

  // Transform to Cartesian coordinates
  Tensor T(2);
  T(1,1) =  ct;
  T(2,1) =  st;
  T(1,2) = -st;
  T(2,2) =  ct;

  return sigma.transform(T);
}


/*!
  \class Pipe

  Expansion of a circular cylinder subjected to a temperature gradient
  over the thickness.

  Reference: SINTEF Report 2014: NEST Project. Isogeometric Finite Element
  Solver For Thermal Stress Analysis. page 11.
*/

Pipe::Pipe (double Ri, double Ro, double Ti, double To, double T0,
            double E, double ny, double alpha, bool use3D, bool usePolar)
  : STensorFunc(use3D ? 3 : 2, true), is3D(use3D), polar(usePolar),
    Tin(Ti), Tex(To), ra(Ri), rb(Ro), nu(ny)
{
  double r = rb/ra;
  ln_rb_ra = log(r);
  rba2m1   = r*r - 1.0;
  C        = E*alpha*(Tex-Tin)*0.5/(1.0-nu*nu);
  // Note: In the reference it is used (1-nu) in the denominator. I suspect
  // that is a mis-print since we get matching results only when using (1-nu^2).
  Ea    = E*alpha;
  T_ref = T0;
}


SymmTensor Pipe::evaluate (const Vec3& X) const
{
  double r     = hypot(X.x,X.y);
  double ln_rb = log(rb/r);
  double rbr2  = (rb/r)*(rb/r);
  double Temp  = Tin + (Tex-Tin)*log(r/ra)/ln_rb_ra;

  // Stress tensor in polar coordinates
  SymmTensor sigma(is3D ? 3 : 2, true);
  sigma(1,1) = C*( ln_rb     /ln_rb_ra - (rbr2-1.0)/rba2m1);
  sigma(2,2) = C*((ln_rb-1.0)/ln_rb_ra + (rbr2+1.0)/rba2m1);
  sigma(3,3) = nu * (sigma(1,1)+sigma(2,2)) - Ea*(Temp-T_ref);
  if (polar) return sigma;

  // Local-to-global transformation
  Tensor T(2);
  T(1,1) =  X.x/r;
  T(2,1) =  X.y/r;
  T(1,2) = -T(2,1);
  T(2,2) =  T(1,1);

  // Transform to global Cartesian coordinates
  return sigma.transform(T);
}


//! \brief Enum defining thin-plate solution derivatives.
enum Derivative { W = 0, dWdx = 1, dWdy = 2,
                  d2Wdx2 = 11, d2Wdxdy = 12, d2Wdydx = 21, d2Wdy2 = 22 };


ThinPlateSol::ThinPlateSol (double E, double v, double t) : nu(v)
{
  D = E*t*t*t / (12.0 - 12.0*nu*nu);
}


/*!
  \class NavierPlate

  Navier plate solution for a constant pressure load or
  a point load at an arbitrary position.

  Reference: S. E. Weberg, "Todimensjonal elastisitetsteori", 1975, pages 61-67.
*/

NavierPlate::NavierPlate (double a, double b, double t, double E, double Poiss,
                          double P, int max_mn)
  : ThinPlateSol(E,Poiss,t), STensorFunc(2),
    pz(P), type(0), xi(0.0), eta(0.0), c2(0.0), d2(0.0), mxmn(max_mn), inc(2)
{
  alpha = M_PI/a;
  beta  = M_PI/b;

  scalSol.push_back(new Displ(pz/D,alpha,beta,xi,eta,c2,d2,type,max_mn,inc));
  stressSol = this;

  // Calculate and print the maximum displacement (at the centre x=a/2, y=b/2)
  std::streamsize oldPrec = std::cout.precision(10);
  std::cout <<"\n\nNavierPlate: w_max = "
            << (*scalSol.front())(Vec3(0.5*a,0.5*b,0.0))
            <<"\n             Max. number of terms in Fourier series = "
            << max_mn << std::endl;
  std::cout.precision(oldPrec);
}


NavierPlate::NavierPlate (double a, double b, double t, double E, double Poiss,
                          double P, double xi_, double eta_,
                          double c, double d, int max_mn)
  : ThinPlateSol(E,Poiss,t), STensorFunc(2),
    pz(P), type(2), mxmn(max_mn), inc(1)
{
  alpha = M_PI/a;
  beta  = M_PI/b;
  xi    = xi_*a;
  eta   = eta_*b;
  if (c == 0.0 || d == 0.0) type = 1;
  if (xi_ == 0.5 && eta_ == 0.5) inc = 2;
  c2    = type == 1 ? a : 0.5*c;
  d2    = type == 1 ? b : 0.5*d;

  scalSol.push_back(new Displ(pz/D,alpha,beta,xi,eta,c2,d2,type,max_mn,inc));
  stressSol = this;

  // Calculate and print the displacement at the centre x=a/2, y=b/2
  std::streamsize oldPrec = std::cout.precision(10);
  std::cout <<"\n\nNavierPlate: w_centre = "
            << (*scalSol.front())(Vec3(0.5*a,0.5*b,0.0))
            <<"\n             Max. number of terms in Fourier series = "
            << max_mn << std::endl;
  std::cout.precision(oldPrec);
}


double NavierPlate::Displ::evaluate (const Vec3& X) const
{
  double w = 0.0;
  for (int m = 1; m <= mn; m += inc)
    for (int n = 1; n <= mn; n += inc)
    {
      double am   = alpha*m;
      double bn   = beta*n;
      double abmn = am*am + bn*bn;
      if (type == 0)
	w += sin(am*X.x)*sin(bn*X.y) / (double(m)*double(n)*abmn*abmn);
      else
      {
	double pzmn = sin(am*xi)*sin(bn*eta) / (abmn*abmn);
	if (type == 2)
	  pzmn *= sin(am*c2)*sin(bn*d2) / (double(m)*double(n));
	w += pzmn*sin(am*X.x)*sin(bn*X.y);
      }
    }

  if (type == 1)
    w *=  4.0*pzD / (c2*d2);
  else
    w *= 16.0*pzD / (M_PI*M_PI);

  return w;
}


void NavierPlate::addTerms (std::vector<double>& M,
                            double x, double y, int m, int n, int deriv) const
{
  double am  = alpha*m;
  double bn  = beta*n;
  double am2 = am*am;
  double bn2 = bn*bn;
  double dmn = double(m) * double(n);

  double pzmn = 0.0;
  switch (type) {
  case 0: // uniform pressure
    pzmn = 1.0 / dmn;
    break;
  case 1: // concentrated point load
    pzmn = sin(am*xi)*sin(bn*eta);
    break;
  case 2: // partial load
    pzmn = sin(am*xi)*sin(bn*eta) * sin(am*c2)*sin(bn*d2) / dmn;
    break;
  }

  pzmn /= (am2+bn2)*(am2+bn2);
  switch (deriv) {
  case W:
    M[0] += pzmn*(am2 + nu*bn2)*sin(am*x)*sin(bn*y);
    M[1] += pzmn*(bn2 + nu*am2)*sin(am*x)*sin(bn*y);
    M[2] += pzmn*(nu - 1.0)*am*bn*cos(am*x)*cos(bn*y);
    break;
  case dWdx:
    M[0] += pzmn*(am2 + nu*bn2)*am*cos(am*x)*sin(bn*y);
    M[1] += pzmn*(bn2 + nu*am2)*am*cos(am*x)*sin(bn*y);
    M[2] += pzmn*(1.0 - nu)*am2*bn*sin(am*x)*cos(bn*y);
    break;
  case dWdy:
    M[0] += pzmn*(am2 + nu*bn2)*bn*sin(am*x)*cos(bn*y);
    M[1] += pzmn*(bn2 + nu*am2)*bn*sin(am*x)*cos(bn*y);
    M[2] += pzmn*(1.0 - nu)*am*bn2*cos(am*x)*sin(bn*y);
    break;
  case d2Wdx2:
    M[0] -= pzmn*(am2 + nu*bn2)*am2*sin(am*x)*sin(bn*y);
    M[1] -= pzmn*(bn2 + nu*am2)*am2*sin(am*x)*sin(bn*y);
    M[2] += pzmn*(1.0 - nu)*am2*am*bn*cos(am*x)*cos(bn*y);
    break;
  case d2Wdy2:
    M[0] -= pzmn*(am2 + nu*bn2)*bn2*sin(am*x)*sin(bn*y);
    M[1] -= pzmn*(bn2 + nu*am2)*bn2*sin(am*x)*sin(bn*y);
    M[2] += pzmn*(1.0 - nu)*am*bn2*bn*cos(am*x)*cos(bn*y);
    break;
  case d2Wdxdy:
  case d2Wdydx:
    M[0] += pzmn*(am2 + nu*bn2)*am*bn*cos(am*x)*cos(bn*y);
    M[1] += pzmn*(bn2 + nu*am2)*am*bn*cos(am*x)*cos(bn*y);
    M[2] += pzmn*(nu - 1.0)*am2*bn2*sin(am*x)*sin(bn*y);
    break;
  }
}


SymmTensor NavierPlate::evaluate (const Vec3& X, int deriv) const
{
  SymmTensor M(2);

#ifdef REVERSED_SUMMATION
  const maxMN = mxmn - (inc-1);
  for (int i = maxMN; i > 0; i -= inc)
  {
    this->addTerms(M,X.x,X.y,i,i,deriv);
    for (int j = i-inc; j > 0; j -= inc)
    {
      this->addTerms(M,X.x,X.y,i,j,deriv);
      this->addTerms(M,X.x,X.y,j,i,deriv);
    }
  }
#else
  const double eps = 1.0e-8;
  double prev = 0.0;
  for (int i = 1; i <= mxmn; i += inc)
  {
    for (int j = 1; j < i; j += inc)
    {
      this->addTerms(M,X.x,X.y,i,j,deriv);
      this->addTerms(M,X.x,X.y,j,i,deriv);
    }
    this->addTerms(M,X.x,X.y,i,i,deriv);

    double norm = M.L2norm();
#if SP_DEBUG > 3
    if (i == 1) std::cout <<"\nNavierPlate, X = "<< X.x <<" "<< X.y <<"\n";
    std::cout << i <<": "<< M(1,1) <<" "<< M(2,2) <<" "<< M(1,2)
              <<" -> "<< norm <<" "<< fabs(norm-prev)/norm << std::endl;
#endif
    if (i%2)
    {
      if (fabs(norm-prev) < eps*norm)
        break;
      else
        prev = norm;
    }
  }
#endif

  if (type == 1) // concentrated load
    M *= 4.0*pz * (alpha/M_PI)*(beta/M_PI);
  else // uniform or partial pressure load
    M *= 16.0*pz / (M_PI*M_PI);

  return M;
}


SymmTensor NavierPlate::evaluate (const Vec3& X) const
{
  return this->evaluate(X,0);
}


SymmTensor NavierPlate::deriv (const Vec3& X, int dir) const
{
  return this->evaluate(X,dir);
}


SymmTensor NavierPlate::dderiv (const Vec3& X, int dir1, int dir2) const
{
  return this->evaluate(X,10*dir1+dir2);
}


CircularPlate::CircularPlate (double r, double t, double E, double n, double P)
  : ThinPlateSol(E,n,t), STensorFunc(2), R(r)
{
  M0 = 0.25*P/M_PI;
  scalSol.push_back(new Displ(P,R,D,nu));
  stressSol = this;
}


CircularPlate::Displ::Displ (double P, double r, double D, double nu) : R(r)
{
  U0 = -P*R*R/(16.0*M_PI*D);
  C0 = (3.0+nu)/(1.0+nu);
}


double CircularPlate::Displ::evaluate (const Vec3& X) const
{
  double rR = hypot(X.x,X.y)/R;
  if (rR < 0.0001) rR = 0.0001;

  double r2 = rR*rR;
  return U0*(C0*(1.0-r2) - 2.0*r2*log10(rR));
}


SymmTensor CircularPlate::evaluate (const Vec3& X) const
{
  double r = hypot(X.x,X.y);
  if (r < 0.0001*R) r = 0.0001*R;

  // Bending moments in polar coordinates
  SymmTensor M(2);
  M(1,1) = M0*(1+nu)*log10(r/R);
  M(2,2) = M(1,1) + M0*(1.0-nu);

  // Local-to-global transformation
  Tensor T(2);
  T(1,1) =  X.x/r;
  T(2,1) =  X.y/r;
  T(1,2) = -T(2,1);
  T(2,2) =  T(1,1);

  // Transform to global Cartesian coordinates
  return M.transform(T);
}
