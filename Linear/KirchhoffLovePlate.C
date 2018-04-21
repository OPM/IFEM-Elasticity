// $Id$
//==============================================================================
//!
//! \file KirchhoffLovePlate.C
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for linear Kirchhoff-Love thin plate problems.
//!
//==============================================================================

#include "KirchhoffLovePlate.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "IFEM.h"


KirchhoffLovePlate::KirchhoffLovePlate (unsigned short int n, short int v)
  : KirchhoffLove(n)
{
  version = v;
  nOrder = 0;
}


void KirchhoffLovePlate::printLog () const
{
  IFEM::cout <<"KirchhoffLovePlate: thickness = "<< thickness
             <<", gravity = "<< gravity << std::endl;
  if (version > 1)
    IFEM::cout <<"\tUsing tensorial formulation (constant D)."<< std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLovePlate*>(this)->material = &defaultMat;
  }

  material->printLog();
  if (version > 1)
    IFEM::cout <<"\tPlate stiffness: D = "
               << this->getStiffness(Vec3()) << std::endl;
}


double KirchhoffLovePlate::getStiffness (const Vec3& X) const
{
  return material->getPlateStiffness(X,thickness);
}


/*!
  The strain-displacement matrix for a Kirchhoff-Love plate element is formally
  defined as:
  \f[
  [B] = \left[\begin{array}{c}
   \frac{\partial^2}{\partial x^2} \\
   \frac{\partial^2}{\partial y^2} \\
  2\frac{\partial^2}{\partial x\partial y}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [NENOD] row-vector.
*/

bool KirchhoffLovePlate::formBmatrix (Matrix& Bmat,
				      const Matrix3D& d2NdX2) const
{
  const size_t nenod = d2NdX2.dim(1);
  const size_t nstrc = nsd*(nsd+1)/2;

  Bmat.resize(nstrc,nenod,true);

  if (d2NdX2.dim(2) != nsd || d2NdX2.dim(3) != nsd)
  {
    std::cerr <<" *** KirchhoffLovePlate::formBmatrix: Invalid dimension on"
	      <<" d2NdX2, "<< d2NdX2.dim(1) <<"x"<< d2NdX2.dim(2)
	      <<"x"<< d2NdX2.dim(3) <<"."<< std::endl;
    return false;
  }

  for (size_t i = 1; i <= nenod; i++)
    if (nsd == 1)
      Bmat(1,i) = d2NdX2(i,1,1);
    else
    {
      Bmat(1,i) = d2NdX2(i,1,1);
      Bmat(2,i) = d2NdX2(i,2,2);
      Bmat(3,i) = d2NdX2(i,1,2)*2.0;
    }

  return true;
}


bool KirchhoffLovePlate::formCmatrix (Matrix& C, const FiniteElement& fe,
                                      const Vec3& X, bool invers) const
{
  SymmTensor dummy(nsd); double U;
  if (!material->evaluate(C,dummy,U,fe,X,dummy,dummy, invers ? -1 : 1))
    return false;

  double factor = thickness*thickness*thickness/12.0;
  C.multiply(invers ? 1.0/factor : factor);
  return true;
}


bool KirchhoffLovePlate::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eS) // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,fe.detJxW);

  if (!eK)
    return true;

  // Integrate the stiffness matrix
  if (version == 1)
    return this->evalK1(elMat.A[eK-1],fe,X);
  else
    return this->evalK2(elMat.A[eK-1],fe,X);
}


bool KirchhoffLovePlate::evalK1 (Matrix& EK,
                                 const FiniteElement& fe, const Vec3& X) const
{
  // Compute the strain-displacement matrix B from d2NdX2
  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.d2NdX2))
    return false;

  // Evaluate the constitutive matrix at this point
  Matrix Cmat;
  if (!this->formCmatrix(Cmat,fe,X))
    return false;

  // Integrate the stiffness matrix
  Matrix CB;
  CB.multiply(Cmat,Bmat).multiply(fe.detJxW); // CB = C*B*|J|*w
  EK.multiply(Bmat,CB,true,false,true); // EK += B^T * CB

  return true;
}


bool KirchhoffLovePlate::evalK2 (Matrix& EK,
                                 const FiniteElement& fe, const Vec3& X) const
{
  // Evaluate the scaled plate stiffness at this point
  double DdJxW = material->getPlateStiffness(X,thickness)*fe.detJxW;

  // Integrate the stiffness matrix
  Vector d2NdX2 = fe.d2NdX2.getColumn(1,1);
  if (!EK.outer_product(d2NdX2,d2NdX2*DdJxW,true)) // EK += N,xx*N,xx^t*|J|*w
    return false;
  else if (nsd < 2)
    return true;

  Vector d2NdY2 = fe.d2NdX2.getColumn(2,2);
  EK.outer_product(d2NdY2,d2NdY2*DdJxW,true); // EK += N,yy*N,yy^t*|J|*w
  if (version == 2)
  {
    Vector d2NdXY = fe.d2NdX2.getColumn(1,2);
    EK.outer_product(d2NdXY,d2NdXY*DdJxW*2.0,true); // EK += 2*N,xy*N,xy*|J|*w
  }
  else
  {
    EK.outer_product(d2NdX2,d2NdY2*DdJxW,true); // EK += N,xx*N,yy^t*|J|*w
    EK.outer_product(d2NdY2,d2NdX2*DdJxW,true); // EK += N,yy*N,xx^t*|J|*w
  }

  return true;
}


bool KirchhoffLovePlate::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& normal) const
{
  std::cerr <<" *** KirchhoffLovePlate::evalBou not implemented."<< std::endl;
  return false;
}


bool KirchhoffLovePlate::evalSol (Vector& s, const Vector& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  if (eV.empty())
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: No displacement vector."
	      << std::endl;
    return false;
  }
  else if (eV.size() != fe.d2NdX2.dim(1))
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: Invalid displacement vector."
              <<"\n     size(eV) = "<< eV.size() <<"   size(d2NdX2) = "
              << fe.d2NdX2.dim(1) <<","<< fe.d2NdX2.dim(2)*fe.d2NdX2.dim(3)
              << std::endl;
    return false;
  }

  if (version == 1)
  {
    // Compute the strain-displacement matrix B from d2NdX2
    Matrix Bmat;
    if (!this->formBmatrix(Bmat,fe.d2NdX2))
      return false;

    // Evaluate the constitutive matrix at this point
    Matrix Cmat;
    if (!this->formCmatrix(Cmat,fe,X))
      return false;

    // Evaluate the curvature tensor
    SymmTensor kappa(nsd), m(nsd);
    if (!Bmat.multiply(eV,kappa)) // kappa = B*eV
      return false;

    // Evaluate the stress resultant tensor
    if (!Cmat.multiply(-1.0*kappa,m)) // m = -C*kappa
      return false;

    // Congruence transformation to local coordinate system at current point
    if (toLocal && locSys) m.transform(locSys->getTmat(X));

    s = m;
  }
  else
  {
    // Compute the Laplacian components, w,xx w,yy and w,xy
    s.resize(this->getNoFields(2));
    for (size_t i = 1; i <= fe.d2NdX2.dim(2); i++)
      s(i) = eV.dot(fe.d2NdX2.getColumn(i,i));
    if (s.size() > 2)
      s(3) = eV.dot(fe.d2NdX2.getColumn(1,2));
  }

#if INT_DEBUG > 3
  std::cout <<"KirchhoffLovePlate::evalSol("<< fe.iel <<", "<< X <<"):\n\ts =";
  for (double v : s) std::cout <<" "<< v;
  std::cout << std::endl;
#endif
  return true;
}


size_t KirchhoffLovePlate::getNoFields (int fld) const
{
  return fld < 2 ? 1 : (version > 2 ? nsd : nsd*(nsd+1)/2);
}


std::string KirchhoffLovePlate::getField1Name (size_t,
					       const char* prefix) const
{
  if (!prefix) return "w";

  return prefix + std::string(" w");
}


std::string KirchhoffLovePlate::getField2Name (size_t i,
					       const char* prefix) const
{
  if (i >= 3) return "";

  static const char* s1[3] = { "m_xx", "m_yy", "m_xy" };
  static const char* s2[3] = { "d2w/dx2", "d2w/dy2", "d2w/dxy" };
  const char** s = version == 1 ? s1 : s2;

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* KirchhoffLovePlate::getNormIntegrand (AnaSol* asol) const
{
  if (!asol)
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this));
  else if (version == 1)
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this),
                                      asol->getStressSol());
  else
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this),
                                      asol->getScalarSecSol());
}


KirchhoffLovePlateNorm::KirchhoffLovePlateNorm (KirchhoffLovePlate& p,
                                                STensorFunc* a)
  : NormBase(p), anasol(a), ana2nd(nullptr)
{
  nrcmp = myProblem.getNoFields(2);
  projBou = true;
}


KirchhoffLovePlateNorm::KirchhoffLovePlateNorm (KirchhoffLovePlate& p,
                                                VecFunc* a)
  : NormBase(p), anasol(nullptr), ana2nd(a)
{
  nrcmp = myProblem.getNoFields(2);
  projBou = true;
}


bool KirchhoffLovePlateNorm::evalInt (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X) const
{
  KirchhoffLovePlate& problem = static_cast<KirchhoffLovePlate&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  int version = problem.getVersion();

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (version == 1 && !problem.formCmatrix(Cinv,fe,X,true))
    return false;

  // Evaluate the finite element stress field
  Vector mh, m, error;
  if (!problem.evalSol(mh,pnorm.vec.front(),fe,X))
    return false;
  else if (mh.size() == 3 && version > 1)
    mh.push_back(mh(3));

  double hk4 = fe.h*fe.h*fe.h*fe.h;
  double Res = 0.0;
  size_t ip = 0;

  // Integrate the energy norm a(w^h,w^h)
  if (version == 1)
    pnorm[ip++] += mh.dot(Cinv*mh)*fe.detJxW;
  else
    pnorm[ip++] += mh.dot(mh)*fe.detJxW;

  // Evaluate the body load
  double p = problem.getPressure(X);
  if (version > 1) p /= problem.getStiffness(X);
  // Evaluate the displacement field
  double w = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (p,w^h)
  pnorm[ip++] += p*w*fe.detJxW;
  if (version > 1) p = -p;

#if INT_DEBUG > 3
  std::cout <<"KirchhoffLovePlateNorm::evalInt("<< fe.iel <<", "<< X <<"):";
#endif

  if (version == 1 && anasol)
  {
    // Evaluate the analytical stress resultant field
    m = (*anasol)(X);

    // Integrate the energy norm a(w,w)
    pnorm[ip++] += m.dot(Cinv*m)*fe.detJxW;
    // Integrate the error in energy norm a(w-w^h,w-w^h)
    error = m - mh;
    pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

    // Residual of analytical solution (should be zero)
    if (nrcmp == 1)
      Res = anasol->dderiv(X,1,1)(1,1) + p;
    else
    {
      double d2MxxdX2 = anasol->dderiv(X,1,1)(1,1);
      double d2MyydY2 = anasol->dderiv(X,2,2)(2,2);
      double d2MxydXY = anasol->dderiv(X,1,2)(1,2);
      Res = d2MxxdX2 + d2MxydXY + d2MxydXY + d2MyydY2 + p;
    }
    // Integrate the residual error in the analytical solution
    pnorm[ip++] += hk4*Res*Res*fe.detJxW;
  }
  else if (version > 1 && ana2nd)
  {
    // Evaluate the analytical Laplacian
    m = (*ana2nd)(X).vec(nrcmp);
    if (nrcmp == 3) m.push_back(m(3));

    // Integrate the energy norm a(w,w)
    pnorm[ip++] += m.dot(m)*fe.detJxW;
    // Integrate the error in energy norm a(w-w^h,w-w^h)
    error = m - mh;
    pnorm[ip++] += error.dot(error)*fe.detJxW;

    // Residual of analytical solution (should be zero)
    if (nrcmp == 1)
      Res = ana2nd->dderiv(X,1,1).x + p;
    else
    {
      double wxxxx = ana2nd->dderiv(X,1,1).x;
      double wxxyy = ana2nd->dderiv(X,2,2).x;
      double wyyyy = ana2nd->dderiv(X,2,2).y;
      Res = wxxxx + wxxyy + wxxyy + wyyyy + p;
    }
    // Integrate the residual error in the analytical solution
    pnorm[ip++] += hk4*Res*Res*fe.detJxW;
  }

  // Integrate the area
  pnorm[ip++] += fe.detJxW;

  size_t nen = fe.N.size();
  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate the projected solution
      Vector mr(nrcmp);
      for (unsigned short int j = 0; j < nrcmp; j++)
        mr[j] = psol.dot(fe.N,j,nrcmp);
      if (nrcmp == 3 && version > 1)
        mr.push_back(mr(3));

#if INT_DEBUG > 3
      std::cout <<"\n\t"<< (version > 1 ? "Laplace{w^r} =" : "m^r =");
      for (double  v : mr) std::cout <<" "<< v;
#endif

      // Integrate the energy norm a(w^r,w^r)
      if (version == 1)
        pnorm[ip++] += mr.dot(Cinv*mr)*fe.detJxW;
      else
        pnorm[ip++] += mr.dot(mr)*fe.detJxW;

      // Integrate the error in energy norm a(w^r-w^h,w^r-w^h)
      error = mr - mh;
      if (version == 1)
        pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
      else
        pnorm[ip++] += error.dot(error)*fe.detJxW;

      // Integrate the L2-norm (m^r,m^r)
      pnorm[ip++] += mr.dot(mr.ptr(),nrcmp)*fe.detJxW;
      // Integrate the error in L2-norm (m^r-m^h,m^r-m^h)
      pnorm[ip++] += error.dot(error.ptr(),nrcmp)*fe.detJxW;

      // Evaluate the interior residual of the projected solution
      if (nrcmp == 1)
      {
        Res = psol.dot(fe.d2NdX2) + p;
#if INT_DEBUG > 3
        if (version > 1) std::cout <<"\n\tw,xxxx^r = "<< Res-p;
#endif
      }
      else if (version == 1)
      {
        double d2mxxdx2 = psol.dot(fe.d2NdX2,0,3);
        double d2myydy2 = psol.dot(fe.d2NdX2,1,3,nen*3);
        double d2mxydxy = psol.dot(fe.d2NdX2,2,3,nen);
        Res = d2mxxdx2 + d2mxydxy + d2mxydxy + d2myydy2 + p;
#if INT_DEBUG > 3
        std::cout <<"\n\tmxx,xx^r = "<< d2mxxdx2
                  <<"\n\tmyy,yy^r = "<< d2myydy2
                  <<"\n\tmxy,xy^r = "<< d2mxydxy;
#endif
      }
      else
      {
        double wxxxx = psol.dot(fe.d2NdX2,0,nrcmp);
        double wyyyy = psol.dot(fe.d2NdX2,1,nrcmp,nen*3);
        double wxxyy = psol.dot(fe.d2NdX2,0,nrcmp,nen*3);
#if INT_DEBUG > 3
        std::cout <<"\n\tw,xxxx^r = "<< wxxxx
                  <<"\n\tw,yyyy^r = "<< wyyyy
                  <<"\n\tw,xxyy^r = "<< wxxyy;
#endif
        Res = wxxxx + wxxyy + wxxyy + wyyyy + p;
      }
      // Integrate the residual error in the projected solution
#if INT_DEBUG > 3
      std::cout <<"\n\tResidual h^4*|Laplace{m^r}-p|^2 = "
                << hk4*Res*Res << std::endl;
#endif
      pnorm[ip++] += hk4*Res*Res*fe.detJxW;
      ip++; // Make room for Jump contributions here

      if (version == 1 && anasol)
      {
        // Integrate the error in the projected solution a(w-w^r,w-w^r)
        error = m - mr;
        pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
        ip += 2; // Make room for the local effectivity indices here
      }
      else if (version > 1 && ana2nd)
      {
        // Integrate the error in the projected solution a(w-w^r,w-w^r)
        error = m - mr;
        pnorm[ip++] += error.dot(error)*fe.detJxW;
        ip += 2; // Make room for the local effectivity indices here
      }
    }

  if (ip == pnorm.size())
    return true;

  std::cerr <<" *** KirchhoffLovePlateNorm::evalInt: Internal error, ip="
            << ip <<" != pnorm.size()="<< pnorm.size() << std::endl;
  return false;
}


bool KirchhoffLovePlateNorm::evalBou (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X, const Vec3& normal) const
{
  if (nrcmp <= 1) return true; // Nothing for 1D problems (beams)

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  int version = static_cast<KirchhoffLovePlate&>(myProblem).getVersion();

  double Jmp, m1, m2, hk3 = fe.h*fe.h*fe.h;

  size_t ip = 2;
  if (version == 1 && anasol)
  {
    ip += 3;
    Vector mx, my;
    if (nOrder%2)
    {
      // Evaluate the analytical moment field
      mx = (*anasol)(X);

      // Jump in the analytical moment (should be zero)
      m1 = mx[0]*normal.x + mx[2]*normal.y;
      m2 = mx[2]*normal.x + mx[1]*normal.y;
      Jmp = m1*m1 + m2*m2;

      // Integrate the edge jump in the analytical solution
      pnorm[ip-1] += fe.h*Jmp*fe.detJxW;
    }
    if (nOrder/2)
    {
      // Evaluate the analytical shear force field, q = {n}*grad{m}
      mx = anasol->deriv(X,1);
      my = anasol->deriv(X,2);
      Jmp = (mx[0]+my[2])*normal.x + (mx[2]+my[1])*normal.y;

      // Integrate the edge jump in the analytical solution
      pnorm[ip-1] += hk3*Jmp*Jmp*fe.detJxW;
    }
  }
  else if (version > 1 && ana2nd)
  {
    ip += 3;
    if (nOrder%2)
    {
      // Evaluate the analytical Laplacian components
      Vec3 wxx = (*ana2nd)(X);
      Jmp = wxx.x*wxx.x + wxx.y*wxx.y;

      // Integrate the edge jump in the analytical solution
      pnorm[ip-1] += fe.h*Jmp*fe.detJxW;
    }
    if (nOrder/2)
    {
      // Evaluate the gradient of the analytical Laplacian components
      Vec3 wxxx = ana2nd->deriv(X,1);
      Vec3 wxxy = ana2nd->deriv(X,2);
      Jmp = (wxxx.x+wxxy.x)*normal.x + (wxxx.y+wxxy.y)*normal.y;

      // Integrate the edge jump in the analytical solution
      pnorm[ip-1] += hk3*Jmp*Jmp*fe.detJxW;
    }
  }

  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      ip += 6;
      if (nOrder%2)
      {
        // Evaluate the projected solution
        Vector mr(nrcmp);
        for (unsigned short int j = 0; j < nrcmp; j++)
          mr[j] = psol.dot(fe.N,j,nrcmp);

        // Jump in the projected solution
        if (version == 1)
        {
          m1 = mr[0]*normal.x + mr[2]*normal.y;
          m2 = mr[2]*normal.x + mr[1]*normal.y;
          Jmp = m1*m1 + m2*m2;
        }
        else
          Jmp = mr[0]*normal.x + mr[1]*normal.y;

        // Integrate the edge jump in the recovered solution
        pnorm[ip-1] += fe.h*Jmp*fe.detJxW;
        pnorm[ip]   += fe.h*Jmp*fe.detJxW;
      }
      if (nOrder/2)
      {
        // Evaluate the gradient of the projected solution.
        // Notice that the matrix multiplication method used here treats
        // the element vector, psol, as a matrix whose number
        // of columns equals the number of rows in the matrix fe.dNdX.
        Matrix dmdX;
        if (!dmdX.multiplyMat(psol,fe.dNdX)) // dmdX = psol*dNdX
          return false;

        if (version == 1) // Shear force q^r = {n}*grad{m^r}
          Jmp = (dmdX(1,1)+dmdX(3,2))*normal.x + (dmdX(3,1)+dmdX(2,2))*normal.y;
        else // {n}*grad{Laplace{w^r}}
          Jmp = (dmdX(1,1)+dmdX(2,1))*normal.x + (dmdX(1,2)+dmdX(2,2))*normal.y;

        // Integrate the residual error in the analytical solution
        pnorm[ip-1] += hk3*Jmp*Jmp*fe.detJxW;
        pnorm[ip]   += hk3*Jmp*Jmp*fe.detJxW;
      }
      if (anasol) ip += 3;
    }

  return true;
}


bool KirchhoffLovePlateNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol && !ana2nd)
    return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as a(e^r,e^r)/a(e,e)
  // with e^r = w^r - w^h  and  e = w - w^h,
  // and (a(e^r,e^r)+res(w^r))/a(e,e)
  for (size_t ip = 14; ip < pnorm.size(); ip += 9)
  {
    pnorm[ip-1] = pnorm[ip-7] / pnorm[3];
    pnorm[ip] = (pnorm[ip-7]+pnorm[ip-4]) / pnorm[3];
  }

  return true;
}


int KirchhoffLovePlateNorm::getIntegrandType () const
{
  return SECOND_DERIVATIVES | ELEMENT_CORNERS;
}


size_t KirchhoffLovePlateNorm::getNoFields (int group) const
{
  if (group == 0)
    return this->NormBase::getNoFields();
  else if (group == 1 || group == -1)
    return anasol || ana2nd ? 6 : 3;
  else if (group > 0 || !prjsol[-group-2].empty())
    return anasol || ana2nd ? 9 : 6;
  else
    return 0;
}


std::string KirchhoffLovePlateNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  if (i == 0 || j == 0 || j > 9 || (i == 1 && j > 6))
    return this->NormBase::getName(i,j,prefix);

  static const char* u[6] = {
    "a(w^h,w^h)^0.5",
    "(p,w^h)^0.5",
    "a(w,w)^0.5",
    "a(e,e)^0.5, e=w-w^h",
    "res(w)^0.5",
    "area"
  };

  static const char* p[9] = {
    "a(w^r,w^r)^0.5",
    "a(e,e)^0.5, e=w^r-w^h",
    "(w^r,w^r)^0.5",
    "(e,e)^0.5, e=w^r-w^h",
    "res(w^r)^0.5",
    "jump(w^r)^0.5",
    "a(e,e)^0.5, e=w-w^r",
    "effectivity index^*",
    "effectivity index^RES"
  };

  if (!anasol && !ana2nd && i == 1 && j == 3) j = 6;
  std::string name(i > 1 ? p[j-1] : u[j-1]);

  return prefix ? prefix + std::string(" ") + name : name;
}
