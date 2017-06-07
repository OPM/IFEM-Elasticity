// $Id$
//==============================================================================
//!
//! \file ElasticCable.C
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nonlinear elastic cable.
//!
//==============================================================================

#include "ElasticCable.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"


void ElasticCable::printLog () const
{
  IFEM::cout <<"ElasticCable: Stiffness = "<< stiffness;
  if (EI > 0.0) IFEM::cout <<" "<< EI;
  IFEM::cout <<", Mass density = "<< lineMass << std::endl;
}


/*!
  \brief Evaluates the unit binormal and normal vectors at current point.
*/

static bool evalLocalAxes (const Vec3& dX, const Vec3& ddX,
                           Vec3& B, Vec3& N, double& B_len, double& N_len)
{
  B.cross(dX,ddX); // {B} = d{X} x dd{X}
  B_len = B.length();
  if (B_len < 1.0e-12)
  {
    // The cable is straigth, need to modify vector ddX

    Vec3 G1(dX);
    double normG1 = dX.length();
    G1 *= 1.0/(normG1*normG1); // Contra-variant basis vector

    B.cross(dX, ddX - (ddX*G1) * G1);
    B_len = B.length();
    if (B_len < 1.0e-10)
    {
      B.cross(dX, Vec3(-dX.y,dX.x,0.0));
      B_len = B.length();
      if (B_len < 1.0e-10)
      {
        // Last resort: Use algorithm from Tensor(const Vec3&) constructor
        if (fabs(dX.y) < fabs(dX.z))
        {
          // Define the normal vector, N, by projecting the global Y-axis
          // onto the normal plane of the tangent direction, dX
          N.x = -dX.y*dX.x;
          N.y =  dX.x*dX.x + dX.z*dX.z;
          N.z = -dX.y*dX.z;
          // Define the binormal vector B as the cross product of dX and N
          B.cross(dX,N);
        }
        else
        {
          // Define the binormal vector by projecting the global Z-axis
          // onto the normal plane of the tangent direction, dX
          B.x = -dX.z*dX.x;
          B.y = -dX.z*dX.y;
          B.z =  dX.x*dX.x + dX.y*dX.y;
        }
        B_len = B.length();
        if (B_len < 1.0e-10)
        {
          std::cerr <<" *** ElasticCable: Degenerated element, dX="<< dX
                    << std::endl;
          return false;
        }
      }
    }
  }

  B *= 1.0/B_len; // Unit binormal vector: {B}
  N.cross(B,dX);  // {N} = {B} x d{X}
  N_len = N.normalize(); // Unit normal vector: {N}

  return true;
}


bool ElasticCable::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const Vec3& X) const
{
  size_t a, aa, b, bb;
  unsigned char i, j, k, l, o;
  const size_t nen = fe.N.size();

  // Set up reference configuration

  Vec3 dX(fe.G.getColumn(1));
  Vec3 ddX(fe.G.getColumn(2));
#if INT_DEBUG > 1
  std::cout <<"ElasticCable: X = "<< X <<" dX = "<< dX <<" ddX = "<< ddX <<"\n";
#endif

  // Compute current configuration

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const Vector& eV = elMat.vec.front(); // Current displacement vector

  Vec3   x(X);
  Vec3  dx(dX);
  Vec3 ddx(ddX);
  for (i = 0; i < 3; i++)
  {
      x[i] += eV.dot(fe.N,i,3);
     dx[i] += eV.dot(fe.dNdX,i,3);
    ddx[i] += eV.dot(fe.d2NdX2,i,3);
  }
#if INT_DEBUG > 1
  std::cout <<"ElasticCable: x = "<< x <<" dx = "<< dx <<" ddx = "<< ddx <<"\n";
#endif

  // Compute local coordinate systems of the reference and current configuration

  Vec3 B_unit, N_unit;
  double B_len, N_len;
  if (!evalLocalAxes(dX,ddX,B_unit,N_unit,B_len,N_len)) return false;
#if INT_DEBUG > 1
  std::cout <<"ElasticCable: B_unit = "<< B_unit <<" N_unit = "<< N_unit <<"\n";
#endif

  Vec3 b_unit, n_unit;
  double b_len, n_len;
  if (!evalLocalAxes(dx,ddx,b_unit,n_unit,b_len,n_len)) return false;
  Vec3   bin    = b_unit * b_len;
  double b_len2 = b_len  * b_len;
  Vec3   n      = n_unit * n_len;
  double n_len2 = n_len  * n_len;
#if INT_DEBUG > 1
  std::cout <<"ElasticCable: b = "<< bin <<" b_unit = "<< b_unit
            <<"\n              n = "<< n <<" n_unit = "<< n_unit << std::endl;
#endif

  // Calculate derivative of b_unit

  std::vector<Matrix> db(nen,Matrix(3,3)), db_unit(nen,Matrix(3,3));
  std::vector<Vec3>   db_normal(nen);

  for (i = 1; i <= 3; i++)
    for (k = 1; k <= 3; k++)
      for (l = 1; l <= 3; l++)
      {
        double eps_kli = 0.5*(k-l)*(l-i)*(i-k);
        double eps_kil = 0.5*(k-i)*(i-l)*(l-k);
        for (a = 1; a <= nen; a++)
          db[a-1](k,i) += (eps_kil*fe.dNdX(a,1)*ddx[l-1] +
                           eps_kli*dx[l-1]*fe.d2NdX2(a,1,1));
      }

  for (i = 1; i <= 3; i++)
    for (a = 0; a < nen; a++)
      for (k = 1; k <= 3; k++)
        db_normal[a][i-1] += b_unit[k-1]*db[a](k,i);

  for (i = 1; i <= 3; i++)
    for (k = 1; k <= 3; k++)
      for (a = 0; a < nen; a++)
        db_unit[a](k,i) += (db[a](k,i) - b_unit[k-1]*db_normal[a][i-1])/b_len;

#if INT_DEBUG > 2
  std::cout <<"ElasticCable: db_unit:\n";
  for (a = 0; a < nen; a++)
    std::cout <<"node "<< a+1 << db_unit[a];
#endif

  // Calculate second derivative of b_unit

  std::vector< std::vector<Matrix3D> > ddb(nen), ddb_unit(nen);
  std::vector< std::vector<Matrix>   > ddb_normal(nen);
  for (a = 0; a < nen; a++)
  {
    ddb[a].resize(nen,Matrix3D(3,3,3));
    ddb_unit[a].resize(nen,Matrix3D(3,3,3));
    ddb_normal[a].resize(nen,Matrix(3,3));
  }

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (k = 1; k <= 3; k++)
      {
        double eps_kij = 0.5*(k-i)*(i-j)*(j-k);
        double eps_kji = 0.5*(k-j)*(j-i)*(i-k);
        for (a = 1; a <= nen; a++)
          for (b = 1; b <= nen; b++)
            ddb[a-1][b-1](k,i,j) = (eps_kji*fe.d2NdX2(a,1,1)*fe.dNdX(b,1) +
                                    eps_kij*fe.d2NdX2(b,1,1)*fe.dNdX(a,1));
      }

#if INT_DEBUG > 3
  std::cout <<"ElasticCable: ddb:\n";
  for (a = 0; a < nen; a++)
    for (b = 0; b < nen; b++)
      std::cout <<"nodes "<< a+1 <<","<< b+1 << ddb[a][b];
#endif

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
          for (k = 1; k <= 3; k++)
            ddb_normal[a][b](i,j) += (ddb[a][b](k,i,j)*bin[k-1] +
                                      db[a](k,i)*db[b](k,j) -
                                      bin[k-1]*db[a](k,i)*bin[k-1]*db[b](k,j) /
                                      b_len2) / b_len;

#if INT_DEBUG > 3
  std::cout <<"ElasticCable: ddb_normal:\n";
  for (a = 0; a < nen; a++)
    for (b = 0; b < nen; b++)
      std::cout <<"nodes "<< a+1 <<","<< b+1 << ddb_normal[a][b];
#endif

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
          for (k = 1; k <= 3; k++)
            ddb_unit[a][b](k,i,j) = (ddb[a][b](k,i,j)/b_len -
                                     db[a](k,i)*db_normal[b][j-1]/b_len2 -
                                     db[b](k,j)*db_normal[a][i-1]/b_len2 -
                                     bin[k-1]*(ddb_normal[a][b](i,j) -
                                               db_normal[a][i-1]*
                                               db_normal[b][j-1]*2.0 /
                                               b_len) / b_len2);

#if INT_DEBUG > 2
  std::cout <<"ElasticCable: ddb_unit:\n";
  for (a = 0; a < nen; a++)
    for (b = 0; b < nen; b++)
      std::cout <<"nodes "<< a+1 <<","<< b+1 << ddb_unit[a][b];
#endif

  // Calculate derivative of n_unit

  std::vector<Matrix> dn(nen,Matrix(3,3)), dn_unit(nen,Matrix(3,3));
  std::vector<Vec3>   dn_normal(nen);

  for (i = 1; i <= 3; i++)
    for (k = 1; k <= 3; k++)
      for (l = 1; l <= 3; l++)
      {
        double eps_kli = 0.5*(k-l)*(l-i)*(i-k);
        for (a = 0; a < nen; a++)
        {
          dn[a](k,i) += eps_kli*b_unit[l-1]*fe.dNdX(1+a,1);
          for (o = 1; o <= 3; o++)
          {
            double eps_kol = 0.5*(k-o)*(o-l)*(l-k);
            dn[a](k,i) += eps_kol*db_unit[a](o,i)*dx[l-1];
          }
        }
      }

  for (i = 1; i <= 3; i++)
    for (a = 0; a < nen; a++)
      for (k = 1; k <= 3; k++)
        dn_normal[a][i-1] += n_unit[k-1]*dn[a](k,i);

  for (i = 1; i <= 3; i++)
    for (k = 1; k <= 3; k++)
      for (a = 0; a < nen; a++)
        dn_unit[a](k,i) += (dn[a](k,i) - n_unit[k-1]*dn_normal[a][i-1])/n_len;

#if INT_DEBUG > 2
  std::cout <<"\nElasticCable: dn_unit:\n";
  for (a = 0; a < nen; a++)
    std::cout <<"node "<< a+1 << dn_unit[a];
#endif

  // Calculate second derivative of n_unit

  std::vector< std::vector<Matrix3D> > ddn(nen), ddn_unit(nen);
  std::vector< std::vector<Matrix>   > ddn_normal(nen);
  for (a = 0; a < nen; a++)
  {
    ddn[a].resize(nen,Matrix3D(3,3,3));
    ddn_unit[a].resize(nen,Matrix3D(3,3,3));
    ddn_normal[a].resize(nen,Matrix(3,3));
  }

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
          for (k = 1; k <= 3; k++)
            for (o = 1; o <= 3; o++)
            {
              double eps_koj = 0.5*(k-o)*(o-j)*(j-k);
              double eps_koi = 0.5*(k-o)*(o-i)*(i-k);
              ddn[a][b](k,i,j) += (eps_koj*db_unit[a](o,i)*fe.dNdX(1+b,1) +
                                   eps_koi*db_unit[b](o,j)*fe.dNdX(1+a,1));
              for (l = 1; l <= 3; l++)
              {
                double eps_kol = 0.5*(k-o)*(o-l)*(l-k);
                ddn[a][b](k,i,j) += eps_kol*ddb_unit[a][b](o,i,j)*dx[l-1];
              }
            }

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
          for (k = 1; k <= 3; k++)
            ddn_normal[a][b](i,j) += (ddn[a][b](k,i,j)*n[k-1] +
                                      dn[a](k,i)*dn[b](k,j) -
                                      n[k-1]*dn[a](k,i)*
                                      n[k-1]*dn[b](k,j)/n_len2) / n_len;

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
          for (k = 1; k <= 3; k++)
            ddn_unit[a][b](k,i,j) = (ddn[a][b](k,i,j)/n_len -
                                     dn[a](k,i)*dn_normal[b][j-1]/n_len2 -
                                     dn[b](k,j)*dn_normal[a][i-1]/n_len2 -
                                     n[k-1]*(ddn_normal[a][b](i,j) -
                                             dn_normal[a][i-1]*
                                             dn_normal[b][j-1]*2.0 /
                                             n_len) / n_len2);

#if INT_DEBUG > 2
  std::cout <<"ElasticCable: ddn_unit:\n";
  for (a = 0; a < nen; a++)
    for (b = 0; b < nen; b++)
      std::cout <<"nodes "<< a+1 <<","<< b+1 << ddn_unit[a][b];
#endif

  // Axial strain
  double eps = 0.5*(dx*dx - dX*dX);

  // Derivative of the axial strain
  Vector deps(3*nen);
  for (a = aa = 1; a <= nen; a++)
    for (i = 1; i <= 3; i++, aa++)
      deps(aa) = fe.dNdX(a,1)*dx[i-1];

  // Second derivative of the axial strain
  Matrix ddeps(3*nen,3*nen);
  for (a = 1; a <= nen; a++)
    for (b = 1; b <= nen; b++)
      for (i = 1; i <= 3; i++)
        ddeps(3*(a-1)+i,3*(b-1)+i) = fe.dNdX(a,1)*fe.dNdX(b,1);

  // Curvature
  double kappa = (ddx*n_unit - ddX*N_unit);

  // Derivative of the curvature
  Vector dkappa(3*nen);
  for (a = aa = 1; a <= nen; a++)
    for (i = 1; i <= 3; i++, aa++)
    {
      dkappa(aa) = fe.d2NdX2(a,1,1)*n_unit[i-1];
      for (k = 1; k <= 3; k++)
        dkappa(aa) += ddx[k-1]*dn_unit[a-1](k,i);
    }

  // Second derivative of the curvature
  Matrix ddkappa(3*nen,3*nen);
  for (a = 0, aa = 1; a < nen; a++)
    for (i = 1; i <= 3; i++, aa++)
      for (b = 0, bb = 1; b < nen; b++)
        for (j = 1; j <= 3; j++, bb++)
        {
          ddkappa(aa,bb) = (fe.d2NdX2(1+a,1,1)*dn_unit[b](i,j) +
                            fe.d2NdX2(1+b,1,1)*dn_unit[a](j,i));
          for (k = 1; k <= 3; k++)
            ddkappa(aa,bb) += ddx[k-1]*ddn_unit[a][b](k,i,j);
        }

#if INT_DEBUG > 1
  std::cout <<"ElasticCable: eps = "<< eps <<" kappa = "<< kappa
            <<"\ndeps:"<< deps <<"dkappa:"<< dkappa
            <<"ddeps:"<< ddeps <<"ddkappa:"<< ddkappa;
#endif

  // Norm of initial contravariant basis (G^1)
  double normG1contr2   = 1.0 / (dX.x*dX.x + dX.y*dX.y + dX.z*dX.z);
  double normG1contr4JW = normG1contr2 * normG1contr2 * fe.detJxW;

  double EAxJW = EA * normG1contr4JW; // volume-weighted axial stiffness
  double EIxJW = EI * normG1contr4JW; // volume-weighted bending stiffness

  if (iS)
  {
    // Integrate the internal forces (note the negative sign here)
    elMat.b[iS-1].add(deps,-eps*EAxJW);
    elMat.b[iS-1].add(dkappa,-kappa*EIxJW);
  }

  if (eKm)
  {
    // Integrate the material stiffness matrix
    elMat.A[eKm-1].outer_product(deps,deps*EAxJW,true);
    elMat.A[eKm-1].outer_product(dkappa,dkappa*EIxJW,true);
  }

  if (eKg)
  {
    // Integrate the geometric stiffness matrix
    elMat.A[eKg-1].add(ddeps,eps*EAxJW);
    elMat.A[eKg-1].add(ddkappa,kappa*EIxJW);
  }

  if (lineMass > 0.0)
  {
    double dMass = lineMass*fe.detJxW;
    if (eM)
    {
      // Integrate the mass matrix
      Matrix& M = elMat.A[eM-1];
      for (a = 1; a <= nen; a++)
        for (b = 1; b <= nen; b++)
          for (i = 1; i <= 3; i++)
            M(3*(a-1)+i,3*(b-1)+i) += fe.N(a)*fe.N(b)*dMass;
    }

    if (eS && !gravity.isZero())
    {
      // Integrate the external (gravitation) forces
      Vector& S = elMat.b[eS-1];
      for (a = 1; a <= nen; a++)
        for (i = 1; i <= 3; i++)
          S(3*(a-1)+i) += fe.N(a)*gravity[i-1]*dMass;
    }
  }

  return true;
}


bool ElasticCable::evalSol (Vector& s, const FiniteElement& fe, const Vec3& X,
                            const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
  {
    int ierr = utl::gather(MNPC,3,primsol.front(),eV);
    if (ierr > 0)
    {
      std::cerr <<" *** ElasticCable::evalSol: Detected "<< ierr
                <<" node numbers out of range."<< std::endl;
      return false;
    }
  }

  // Set up reference and current configuration

  Vec3  dX(fe.G.getColumn(1));
  Vec3 ddX(fe.G.getColumn(2));

  Vec3   x(X);
  Vec3  dx(dX);
  Vec3 ddx(ddX);
  for (size_t i = 0; i < 3; i++)
  {
      x[i] += eV.dot(fe.N,i,3);
     dx[i] += eV.dot(fe.dNdX,i,3);
    ddx[i] += eV.dot(fe.d2NdX2,i,3);
  }
#if INT_DEBUG > 1
  std::cout <<"ElasticCable: X = "<< X <<" u = "<< X-x <<"\n";
  std::cout <<"ElasticCable: x = "<< x <<" dx = "<< dx <<" ddx = "<< ddx <<"\n";
#endif

  // Compute local coordinate systems of the reference and current configuration

  Vec3 B_unit, N_unit;
  double B_len, N_len;
  if (!evalLocalAxes(dX,ddX,B_unit,N_unit,B_len,N_len))
    return false;

  Vec3 b_unit, n_unit;
  double b_len, n_len;
  if (!evalLocalAxes(dx,ddx,b_unit,n_unit,b_len,n_len))
    return false;

#if INT_DEBUG > 1
  std::cout <<"ElasticCable: N_unit = "<< N_unit
            <<" n_unit = "<< n_unit << std::endl;
#endif

  s.resize(2);
  s[0] = EA*0.5*(dx*dx - dX*dX);       // Axial force
  s[1] = EI*(ddx*n_unit - ddX*N_unit); // Bending moment
  return true;
}


std::string ElasticCable::getField2Name (size_t i, const char* prefix) const
{
  std::string name;
  if (prefix)
    name = std::string(prefix) + " ";
  else
    name.clear();

  switch (i) {
  case 0: name += "Axial force"; break;
  case 1: name += "Bending moment"; break;
  default: return "";
  }

  return name;
}
