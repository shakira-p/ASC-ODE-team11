#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <autodiff.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};

class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

class DistanceConstraint
{
public:
  double length;
  std::array<Connector,2> connectors;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  Vec<D> m_gravity=0.0;
  std::vector<DistanceConstraint> m_constraints;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  size_t addConstraint(DistanceConstraint c)
  {
    m_constraints.push_back(c);
    return m_constraints.size()-1;
  }

  auto & constraints() { return m_constraints; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const { return D*mss.masses().size(); }
  virtual size_t dimF() const { return D*mss.masses().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    evaluateT(x, f);
  }

  template <typename T>
  void evaluateT (VectorView<T> x, VectorView<T> f) const
  {
    f = 0.0;

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto fmat = f.asMatrix(mss.masses().size(), D);

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D, T> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        T dist = norm(p1 - p2);
        T force = spring.stiffness * (dist - spring.length);
        Vec<D, T> dir12 = T(1.0)/ dist * (p2-p1);

        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) *= T(1.0)/mss.masses()[i].mass;
  }

  // template <typename T>
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    // TODO: exact differentiation
    const size_t N = dimX();

    Vector<AutoDiff<N>> xad(N);
    for (size_t i=0; i<N; i++)
    {
      xad(i) = AutoDiff<N>(x(i));  // no {} because one-liner
      xad(i).deriv()[i] = 1.0;  // wtf really setting manually?
    }

    Vector<AutoDiff<N>> fad(dimF());
    evaluateT(xad.view(), fad.view());

    for (size_t i = 0; i < dimF(); i++)
      for (size_t j = 0; j < N; j++)
        df(i,j) = derivative(fad(i), j);

    //////
    // double eps = 1e-8;
    // Vector<> xl(dimX()), xr(dimX()), fl(dimF()), fr(dimF());
    // for (size_t i = 0; i < dimX(); i++)
    //   {
    //     xl = x;
    //     xl(i) -= eps;
    //     xr = x;
    //     xr(i) += eps;
    //     evaluateT (xl, fl);
    //     evaluateT (xr, fr);
    //     df.col(i) = 1/(2*eps) * (fr-fl);
    //   }
  }
  
};

#endif
