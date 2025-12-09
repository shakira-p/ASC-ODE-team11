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
  for (auto sp : mss.springs())et_initial_velocities({
    mA: (0, 0, 1)})
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

  virtual size_t dimX() const { return D*mss.masses().size() + mss.constraints().size(); }
  virtual size_t dimF() const { return dimX(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const
  {
    evaluateT(x, f);
  }

  // Helper to get position of a connector from input vector x
  template <typename T>
  Vec<D, T> getConnectorPos(VectorView<T> x, const Connector& c) const
  {
    Vec<D, T> pos;
    if (c.type == Connector::FIX) {
      for (int d = 0; d < D; d++)
        pos(d) = T(mss.fixes()[c.nr].pos(d));
    } else {
      // Read position from input vector x (first D*n_masses entries are positions)
      for (int d = 0; d < D; d++)
        pos(d) = x(c.nr * D + d);
    }
    return pos;
  }

  template <typename T>
  void evaluateT (VectorView<T> x, VectorView<T> f) const
  {
    size_t n_masses = mss.masses().size();
    size_t n_constraints = mss.constraints().size();

    for (size_t i = 0; i < f.size(); i++)
      f(i) = T(0.0);

    // Helper lambda to access f as a 2D array (row-major: f[i*D + d])
    auto fmat = [&](size_t i, int d) -> T& { return f(i * D + d); };

    // Gravity forces (converted to accelerations later if no constraints)
    for (size_t i = 0; i < n_masses; i++)
      for (int d = 0; d < D; d++)
        fmat(i, d) = T(mss.masses()[i].mass * mss.getGravity()(d));

    // Spring forces (elastic) - now using positions from input vector x
    for (auto spring : mss.springs())
    {
      auto [c1, c2] = spring.connectors;
      Vec<D, T> p1 = getConnectorPos(x, c1);
      Vec<D, T> p2 = getConnectorPos(x, c2);

      Vec<D, T> diff;
      for (int d = 0; d < D; d++)
        diff(d) = p1(d) - p2(d);

      T dist = vecNorm(diff);
      T force = T(spring.stiffness) * (dist - T(spring.length));
      Vec<D, T> dir;
      for (int d = 0; d < D; d++)
        dir(d) = (p2(d) - p1(d)) / dist;

      if (c1.type == Connector::MASS)
        for (int d = 0; d < D; d++)
          fmat(c1.nr, d) = fmat(c1.nr, d) + force * dir(d);
      if (c2.type == Connector::MASS)
        for (int d = 0; d < D; d++)
          fmat(c2.nr, d) = fmat(c2.nr, d) - force * dir(d);
    }

    if (n_constraints == 0)
    {
      // No constraints: just return accelerations (f = m*a -> a = f/m)
      for (size_t i = 0; i < n_masses; i++)
        for (int d = 0; d < D; d++)
          fmat(i, d) = fmat(i, d) / T(mss.masses()[i].mass);
    }
    else
    {
      // With constraints: this function is called with x = [positions, lambdas]
      // and should return f = [forces, constraint_values]
      // The time integrator composes this with the Newmark update

      // For constraint forces: add lambda * grad(g) to forces
      for (size_t c = 0; c < n_constraints; c++) {
        auto& con = mss.constraints()[c];
        auto [c1, c2] = con.connectors;
        T lambda = x(D * n_masses + c);

        Vec<D, T> p1 = getConnectorPos(x, c1);
        Vec<D, T> p2 = getConnectorPos(x, c2);

        Vec<D, T> diff;
        for (int d = 0; d < D; d++) diff(d) = p2(d) - p1(d);

        // grad(g) = 2*(p2-p1) for g = |p2-p1|^2 - L^2
        // Force on mass = -lambda * grad_mass(g)
        if (c1.type == Connector::MASS)
          for (int d = 0; d < D; d++)
            fmat(c1.nr, d) = fmat(c1.nr, d) + lambda * T(2.0) * diff(d);
        if (c2.type == Connector::MASS)
          for (int d = 0; d < D; d++)
            fmat(c2.nr, d) = fmat(c2.nr, d) - lambda * T(2.0) * diff(d);
      }

      // Convert forces to accelerations: a = F/m
      for (size_t i = 0; i < n_masses; i++) {
        T mass_val = T(mss.masses()[i].mass);
        for (int d = 0; d < D; d++)
          fmat(i, d) = fmat(i, d) / mass_val;
      }

      // Constraint equations: g(x) = |p2-p1|^2 - L^2 = 0
      for (size_t c = 0; c < n_constraints; c++) {
        auto& con = mss.constraints()[c];
        auto [c1, c2] = con.connectors;

        Vec<D, T> p1 = getConnectorPos(x, c1);
        Vec<D, T> p2 = getConnectorPos(x, c2);

        T dist_sq = T(0.0);
        for (int d = 0; d < D; d++)
          dist_sq = dist_sq + (p2(d) - p1(d)) * (p2(d) - p1(d));

        f(D * n_masses + c) = dist_sq - T(con.length * con.length);
      }
    }
  }

  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    const size_t N = dimX();

    // Use dynamic AutoDiff - size determined at runtime
    Vector<AutoDiff<>> xad(N);
    for (size_t i = 0; i < N; i++)
    {
      xad(i) = AutoDiff<>(x(i), i, N);  // value, derivIndex, size
    }

    Vector<AutoDiff<>> fad(dimF());
    VectorView<AutoDiff<>> xad_view(N, xad.data());
    VectorView<AutoDiff<>> fad_view(dimF(), fad.data());

    evaluateT(xad_view, fad_view);

    for (size_t i = 0; i < dimF(); i++)
      for (size_t j = 0; j < N; j++)
        df(i, j) = derivative(fad(i), j);

    //// Numerical differentiation
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
