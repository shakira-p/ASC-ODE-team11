#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ImprovedEuler : public TimeStepper
  {
    Vector<> m_vecf;
    Vector<> m_y_hat;  // adding variable storage
  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
    : TimeStepper(rhs), m_vecf(rhs->dimF()), m_y_hat(rhs->dimX()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      //evaluating f at current y
      this->m_rhs->evaluate(y, m_vecf);

      // intermediate point
      m_y_hat = y;
      m_y_hat += tau/2 * m_vecf;

      // evaluating at y_hat
      this->m_rhs->evaluate(m_y_hat, m_vecf);

      // actual update: y += tau * f(y_hat)
      y += tau * m_vecf;
    }
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };


class CrankNicolson : public TimeStepper
{
  std::shared_ptr<NonlinearFunction> m_equ;
  std::shared_ptr<Parameter> m_tau;
  std::shared_ptr<ConstantFunction> m_yold;
  std::shared_ptr<ConstantFunction> m_fold;
  Vector<> m_vecf_old;
public:
  CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
  : TimeStepper(rhs),
    m_tau(std::make_shared<Parameter>(0.0)),
    m_yold(std::make_shared<ConstantFunction>(rhs->dimX())),
    m_fold(std::make_shared<ConstantFunction>(rhs->dimF())),
    m_vecf_old(rhs->dimF())
  {
    // m_equ will be constructed in DoStep,  cause tau and f_old are first known there
  }

  void DoStep(double tau, VectorView<double> y) override
  {
    // save old value y_old = y
    m_yold->set(y);

    // f_old = f(y_old)
    this->m_rhs->evaluate(y, m_vecf_old);
    m_fold->set(m_vecf_old);

    //  R(y_new) = y_new - y_old - (tau/2)*(f_old + f(y_new))
    auto ynew = std::make_shared<IdentityFunction>(this->m_rhs->dimX());
    auto tau_param = std::make_shared<Parameter>(0.5 * tau);
    m_equ = ynew - m_yold - tau_param * (m_fold + this->m_rhs);

    // Newton solves R(y_new)=0, start value is current y
    NewtonSolver(m_equ, y);
  }
};






  

}


#endif
