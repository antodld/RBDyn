/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "RBDyn/Momentum.h"

// includes
// RBDyn
#include "RBDyn/Jacobian.h"
#include "RBDyn/MultiBody.h"
#include "RBDyn/MultiBodyConfig.h"

namespace rbd
{

sva::ForceVecd computeCentroidalMomentum(const MultiBody & mb, const MultiBodyConfig & mbc, const Eigen::Vector3d & com)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();
  Vector6d cm(Vector6d::Zero());

  sva::PTransformd X_com_0(Vector3d(-com));
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    // body inertia in body coordinate
    sva::ForceVecd hi = bodies[i].inertia() * mbc.bodyVelB[i];

    // momentum at CoM for link i : {}^iX_{com}^T {}^iI_i {}^iV_i
    cm += (mbc.bodyPosW[i] * X_com_0).transMul(hi).vector();
  }

  return sva::ForceVecd(cm);
}

sva::ForceVecd computeCentroidalMomentumDot(const MultiBody & mb,
                                            const MultiBodyConfig & mbc,
                                            const Eigen::Vector3d & com,
                                            const Eigen::Vector3d & comDot)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();
  Vector6d cm(Vector6d::Zero());

  sva::PTransformd X_0_com(com);
  sva::MotionVecd com_Vel(Vector3d::Zero(), comDot);

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    sva::MotionVecd body_i_Vel(mbc.bodyVelB[i]);
    sva::PTransformd X_com_i(mbc.bodyPosW[i] * X_0_com.inv());
    sva::PTransformd X_i_com(X_com_i.inv());

    sva::ForceVecd body_i_Momentum(bodies[i].inertia() * body_i_Vel);

    const auto I = bodies[i].inertia();

    // transform in com coordinate
    cm += X_i_com.dualMul( body_i_Vel.crossDual(I * body_i_Vel) + I * mbc.bodyAccB[i]).vector();
  }

  return sva::ForceVecd(cm);
}

sva::ForceVecd sComputeCentroidalMomentum(const MultiBody & mb,
                                          const MultiBodyConfig & mbc,
                                          const Eigen::Vector3d & com)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  return computeCentroidalMomentum(mb, mbc, com);
}

sva::ForceVecd sComputeCentroidalMomentumDot(const MultiBody & mb,
                                             const MultiBodyConfig & mbc,
                                             const Eigen::Vector3d & com,
                                             const Eigen::Vector3d & comDot)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodyAcc(mb, mbc);
  return computeCentroidalMomentumDot(mb, mbc, com, comDot);
}

Eigen::Matrix6d jacProjector(const sva::PTransformd & X_i_com, const sva::RBInertiad & I_i)
{
  return Eigen::Matrix6d(X_i_com.dualMatrix() * I_i.matrix());
}

Eigen::Matrix6d jacProjectorDot(const sva::PTransformd & X_i_com,
                                const sva::RBInertiad & I_i,
                                const sva::MotionVecd & V_i,
                                const sva::MotionVecd & V_com)
{
  Eigen::Matrix6d X_i_com_d = X_i_com.dualMatrix() * sva::vector6ToCrossDualMatrix(V_i.vector())
                              - sva::vector6ToCrossDualMatrix(V_com.vector()) * X_i_com.dualMatrix();
  return Eigen::Matrix6d(X_i_com_d * I_i.matrix());
}

CentroidalMomentumMatrix::CentroidalMomentumMatrix()
: cmMat_(), cmMatDot_(), jacVec_(), blocksVec_(), jacWork_(), bodiesWeight_()
{
}

CentroidalMomentumMatrix::CentroidalMomentumMatrix(const MultiBody & mb)
: cmMat_(6, mb.nrDof()), cmMatDot_(6, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
  blocksVec_(static_cast<size_t>(mb.nrBodies())), jacWork_(static_cast<size_t>(mb.nrBodies())),
  bodiesWeight_(static_cast<size_t>(mb.nrBodies()), 1.), normalAcc_(static_cast<size_t>(mb.nrBodies()))
{
  init(mb);
}

CentroidalMomentumMatrix::CentroidalMomentumMatrix(const MultiBody & mb, std::vector<double> weight)
: cmMat_(6, mb.nrDof()), cmMatDot_(6, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
  blocksVec_(static_cast<size_t>(mb.nrBodies())), jacWork_(static_cast<size_t>(mb.nrBodies())),
  bodiesWeight_(std::move(weight)), normalAcc_(static_cast<size_t>(mb.nrBodies()))
{
  init(mb);

  if(int(bodiesWeight_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "weight vector must be of size " << mb.nrBodies() << " not " << bodiesWeight_.size() << std::endl;
    throw std::domain_error(ss.str());
  }
}

void CentroidalMomentumMatrix::computeMatrix(const MultiBody & mb,
                                             const MultiBodyConfig & mbc,
                                             const Eigen::Vector3d & com)
{
  using namespace Eigen;
  const std::vector<Body> & bodies = mb.bodies();
  cmMat_.setZero();

  sva::PTransformd X_0_com(com);
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const MatrixXd & jac = jacVec_[i].bodyJacobian(mb, mbc);
    sva::PTransformd X_i_com(X_0_com * (mbc.bodyPosW[i].inv()));
    Eigen::MatrixXd j_full = Eigen::MatrixXd::Zero(cmMat_.rows(),cmMat_.cols());
    jacVec_[i].fullJacobian(mb,jac,j_full);
    const Matrix6d proj = X_i_com.dualMatrix() * bodies[i].inertia().matrix();

    cmMat_ += proj * j_full;
  }
}

void CentroidalMomentumMatrix::computeMatrixDot(const MultiBody & mb,
                                                const MultiBodyConfig & mbc,
                                                const Eigen::Vector3d & com,
                                                const Eigen::Vector3d & comDot)
{
  using namespace Eigen;
  const std::vector<Body> & bodies = mb.bodies();
  cmMatDot_.setZero();

  sva::PTransformd X_0_com(com);
  sva::MotionVecd com_Vel(Vector3d::Zero(), comDot);
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const MatrixXd & jac = jacVec_[i].bodyJacobian(mb, mbc);
    const MatrixXd & jacDot = jacVec_[i].bodyJacobianDot(mb, mbc);

    Eigen::MatrixXd j_full = Eigen::MatrixXd::Zero(6,mb.nrDof());
    Eigen::MatrixXd jdot_full = Eigen::MatrixXd::Zero(6,mb.nrDof());

    const sva::PTransformd X_i_com(X_0_com * (mbc.bodyPosW[i].inv()));
    const sva::MotionVecd V_i = mbc.bodyVelB[i];
    const auto I = bodies[i].inertia().matrix();

    const auto Idot = sva::vector6ToCrossDualMatrix(V_i.vector()) * I  ;

    jacVec_[i].fullJacobian(mb,jac,j_full);
    jacVec_[i].fullJacobian(mb,jacDot,jdot_full);

    cmMatDot_ +=  X_i_com.dualMatrix() * (Idot * j_full + I * jdot_full);

  }
}

void CentroidalMomentumMatrix::computeMatrixAndMatrixDot(const MultiBody & mb,
                                                         const MultiBodyConfig & mbc,
                                                         const Eigen::Vector3d & com,
                                                         const Eigen::Vector3d & comDot)
{

  computeMatrix(mb,mbc,com);
  computeMatrixDot(mb,mbc,com,comDot);

}

const Eigen::MatrixXd & CentroidalMomentumMatrix::matrix() const
{
  return cmMat_;
}

const Eigen::MatrixXd & CentroidalMomentumMatrix::matrixDot() const
{
  return cmMatDot_;
}

sva::ForceVecd CentroidalMomentumMatrix::momentum(const MultiBody & mb,
                                                  const MultiBodyConfig & mbc,
                                                  const Eigen::Vector3d & com) const
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();
  Vector6d cm(Vector6d::Zero());

  sva::PTransformd X_com_0(Vector3d(-com));
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    // body inertia in body coordinate
    sva::ForceVecd hi = bodies[i].inertia() * mbc.bodyVelB[i];

    // momentum at CoM for link i : {}^iX_{com}^T {}^iI_i {}^iV_i
    cm += ((mbc.bodyPosW[i] * X_com_0).transMul(hi).vector());
  }

  return sva::ForceVecd(cm);
}

sva::ForceVecd CentroidalMomentumMatrix::normalMomentumDot(const MultiBody & mb,
                                                           const MultiBodyConfig & mbc,
                                                           const Eigen::Vector3d & com,
                                                           const Eigen::Vector3d & comDot)
{
  using namespace Eigen;

  const std::vector<int> & pred = mb.predecessors();
  const std::vector<int> & succ = mb.successors();

  for(size_t i = 0; i < static_cast<size_t>(mb.nrJoints()); ++i)
  {
    const auto pred_index = static_cast<size_t>(pred[i]);
    const auto succ_index = static_cast<size_t>(succ[i]);

    const sva::PTransformd & X_p_i = mbc.parentToSon[i];
    const sva::MotionVecd & vj_i = mbc.jointVelocity[i];
    const sva::MotionVecd & vb_i = mbc.bodyVelB[i];

    if(pred[i] != -1)
      normalAcc_[succ_index] = X_p_i * normalAcc_[pred_index] + vb_i.cross(vj_i);
    else
      normalAcc_[succ_index] = vb_i.cross(vj_i);
  }

  return normalMomentumDot(mb, mbc, com, comDot, normalAcc_);
}

sva::ForceVecd CentroidalMomentumMatrix::normalMomentumDot(const MultiBody & mb,
                                                           const MultiBodyConfig & mbc,
                                                           const Eigen::Vector3d & com,
                                                           const Eigen::Vector3d & comDot,
                                                           const std::vector<sva::MotionVecd> & normalAccB) const
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();
  Vector6d cm(Vector6d::Zero());

  sva::PTransformd X_com_0(Vector3d(-com));
  sva::MotionVecd com_Vel(Vector3d::Zero(), comDot);

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    sva::MotionVecd body_i_Vel(mbc.bodyVelB[i]);
    sva::PTransformd X_com_i(mbc.bodyPosW[i] * X_com_0);

    sva::ForceVecd body_i_Momentum(bodies[i].inertia() * body_i_Vel);

    const auto I = bodies[i].inertia();

    // transform in com coordinate
    cm += X_com_i.transMul( I * normalAccB[i] + body_i_Vel.crossDual(I * body_i_Vel)).vector();
  }

  return sva::ForceVecd(cm);
}

void CentroidalMomentumMatrix::sComputeMatrix(const MultiBody & mb,
                                              const MultiBodyConfig & mbc,
                                              const Eigen::Vector3d & com)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);
  computeMatrix(mb, mbc, com);
}

void CentroidalMomentumMatrix::sComputeMatrixDot(const MultiBody & mb,
                                                 const MultiBodyConfig & mbc,
                                                 const Eigen::Vector3d & com,
                                                 const Eigen::Vector3d & comDot)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);
  computeMatrixDot(mb, mbc, com, comDot);
}

void CentroidalMomentumMatrix::sComputeMatrixAndMatrixDot(const MultiBody & mb,
                                                          const MultiBodyConfig & mbc,
                                                          const Eigen::Vector3d & com,
                                                          const Eigen::Vector3d & comDot)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);
  computeMatrixAndMatrixDot(mb, mbc, com, comDot);
}

sva::ForceVecd CentroidalMomentumMatrix::sMomentum(const MultiBody & mb,
                                                   const MultiBodyConfig & mbc,
                                                   const Eigen::Vector3d & com) const
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);

  return momentum(mb, mbc, com);
}

sva::ForceVecd CentroidalMomentumMatrix::sNormalMomentumDot(const MultiBody & mb,
                                                            const MultiBodyConfig & mbc,
                                                            const Eigen::Vector3d & com,
                                                            const Eigen::Vector3d & comDot)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchJointConf(mb, mbc);
  checkMatchParentToSon(mb, mbc);

  return normalMomentumDot(mb, mbc, com, comDot);
}

sva::ForceVecd CentroidalMomentumMatrix::sNormalMomentumDot(const MultiBody & mb,
                                                            const MultiBodyConfig & mbc,
                                                            const Eigen::Vector3d & com,
                                                            const Eigen::Vector3d & comDot,
                                                            const std::vector<sva::MotionVecd> & normalAccB) const
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodiesVector(mb, normalAccB, "normalAccB");

  return normalMomentumDot(mb, mbc, com, comDot, normalAccB);
}

void CentroidalMomentumMatrix::init(const rbd::MultiBody & mb)
{
  using namespace Eigen;
  for(int i = 0; i < mb.nrBodies(); ++i)
  {
    const auto ui = static_cast<size_t>(i);
    jacVec_[ui] = Jacobian(mb, mb.body(i).name());
    blocksVec_[ui] = jacVec_[ui].compactPath(mb);
    jacWork_[ui].resize(6, jacVec_[ui].dof());
  }
}

} // namespace rbd
