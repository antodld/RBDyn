/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "RBDyn/CoM.h"

// RBDyn
#include "RBDyn/MultiBody.h"
#include "RBDyn/MultiBodyConfig.h"
#include "RBDyn/Momentum.h"

namespace rbd
{

Eigen::Vector3d computeCoM(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();

  Vector3d com = Vector3d::Zero();
  double totalMass = 0.;

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    double mass = bodies[i].inertia().mass();

    totalMass += mass;
    sva::PTransformd scaledBobyPosW(mbc.bodyPosW[i].rotation(), mass * mbc.bodyPosW[i].translation());
    com += (sva::PTransformd(bodies[i].inertia().momentum()) * scaledBobyPosW).translation();
  }

  assert(totalMass > 0 && "Invalid multibody. Totalmass must be strictly positive");
  return com / totalMass;
}

sva::MotionVecd computeCoMVelocity6D(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const Vector3d com = computeCoM(mb,mbc);
  const sva::ForceVecd hc_ = rbd::computeCentroidalMomentum(mb,mbc,com);

  Matrix6d Ic = centroidalInertia(mb,mbc);

  return sva::MotionVecd(Ic.inverse() * hc_.vector());
}

Eigen::Matrix6d centroidalInertia(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;
  const Vector3d com = computeCoM(mb,mbc);
  Matrix6d Ic = Matrix6d::Zero();
  const sva::PTransformd X_0_c = sva::PTransformd(com);

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const MatrixXd I = mb.bodies()[i].inertia().matrix();
    const sva::PTransformd X_0_b = mbc.bodyPosW[i];
    const sva::PTransformd X_b_c = X_0_c * X_0_b.inv();

    Ic += (X_b_c.dualMatrix() * I * X_b_c.inv().matrix()) ;

  }

  return Ic;
}

Eigen::Vector3d computeCoMVelocity(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  return computeCoMVelocity6D(mb,mbc).linear();
}

Eigen::Vector3d computeCoMAcceleration(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  return computeCoMAcceleration6D(mb,mbc).linear();
}

sva::MotionVecd computeCoMAcceleration6D(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const Vector3d com = computeCoM(mb,mbc);
  const sva::MotionVecd comVel = computeCoMVelocity6D(mb,mbc);
  const auto hc_ = rbd::computeCentroidalMomentum(mb,mbc,com);
  const auto hcVel_ = rbd::computeCentroidalMomentumDot(mb,mbc,com,comVel.linear());

  const auto Ic = centroidalInertia(mb,mbc);

  return sva::MotionVecd(Ic.inverse() * (hcVel_ - comVel.crossDual(hc_)).vector());
}

Eigen::Vector3d sComputeCoM(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  return computeCoM(mb, mbc);
}

Eigen::Vector3d sComputeCoMVelocity(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  return computeCoMVelocity(mb, mbc);
}

Eigen::Vector3d sComputeCoMAcceleration(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodyAcc(mb, mbc);
  return computeCoMAcceleration(mb, mbc);
}

/**
 *														CoMJacobianDummy
 */

CoMJacobianDummy::CoMJacobianDummy() {}

CoMJacobianDummy::CoMJacobianDummy(const MultiBody & mb)
: jac_(6, mb.nrDof()), jacDot_(6, mb.nrDof()), jacFull_(6, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
  totalMass_(0.), bodiesWeight_(static_cast<size_t>(mb.nrBodies()), 1.)
{
  init(mb);
}

CoMJacobianDummy::CoMJacobianDummy(const MultiBody & mb, std::vector<double> weight)
: jac_(6, mb.nrDof()), jacDot_(6, mb.nrDof()), jacFull_(6, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
  totalMass_(0.), bodiesWeight_(std::move(weight))
{
  init(mb);

  if(int(bodiesWeight_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "weight vector must be of size " << mb.nrBodies() << " not " << bodiesWeight_.size() << std::endl;
    throw std::domain_error(ss.str());
  }
}

CoMJacobianDummy::~CoMJacobianDummy() {}

Eigen::MatrixXd CoMJacobianDummy::jacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  rbd::CentroidalMomentumMatrix centroidalMat = rbd::CentroidalMomentumMatrix(mb);
  const Vector3d com = computeCoM(mb,mbc);
  centroidalMat.computeMatrix(mb,mbc,com);

  const MatrixXd A = centroidalMat.matrix();
  const MatrixXd Ic = centroidalInertia(mb,mbc);
  jac_ = (Ic.inverse() * A);

  return jac_.block(3,0,3,mb.nrDof());
}

Eigen::MatrixXd CoMJacobianDummy::jacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  jacobian(mb,mbc);
  rbd::CentroidalMomentumMatrix centroidalMat = rbd::CentroidalMomentumMatrix(mb);
  const Vector3d com = computeCoM(mb,mbc);
  const sva::MotionVecd comVel = computeCoMVelocity6D(mb,mbc);

  centroidalMat.computeMatrixDot(mb,mbc,com,comVel.linear());

  const MatrixXd Adot = centroidalMat.matrixDot();

  const MatrixXd Ic = centroidalInertia(mb,mbc);
  jacDot_ = Ic.inverse() * (Adot - sva::vector6ToCrossDualMatrix<double>(comVel.vector()) * Ic * jac_);

  return jacDot_.block(3,0,3,mb.nrDof());
}

Eigen::MatrixXd CoMJacobianDummy::sJacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobian(mb, mbc);
}

Eigen::MatrixXd CoMJacobianDummy::sJacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobianDot(mb, mbc);
}

void CoMJacobianDummy::init(const rbd::MultiBody & mb)
{
  using namespace Eigen;
  for(int i = 0; i < mb.nrBodies(); ++i)
  {
    double bodyMass = mb.body(i).inertia().mass();
    Vector3d comT(0, 0, 0);
    if(bodyMass > 0) comT = mb.body(i).inertia().momentum() / bodyMass;

    jacVec_[static_cast<size_t>(i)] = Jacobian(mb, mb.body(i).name(), comT);
    totalMass_ += mb.body(i).inertia().mass();
  }
}

/**
 *														CoMJacobian
 */

CoMJacobian::CoMJacobian() {}

CoMJacobian::CoMJacobian(const MultiBody & mb)
: jac_(6, mb.nrDof()), jacDot_(6, mb.nrDof()), bodiesCoeff_(static_cast<size_t>(mb.nrBodies())),
  bodiesCoM_(static_cast<size_t>(mb.nrBodies())), jointsSubBodies_(static_cast<size_t>(mb.nrJoints())),
  bodiesCoMWorld_(static_cast<size_t>(mb.nrBodies())), bodiesCoMVelB_(static_cast<size_t>(mb.nrBodies())),
  normalAcc_(static_cast<size_t>(mb.nrJoints())), weight_(static_cast<size_t>(mb.nrBodies()), 1.)
{
  init(mb);
}

CoMJacobian::CoMJacobian(const MultiBody & mb, std::vector<double> weight)
: jac_(6, mb.nrDof()), jacDot_(6, mb.nrDof()), bodiesCoeff_(static_cast<size_t>(mb.nrBodies())),
  bodiesCoM_(static_cast<size_t>(mb.nrBodies())), jointsSubBodies_(static_cast<size_t>(mb.nrJoints())),
  bodiesCoMWorld_(static_cast<size_t>(mb.nrBodies())), bodiesCoMVelB_(static_cast<size_t>(mb.nrBodies())),
  normalAcc_(static_cast<size_t>(mb.nrJoints())), weight_(std::move(weight))
{
  if(int(weight_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "weight vector must be of size " << mb.nrBodies() << " not " << weight_.size() << std::endl;
    throw std::domain_error(ss.str());
  }

  init(mb);
}

void CoMJacobian::updateInertialParameters(const MultiBody & mb)
{
  double mass = 0.;

  for(int i = 0; i < mb.nrBodies(); ++i)
  {
    mass += mb.body(i).inertia().mass();
  }

  for(int i = 0; i < mb.nrBodies(); ++i)
  {
    double bodyMass = mb.body(i).inertia().mass();
    bodiesCoeff_[static_cast<size_t>(i)] = (bodyMass * weight_[static_cast<size_t>(i)]) / mass;
    if(bodyMass > 0)
      bodiesCoM_[static_cast<size_t>(i)] = sva::PTransformd((mb.body(i).inertia().momentum() / bodyMass).eval());
    else
      bodiesCoM_[static_cast<size_t>(i)] = sva::PTransformd::Identity();
  }
}

const std::vector<double> & CoMJacobian::weight() const
{
  return weight_;
}

void CoMJacobian::weight(const MultiBody & mb, std::vector<double> w)
{
  weight_ = std::move(w);
  updateInertialParameters(mb);
}

Eigen::MatrixXd CoMJacobian::jacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;
  rbd::CentroidalMomentumMatrix centroidalMat = rbd::CentroidalMomentumMatrix(mb);
  const Vector3d com = computeCoM(mb,mbc);
  centroidalMat.computeMatrix(mb,mbc,com);

  const MatrixXd A = centroidalMat.matrix();
  const MatrixXd Ic = centroidalInertia(mb,mbc);
  jac_ = Ic.inverse() * A;
  return jac_.block(3,0,3,mb.nrDof());
}

Eigen::MatrixXd CoMJacobian::jacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  jacDot_.setZero();
  jacobian(mb,mbc);
  rbd::CentroidalMomentumMatrix centroidalMat = rbd::CentroidalMomentumMatrix(mb);
  const Vector3d com = computeCoM(mb,mbc);
  const sva::MotionVecd comVel = computeCoMVelocity6D(mb,mbc);

  centroidalMat.computeMatrixDot(mb,mbc,com,comVel.linear());

  const MatrixXd Adot = centroidalMat.matrixDot();

  const MatrixXd Ic = centroidalInertia(mb,mbc);
  const MatrixXd Jc = jacobian(mb,mbc);
  jacDot_ = Ic.inverse() * (Adot - sva::vector6ToCrossDualMatrix<double>(comVel.vector()) * Ic * jac_);

  return jacDot_.block(3,0,3,mb.nrDof());
}

Eigen::Vector3d CoMJacobian::velocity(const MultiBody & mb, const MultiBodyConfig & mbc) const
{
  return computeCoMVelocity(mb,mbc);
}

Eigen::Vector3d CoMJacobian::normalAcceleration(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  const std::vector<int> & pred = mb.predecessors();
  const std::vector<int> & succ = mb.successors();

  for(size_t i = 0; i < static_cast<size_t>(mb.nrJoints()); ++i)
  {
    const sva::PTransformd & X_p_i = mbc.parentToSon[i];
    const sva::MotionVecd & vj_i = mbc.jointVelocity[i];
    const sva::MotionVecd & vb_i = mbc.bodyVelB[i];

    const auto pred_index = static_cast<size_t>(pred[i]);
    const auto succ_index = static_cast<size_t>(succ[i]);

    if(pred[i] != -1)
      normalAcc_[succ_index] = X_p_i * normalAcc_[pred_index] + vb_i.cross(vj_i);
    else
      normalAcc_[succ_index] = vb_i.cross(vj_i);
  }

  return normalAcceleration(mb, mbc, normalAcc_);
}

Eigen::Vector3d CoMJacobian::normalAcceleration(const MultiBody & mb,
                                                const MultiBodyConfig & mbc,
                                                const std::vector<sva::MotionVecd> & normalAccB) const
{
  Eigen::Vector3d comA(Eigen::Vector3d::Zero());
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const Eigen::Vector3d & comT = bodiesCoM_[i].translation();

    // Normal Acceleration at CoM : com_T_b·A_b
    // Normal Acceleration at CoM world frame :
    //    0_R_b·com_T_b·A_b + 0_R_b_d·com_T_b·V_b
    // O_R_b_d : (Angvel_W)_b x 0_R_b
    sva::PTransformd X_0_i(mbc.bodyPosW[i].rotation().transpose(), comT);
    sva::MotionVecd angvel_W(mbc.bodyVelW[i].angular(), Eigen::Vector3d::Zero());
    comA += (X_0_i * normalAccB[i]).linear() * bodiesCoeff_[i];
    comA += (angvel_W.cross(X_0_i * mbc.bodyVelB[i])).linear() * bodiesCoeff_[i];
  }

  return comA;
}

void CoMJacobian::sUpdateInertialParameters(const MultiBody & mb)
{
  if(int(bodiesCoeff_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "mb should have " << bodiesCoeff_.size() << " bodies, not " << mb.nrBodies() << std::endl;
    throw std::domain_error(ss.str());
  }

  updateInertialParameters(mb);
}

void CoMJacobian::sWeight(const MultiBody & mb, std::vector<double> w)
{
  if(int(bodiesCoeff_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "mb should have " << bodiesCoeff_.size() << " bodies, not " << mb.nrBodies() << std::endl;
    throw std::domain_error(ss.str());
  }

  if(int(weight_.size()) != mb.nrBodies())
  {
    std::stringstream ss;
    ss << "weight vector must be of size " << mb.nrBodies() << " not " << weight_.size() << std::endl;
    throw std::domain_error(ss.str());
  }

  weight(mb, w);
}

Eigen::MatrixXd CoMJacobian::sJacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobian(mb, mbc);
}

Eigen::MatrixXd CoMJacobian::sJacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobianDot(mb, mbc);
}

Eigen::Vector3d CoMJacobian::sVelocity(const MultiBody & mb, const MultiBodyConfig & mbc) const
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);

  return velocity(mb, mbc);
}

Eigen::Vector3d CoMJacobian::sNormalAcceleration(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchJointConf(mb, mbc);
  checkMatchParentToSon(mb, mbc);

  return normalAcceleration(mb, mbc);
}

Eigen::Vector3d CoMJacobian::sNormalAcceleration(const MultiBody & mb,
                                                 const MultiBodyConfig & mbc,
                                                 const std::vector<sva::MotionVecd> & normalAccB) const
{
  checkMatchBodyPos(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodiesVector(mb, normalAccB, "normalAccB");

  return normalAcceleration(mb, mbc, normalAccB);
}

// inefficient but the best we can do without mbg
void jointBodiesSuccessors(const MultiBody & mb, int joint, std::vector<int> & subBodies)
{
  int sonBody = mb.successor(joint);
  subBodies.push_back(sonBody);
  for(int i = 0; i < mb.nrJoints(); ++i)
  {
    if(mb.predecessor(i) == sonBody)
    {
      jointBodiesSuccessors(mb, i, subBodies);
    }
  }
}

void CoMJacobian::init(const MultiBody & mb)
{
  updateInertialParameters(mb);

  for(int i = 0; i < mb.nrJoints(); ++i)
  {
    std::vector<int> & subBodies = jointsSubBodies_[static_cast<size_t>(i)];
    jointBodiesSuccessors(mb, i, subBodies);
  }
}

} // namespace rbd
