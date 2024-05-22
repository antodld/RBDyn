/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "RBDyn/CoM.h"

// RBDyn
#include "RBDyn/Momentum.h"
#include "RBDyn/MultiBody.h"
#include "RBDyn/MultiBodyConfig.h"
#include "RBDyn/NumericalIntegration.h"

namespace rbd
{

void resetCoMFrame(const MultiBody & mb, MultiBodyConfig & mbc)
{
  const Eigen::Vector3d com = computeCoM(mb, mbc);
  auto quat = Eigen::Quaterniond(mbc.bodyPosW[0].rotation().transpose());
  mbc.com_ori = {quat.w(), quat.x(), quat.y(), quat.z()};
  mbc.com = sva::PTransformd(QuatToE<double>(mbc.com_ori), com);
  // mbc.com = sva::PTransformd(com);

  const Eigen::Matrix6d Ic = centroidalInertia(mb, mbc, com);

  const auto Ic_inv = Ic.inverse();
  const Eigen::Matrix6d Ic_dot = centroidalInertiaDot(mb, mbc, com, mbc.comVel.linear());
  const auto h = computeCentroidalMomentum(mb, mbc, com).vector();
  mbc.comVel = sva::MotionVecd(Ic_inv * h);
  // mbc.comVel = sva::MotionVecd(Eigen::Vector3d::Zero(),mbc.comVel.linear());

  mbc.comAcc = sva::MotionVecd(Ic_inv
                               * (computeCentroidalMomentumDot(mb, mbc, com, mbc.comVel.linear()).vector()
                                  - sva::vector6ToCrossDualMatrix<double>(mbc.comVel.vector()) * h));
  mbc.comAcc += sva::MotionVecd(Eigen::Vector3d::Zero(), mbc.comVel.angular().cross(mbc.comVel.linear()));
  // mbc.comAcc = sva::MotionVecd(Eigen::Vector3d::Zero(),mbc.comAcc.linear());

  mbc.comAcc = sva::MotionVecd(
      Ic_inv
      * (computeCentroidalMomentumDot(mb, mbc, com, mbc.comVel.linear()).vector() - Ic_dot * mbc.comVel.vector()));
  mbc.comAcc += sva::MotionVecd(Eigen::Vector3d::Zero(), mbc.comVel.angular().cross(mbc.comVel.linear()));

  rbd::CentroidalMomentumMatrix cmm(mb);
  cmm.computeMatrixAndMatrixDot(mb, mbc, com, mbc.comVel.linear());
  mbc.Jcom = Ic_inv * cmm.matrix();
  mbc.Jcomdot = Ic_inv * (cmm.matrixDot() - sva::vector6ToCrossDualMatrix(mbc.comVel.vector()) * Ic * mbc.Jcom);
}

void updateCoMFrame(const MultiBody & mb, MultiBodyConfig & mbc, const double step)
{
  const Eigen::Vector3d com = computeCoM(mb, mbc);

  auto qd = rbd::dofToVector(mb, mbc.alpha);
  auto qdd = rbd::dofToVector(mb, mbc.alphaD);
  mbc.comAcc = sva::MotionVecd(mbc.Jcom * qdd + mbc.Jcomdot * qd);

  mbc.comAcc = mbc.comAcc + sva::MotionVecd(Eigen::Vector3d::Zero(), mbc.comVel.angular().cross(mbc.comVel.linear()));

  const Eigen::Quaterniond com_ori(mbc.com_ori[0], mbc.com_ori[1], mbc.com_ori[2], mbc.com_ori[3]);

  const auto new_ori = rbd::SO3Integration(com_ori, mbc.comVel.angular(), mbc.comAcc.angular(), step).first;
  const double nq = new_ori.norm();
  mbc.com_ori = {new_ori.w() / nq, new_ori.x() / nq, new_ori.y() / nq, new_ori.z() / nq};
  mbc.com = sva::PTransformd(QuatToE(mbc.com_ori), com);
  // mbc.com = sva::PTransformd(com);

  const Eigen::Matrix6d Ic = centroidalInertia(mb, mbc, com);
  const Eigen::Matrix6d Ic_dot = centroidalInertiaDot(mb, mbc, com, mbc.comVel.linear());

  const auto Ic_inv = Ic.inverse();
  mbc.comVel = sva::MotionVecd(Ic_inv * computeCentroidalMomentum(mb, mbc, com).vector());

  // mbc.comVel = sva::MotionVecd(Eigen::Vector3d::Zero(), mbc.comVel.linear());

  mbc.comAcc = sva::MotionVecd(
      Ic_inv
      * (computeCentroidalMomentumDot(mb, mbc, com, mbc.comVel.linear()).vector() - Ic_dot * mbc.comVel.vector()));
  mbc.comAcc += sva::MotionVecd(Eigen::Vector3d::Zero(), mbc.comVel.angular().cross(mbc.comVel.linear()));

  // std::cout << "//" << std::endl;
  // std::cout << Ic_dot << std::endl;
  // std::cout << Ic_dot * mbc.comVel.vector() << std::endl;

  rbd::CentroidalMomentumMatrix cmm(mb);
  cmm.computeMatrixAndMatrixDot(mb, mbc, com, mbc.comVel.linear());
  
  mbc.Jcom = Ic_inv * cmm.matrix();
  mbc.Jcomdot = Ic_inv * (cmm.matrixDot() - Ic_dot * mbc.Jcom);

  const auto IcJac = rbd::centroidalInertiaJacobian(mb,mbc);
  const auto IcJacDot = rbd::centroidalInertiaJacobianDot(mb,mbc);
  mbc.JIdv = Ic_dot * mbc.Jcom;
  mbc.JIdvDot = Ic_dot * mbc.Jcomdot;
  size_t i = 0;
  for(auto & Jic : IcJac)
  {
    mbc.JIdv.block(0,i,6,1) += Jic * mbc.comVel.vector();
    mbc.JIdvDot.block(0,i,6,1) += IcJacDot[i] * mbc.comVel.vector();
    i+=1;
  }
}

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

Eigen::Vector3d computeCoMVelocity(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();

  Vector3d comV = Vector3d::Zero();
  double totalMass = 0.;

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    double mass = bodies[i].inertia().mass();
    totalMass += mass;

    // Velocity at CoM : com_T_b·V_b
    // Velocity at CoM world frame : 0_R_b·com_T_b·V_b
    sva::PTransformd X_0_i(mbc.bodyPosW[i].rotation().transpose(), bodies[i].inertia().momentum());
    sva::MotionVecd scaledBodyVelB(mbc.bodyVelB[i].angular(), mass * mbc.bodyVelB[i].linear());
    comV += (X_0_i * scaledBodyVelB).linear();
  }

  assert(totalMass > 0 && "Invalid multibody. Totalmass must be strictly positive");
  return comV / totalMass;
}

Eigen::Vector3d computeCoMAcceleration(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();

  Vector3d comA = Vector3d::Zero();
  double totalMass = 0.;

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    double mass = bodies[i].inertia().mass();

    totalMass += mass;

    // Acceleration at CoM : com_T_b·A_b
    // Acceleration at CoM world frame :
    //    0_R_b·com_T_b·A_b + 0_R_b_d·com_T_b·V_b
    // O_R_b_d : (Angvel_W)_b x 0_R_b
    sva::PTransformd X_0_iscaled(mbc.bodyPosW[i].rotation().transpose(), bodies[i].inertia().momentum());
    sva::MotionVecd angvel_W(mbc.bodyVelW[i].angular(), Eigen::Vector3d::Zero());

    sva::MotionVecd scaledBodyAccB(mbc.bodyAccB[i].angular(), mass * mbc.bodyAccB[i].linear());
    sva::MotionVecd scaledBodyVelB(mbc.bodyVelB[i].angular(), mass * mbc.bodyVelB[i].linear());

    comA += (X_0_iscaled * scaledBodyAccB).linear();
    comA += (angvel_W.cross(X_0_iscaled * scaledBodyVelB)).linear();
  }

  assert(totalMass > 0 && "Invalid multibody. Totalmass must be strictly positive");
  return comA / totalMass;
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
: jac_(3, mb.nrDof()), jacDot_(3, mb.nrDof()), jacFull_(3, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
  totalMass_(0.), bodiesWeight_(static_cast<size_t>(mb.nrBodies()), 1.)
{
  init(mb);
}

CoMJacobianDummy::CoMJacobianDummy(const MultiBody & mb, std::vector<double> weight)
: jac_(3, mb.nrDof()), jacDot_(3, mb.nrDof()), jacFull_(3, mb.nrDof()), jacVec_(static_cast<size_t>(mb.nrBodies())),
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

const Eigen::MatrixXd & CoMJacobianDummy::jacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();

  jac_.setZero();

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const MatrixXd & jac = jacVec_[i].jacobian(mb, mbc);
    jacVec_[i].fullJacobian(mb, jac.block(3, 0, 3, jac.cols()), jacFull_);
    jac_.noalias() += jacFull_ * (bodies[i].inertia().mass() * bodiesWeight_[i]);
  }

  assert(totalMass_ > 0 && "Invalid multibody. Totalmass must be strictly positive");
  jac_ /= totalMass_;

  return jac_;
}

const Eigen::MatrixXd & CoMJacobianDummy::jacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  using namespace Eigen;

  const std::vector<Body> & bodies = mb.bodies();

  jacDot_.setZero();

  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const MatrixXd & jac = jacVec_[i].jacobianDot(mb, mbc);
    jacVec_[i].fullJacobian(mb, jac.block(3, 0, 3, jac.cols()), jacFull_);
    jacDot_.noalias() += jacFull_ * (bodies[i].inertia().mass() * bodiesWeight_[i]);
  }

  assert(totalMass_ > 0 && "Invalid multibody. Totalmass must be strictly positive");
  jacDot_ /= totalMass_;

  return jacDot_;
}

const Eigen::MatrixXd & CoMJacobianDummy::sJacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobian(mb, mbc);
}

const Eigen::MatrixXd & CoMJacobianDummy::sJacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
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
: jac_(3, mb.nrDof()), jacDot_(3, mb.nrDof()), bodiesCoeff_(static_cast<size_t>(mb.nrBodies())),
  bodiesCoM_(static_cast<size_t>(mb.nrBodies())), jointsSubBodies_(static_cast<size_t>(mb.nrJoints())),
  bodiesCoMWorld_(static_cast<size_t>(mb.nrBodies())), bodiesCoMVelB_(static_cast<size_t>(mb.nrBodies())),
  normalAcc_(static_cast<size_t>(mb.nrJoints())), weight_(static_cast<size_t>(mb.nrBodies()), 1.)
{
  init(mb);
}

CoMJacobian::CoMJacobian(const MultiBody & mb, std::vector<double> weight)
: jac_(3, mb.nrDof()), jacDot_(3, mb.nrDof()), bodiesCoeff_(static_cast<size_t>(mb.nrBodies())),
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

const Eigen::MatrixXd & CoMJacobian::jacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  const std::vector<Joint> & joints = mb.joints();

  jac_.setZero();

  // we pre compute the CoM position of each bodie in world frame
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    // the transformation must be read {}^0E_p {}^pT_N {}^NX_0
    sva::PTransformd X_0_com_w = bodiesCoM_[i] * mbc.bodyPosW[i];
    bodiesCoMWorld_[i] = sva::PTransformd(X_0_com_w.translation());
  }

  int curJ = 0;
  for(size_t i = 0; i < static_cast<size_t>(mb.nrJoints()); ++i)
  {
    std::vector<int> & subBodies = jointsSubBodies_[i];
    sva::PTransformd X_i_0 = mbc.bodyPosW[i].inv();
    for(int b : subBodies)
    {
      const auto body_index = static_cast<size_t>(b);
      sva::PTransformd X_i_com = bodiesCoMWorld_[body_index] * X_i_0;
      for(int dof = 0; dof < joints[i].dof(); ++dof)
      {
        jac_.col(curJ + dof).noalias() +=
            (X_i_com.linearMul(sva::MotionVecd(mbc.motionSubspace[i].col(dof)))) * bodiesCoeff_[body_index];
      }
    }
    curJ += joints[i].dof();
  }

  return jac_;
}

const Eigen::MatrixXd & CoMJacobian::jacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  const std::vector<Joint> & joints = mb.joints();

  jacDot_.setZero();

  // we pre compute the CoM position/velocity of each bodie
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    bodiesCoMWorld_[i] = bodiesCoM_[i] * mbc.bodyPosW[i];
    bodiesCoMVelB_[i] = bodiesCoM_[i] * mbc.bodyVelB[i];
  }

  int curJ = 0;
  for(size_t i = 0; i < static_cast<size_t>(mb.nrJoints()); ++i)
  {
    std::vector<int> & subBodies = jointsSubBodies_[i];
    sva::PTransformd X_i_0 = mbc.bodyPosW[i].inv();

    for(int b : subBodies)
    {
      const auto body_index = static_cast<size_t>(b);
      sva::PTransformd X_i_com = bodiesCoMWorld_[body_index] * X_i_0;
      sva::PTransformd E_b_0(Eigen::Matrix3d(mbc.bodyPosW[body_index].rotation().transpose()));

      // angular velocity of rotation N to O
      sva::MotionVecd E_Vb(mbc.bodyVelW[body_index].angular(), Eigen::Vector3d::Zero());
      sva::MotionVecd X_Vcom_i_com = X_i_com * mbc.bodyVelB[i] - bodiesCoMVelB_[body_index];

      for(int dof = 0; dof < joints[i].dof(); ++dof)
      {
        sva::MotionVecd S_ij(mbc.motionSubspace[i].col(dof));

        // JD_i = (E_com_0_d*X_i_com*S_i + E_com_0*X_i_com_d*S_i)*(mass/totalMass)
        // E_com_0_d = (ANG_Vcom)_0 x E_com_0
        // X_i_com_d = (Vi - Vcom)_com x X_i_com
        jacDot_.col(curJ + dof).noalias() +=
            ((E_Vb.cross(E_b_0 * X_i_com * S_ij)).linear() + (E_b_0 * X_Vcom_i_com.cross(X_i_com * S_ij)).linear())
            * bodiesCoeff_[body_index];
      }
    }
    curJ += joints[i].dof();
  }

  return jacDot_;
}

Eigen::Vector3d CoMJacobian::velocity(const MultiBody & mb, const MultiBodyConfig & mbc) const
{
  Eigen::Vector3d comV = Eigen::Vector3d::Zero();
  for(size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    const Eigen::Vector3d & comT = bodiesCoM_[i].translation();

    // Velocity at CoM : com_T_b·V_b
    // Velocity at CoM world frame : 0_R_b·com_T_b·V_b
    sva::PTransformd X_0_i(mbc.bodyPosW[i].rotation().transpose(), comT);
    comV += (X_0_i * mbc.bodyVelB[i]).linear() * bodiesCoeff_[i];
  }

  return comV;
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

const Eigen::MatrixXd & CoMJacobian::sJacobian(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchBodyPos(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);

  return jacobian(mb, mbc);
}

const Eigen::MatrixXd & CoMJacobian::sJacobianDot(const MultiBody & mb, const MultiBodyConfig & mbc)
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
