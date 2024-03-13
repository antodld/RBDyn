/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_MODULE MomentumTest
#include <boost/test/unit_test.hpp>

// RBDyn
#include "RBDyn/CoM.h"
#include "RBDyn/FA.h"
#include "RBDyn/FK.h"
#include "RBDyn/FV.h"
#include "RBDyn/Momentum.h"
#include "RBDyn/MultiBody.h"
#include "RBDyn/MultiBodyConfig.h"
#include "RBDyn/NumericalIntegration.h"

// arm
#include "XYZSarm.h"

// tree
#include "Tree30Dof.h"

const double TOL = 1e-6;

BOOST_AUTO_TEST_CASE(centroidalMomentum)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;
  std::tie(mb, mbc, mbg) = makeTree30Dof(false);

  VectorXd q(mb.nrParams());
  VectorXd alpha(mb.nrDof());
  CentroidalMomentumMatrix cmm(mb);

  {
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::paramToVector(mbc.alpha, alpha);
    rbd::resetCoMFrame(mb, mbc);

    Vector3d com = rbd::computeCoM(mb, mbc);
    ForceVecd momentum = rbd::computeCentroidalMomentum(mb, mbc, com);
    cmm.computeMatrix(mb, mbc, com);

    const ForceVecd momentumM(cmm.matrix() * alpha);

    BOOST_CHECK_EQUAL(momentum.vector().norm(), 0.);
    BOOST_CHECK_EQUAL(momentumM.vector().norm(), 0.);
  }

  // test J·q against computeCentroidalMomentum
  for(int i = 0; i < 100; ++i)
  {
    q.setRandom();
    q.segment(0, 4).normalize();
    alpha.setRandom();
    rbd::vectorToParam(q, mbc.q);
    rbd::vectorToParam(alpha, mbc.alpha);

    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::resetCoMFrame(mb, mbc);

    const auto E_com_0 = sva::PTransformd(mbc.com.rotation()).inv();
    Vector3d com = rbd::computeCoM(mb, mbc);
    Vector3d comVel = rbd::computeCoMVelocity(mb, mbc);
    ForceVecd momentum = rbd::computeCentroidalMomentum(mb, mbc, com);
    cmm.computeMatrix(mb, mbc, com);

    ForceVecd momentumM(cmm.matrix() * alpha);

    BOOST_CHECK_SMALL((momentum - momentumM).vector().norm(), TOL);
    BOOST_CHECK_SMALL((comVel - (E_com_0 * mbc.comVel).linear()).norm(), TOL);
  }

  // test J·q against CentroidalMomentumMatrix::momentum
  for(int i = 0; i < 50; ++i)
  {
    std::vector<double> weight(static_cast<size_t>(mb.nrBodies()));
    for(std::size_t i = 0; i < weight.size(); ++i)
    {
      weight[i] = Eigen::Matrix<double, 1, 1>::Random()(0);
    }

    CentroidalMomentumMatrix cmmW(mb, weight);

    q.setRandom();
    q.segment(0, 4).normalize();
    alpha.setRandom();

    rbd::vectorToParam(q, mbc.q);
    rbd::vectorToParam(alpha, mbc.alpha);

    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::resetCoMFrame(mb, mbc);

    Vector3d com = rbd::computeCoM(mb, mbc);
    ForceVecd momentum = cmmW.momentum(mb, mbc, com);
    cmmW.computeMatrix(mb, mbc, com);

    ForceVecd momentumM(cmmW.matrix() * alpha);

    BOOST_CHECK_SMALL((momentum - momentumM).vector().norm(), TOL);
  }
}

BOOST_AUTO_TEST_CASE(centroidalMomentumDot)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;
  std::tie(mb, mbc, mbg) = makeTree30Dof(false);

  CentroidalMomentumMatrix cmm(mb);

  VectorXd q(mb.nrParams());
  VectorXd alpha(mb.nrDof());
  VectorXd alphaD(mb.nrDof());

  const double step = 1e-8;

  for(int i = 0; i < 10; ++i)
  {
    q.setRandom();
    q.segment(0, 4).normalize();
    alpha.setRandom();
    alphaD.setRandom();
    rbd::vectorToParam(q, mbc.q);
    rbd::vectorToParam(alpha, mbc.alpha);
    rbd::vectorToParam(alphaD, mbc.alphaD);

    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    rbd::resetCoMFrame(mb, mbc);
    for(int j = 0; j < 10; ++j)
    {
      Eigen::Vector3d oldCom = rbd::computeCoM(mb, mbc);
      ForceVecd oldMomentum = rbd::computeCentroidalMomentum(mb, mbc, oldCom);
      Eigen::Vector3d oldComVel = rbd::computeCoMVelocity(mb, mbc);
      Eigen::Vector3d oldComAcc = rbd::computeCoMAcceleration(mb, mbc);

      sva::MotionVecd oldComVel6D = mbc.comVel;
      sva::MotionVecd oldComAcc6D = mbc.comAcc;

      // Eigen::Matrix6d oldIc = rbd::centroidalInertia(mb, mbc, oldCom);
      // Eigen::Matrix6d IcDot = rbd::centroidalInertiaDot(mb, mbc, oldCom, oldComVel);

      sva::PTransformd E_com_0 = sva::PTransformd(mbc.com.rotation()).inv();

      ForceVecd momentumDot = rbd::computeCentroidalMomentumDot(mb, mbc, oldCom, oldComVel);
      cmm.computeMatrix(mb, mbc, oldCom);
      cmm.computeMatrixDot(mb, mbc, oldCom, oldComVel);
      Eigen::MatrixXd cmmMatrix = cmm.matrix();
      Eigen::MatrixXd cmmMatrixDot = cmm.matrixDot();
      Vector6d momentumDotCMM = cmmMatrix * alphaD + cmmMatrixDot * alpha;

      // check that the momentum are the same
      BOOST_CHECK_SMALL((momentumDot.vector() - momentumDotCMM).norm(), TOL);

      // check that each compute compute the same thing
      cmm.computeMatrixAndMatrixDot(mb, mbc, oldCom, oldComVel);

      BOOST_CHECK_SMALL((cmmMatrix - cmm.matrix()).norm(), TOL);
      BOOST_CHECK_SMALL((cmmMatrixDot - cmm.matrixDot()).norm(), TOL);
      BOOST_CHECK_SMALL(((E_com_0 * oldComVel6D).linear() - oldComVel).norm(), TOL);
      BOOST_CHECK_SMALL(((E_com_0 * oldComAcc6D).linear() - oldComAcc).norm(), TOL);

      oldComAcc6D = E_com_0 * oldComAcc6D;
      oldComVel6D = E_com_0 * oldComVel6D;
      oldMomentum = mbc.com.transMul(oldMomentum);
      momentumDot = mbc.com.transMul(momentumDot);

      rbd::integration(mb, mbc, step);

      rbd::forwardKinematics(mb, mbc);
      rbd::forwardVelocity(mb, mbc);
      rbd::forwardAcceleration(mb, mbc);
      rbd::updateCoMFrame(mb, mbc, step);

      rbd::paramToVector(mbc.alpha, alpha);
      rbd::paramToVector(mbc.alphaD, alphaD);

      E_com_0 = sva::PTransformd(mbc.com.rotation()).inv();

      Vector3d newCom = rbd::computeCoM(mb, mbc);
      // MotionVecd comAcc6DDiff = (E_com_0 * mbc.comVel - oldComVel6D) * (1 / step);
      auto comAccDiff = (rbd::computeCoMVelocity(mb, mbc) - oldComVel) / step;
      // Matrix6d newIc = rbd::centroidalInertia(mb, mbc, newCom);
      // Matrix6d newIcDiff = oldIc + IcDot * step;

      ForceVecd newMomentum = rbd::computeCentroidalMomentum(mb, mbc, newCom);
      ForceVecd momentumDotDiff = (mbc.com.transMul(newMomentum) - oldMomentum) * (1. / step);

      BOOST_CHECK_SMALL((momentumDot - momentumDotDiff).vector().norm(), 10 * TOL);
      // BOOST_CHECK_SMALL((newIcDiff - newIc).norm(), TOL);
      // BOOST_CHECK_SMALL((comAcc6DDiff - oldComAcc6D).angular().norm(), TOL);
      BOOST_CHECK_SMALL((comAccDiff - oldComAcc).norm(), TOL);
    }
  }

  // test JDot·q against CentroidalMomentumMatrix::normalMomentumDot
  for(int i = 0; i < 50; ++i)
  {
    std::vector<double> weight(static_cast<size_t>(mb.nrBodies()));
    for(std::size_t i = 0; i < weight.size(); ++i)
    {
      weight[i] = Eigen::Matrix<double, 1, 1>::Random()(0);
    }

    CentroidalMomentumMatrix cmmW(mb, weight);

    q.setRandom();
    q.segment(0, 4).normalize();
    // q.segment<4>(mb.jointPosInParam(mb.jointIndexByName("j3"))).normalize();
    alpha.setRandom();
    alphaD.setZero();
    rbd::vectorToParam(q, mbc.q);
    rbd::vectorToParam(alpha, mbc.alpha);
    // calcul the normal acceleration since alphaD is zero
    rbd::vectorToParam(alphaD, mbc.alphaD);

    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    rbd::resetCoMFrame(mb, mbc);

    Vector3d com = rbd::computeCoM(mb, mbc);
    Vector3d comDot = rbd::computeCoMVelocity(mb, mbc);
    ForceVecd normalMomentumDot1 = cmmW.normalMomentumDot(mb, mbc, com, comDot);
    ForceVecd normalMomentumDot2 = cmmW.normalMomentumDot(mb, mbc, com, comDot, mbc.bodyAccB);
    cmmW.computeMatrixDot(mb, mbc, com, comDot);

    ForceVecd normalMomentumDotM(cmmW.matrixDot() * alpha);

    BOOST_CHECK_SMALL((normalMomentumDot1 - normalMomentumDotM).vector().norm(), TOL);
    BOOST_CHECK_SMALL((normalMomentumDot2 - normalMomentumDotM).vector().norm(), TOL);
  }
}
