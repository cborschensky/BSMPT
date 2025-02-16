// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include <algorithm>              // for max, copy
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <cstddef>               // for std::size_t

#include <BSMPT/models/ClassPotentialRxSM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT
{
namespace Models
{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of
 * Lagrangian parameters AFTER using the tadpole conditions), nParCT (number of
 * counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_RxSM::Class_RxSM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model = ModelID::ModelIDs::RXSMNOZ2; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 3; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 6; // number of independent input parameters (in the tree-Level Lagrangian)
  nParCT = 12; // number of parameters in the counterterm potential

  nVEV = 2; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 3;
  VevOrder[1] = 4;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_RxSM::~Class_RxSM()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_RxSM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmusq");
  labels.push_back("dlam");
  labels.push_back("da1");
  labels.push_back("da2");
  labels.push_back("db2");
  labels.push_back("db3");
  labels.push_back("db4");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_RxSM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("omega_c"); // Label for the critical vev
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_H(T_c)");
  labels.push_back("omega_S(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_RxSM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "h_SM";
  particles[4] = "h_H";

  std::string out = "Tree_";
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_RxSM::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega_H");
  labels.push_back("omega_S");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_RxSM::ReadAndSet(const std::string &linestr,
                                std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lamIn{0}, a1In{0}, a2In{0}, b3In{0}, b4In{0}, vSIn{0};

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 6; ++k)
  {
    ss >> tmp;
    if (k == 1)
      lamIn = tmp;
    else if (k == 2)
      a1In = tmp;
    else if (k == 3)
      a2In = tmp;
    else if (k == 4)
      b3In = tmp;
    else if (k == 5)
      b4In = tmp;
    else if (k == 6)
      vSIn = tmp;
  }
  par[0] = lamIn;
  par[1] = a1In;
  par[2] = a2In;
  par[3] = b3In;
  par[4] = b4In;
  par[5] = vSIn;

  set_gen(par); // This you have to call so that everything will be set
  return;

}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_RxSM::set_gen(const std::vector<double> &par)
{
  lam = par[0];
  a1  = par[1];
  a2  = par[2];
  b3  = par[3];
  b4  = par[4];
  vS  = par[5];

  const double ZeroThreshold = 1e-5;

  UnbrokenSingletPhase = false;
  if (std::abs(vS) < ZeroThreshold)
  {
    std::cout << "Currently not implemented." << std::endl;
  }

  vH = SMConstants.C_vev0;

  // Tree-level tadpole equations
  musq = (1.0/2.0)*a1*vS + (1.0/2.0)*a2*std::pow(vS, 2) + lam*std::pow(vH, 2);
  b2 = -1.0/4.0*a1*std::pow(vH, 2)/vS - 1.0/2.0*a2*std::pow(vH, 2) - b3*vS - b4*std::pow(vS, 2);

  scale = vH; // Renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  vevTreeMin[0] = vH;
  vevTreeMin[1] = vS;
  vevTree = MinimizeOrderVEV(vevTreeMin);

  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_RxSM::set_CT_Pot_Par(const std::vector<double> &par)
{

  dmusq = par.at(0);
  dlam = par.at(1);
  da1 = par.at(2);
  da2 = par.at(3);
  db2 = par.at(4);
  db3 = par.at(5);
  db4 = par.at(6);
  dT1 = par.at(7);
  dT2 = par.at(8);
  dT3 = par.at(9);
  dT4 = par.at(10);
  dT5 = par.at(11);

  Curvature_Higgs_CT_L1.at(0) = dT1;
  Curvature_Higgs_CT_L1.at(1) = dT2;
  Curvature_Higgs_CT_L1.at(2) = dT3;
  Curvature_Higgs_CT_L1.at(3) = dT4;
  Curvature_Higgs_CT_L1.at(4) = dT5;

  Curvature_Higgs_CT_L2.at(0).at(0) = -dmusq;
  Curvature_Higgs_CT_L2.at(1).at(1) = -dmusq;
  Curvature_Higgs_CT_L2.at(2).at(2) = -dmusq;
  Curvature_Higgs_CT_L2.at(3).at(3) = -dmusq;
  Curvature_Higgs_CT_L2.at(4).at(4) = db2;

  Curvature_Higgs_CT_L3.at(0).at(0).at(4) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(0).at(4).at(0) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(1).at(1).at(4) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(1).at(4).at(1) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(2).at(2).at(4) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(2).at(4).at(2) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(3).at(3).at(4) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(3).at(4).at(3) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(4).at(0).at(0) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(4).at(1).at(1) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(4).at(2).at(2) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(4).at(3).at(3) = (1.0/2.0)*da1;
  Curvature_Higgs_CT_L3.at(4).at(4).at(4) = 2*db3;

  Curvature_Higgs_CT_L4.at(0).at(0).at(0).at(0) = 6*dlam;
  Curvature_Higgs_CT_L4.at(0).at(0).at(1).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(0).at(2).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(0).at(3).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(0).at(4).at(4) = da2;
  Curvature_Higgs_CT_L4.at(0).at(1).at(0).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(1).at(1).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(2).at(0).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(2).at(2).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(3).at(0).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(3).at(3).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(0).at(4).at(0).at(4) = da2;
  Curvature_Higgs_CT_L4.at(0).at(4).at(4).at(0) = da2;
  Curvature_Higgs_CT_L4.at(1).at(0).at(0).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(0).at(1).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(1).at(0).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(1).at(1).at(1) = 6*dlam;
  Curvature_Higgs_CT_L4.at(1).at(1).at(2).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(1).at(3).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(1).at(4).at(4) = da2;
  Curvature_Higgs_CT_L4.at(1).at(2).at(1).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(2).at(2).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(3).at(1).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(3).at(3).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(1).at(4).at(1).at(4) = da2;
  Curvature_Higgs_CT_L4.at(1).at(4).at(4).at(1) = da2;
  Curvature_Higgs_CT_L4.at(2).at(0).at(0).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(0).at(2).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(1).at(1).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(1).at(2).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(2).at(0).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(2).at(1).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(2).at(2).at(2) = 6*dlam;
  Curvature_Higgs_CT_L4.at(2).at(2).at(3).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(2).at(4).at(4) = da2;
  Curvature_Higgs_CT_L4.at(2).at(3).at(2).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(3).at(3).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(2).at(4).at(2).at(4) = da2;
  Curvature_Higgs_CT_L4.at(2).at(4).at(4).at(2) = da2;
  Curvature_Higgs_CT_L4.at(3).at(0).at(0).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(0).at(3).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(1).at(1).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(1).at(3).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(2).at(2).at(3) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(2).at(3).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(3).at(0).at(0) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(3).at(1).at(1) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(3).at(2).at(2) = 2*dlam;
  Curvature_Higgs_CT_L4.at(3).at(3).at(3).at(3) = 6*dlam;
  Curvature_Higgs_CT_L4.at(3).at(3).at(4).at(4) = da2;
  Curvature_Higgs_CT_L4.at(3).at(4).at(3).at(4) = da2;
  Curvature_Higgs_CT_L4.at(3).at(4).at(4).at(3) = da2;
  Curvature_Higgs_CT_L4.at(4).at(0).at(0).at(4) = da2;
  Curvature_Higgs_CT_L4.at(4).at(0).at(4).at(0) = da2;
  Curvature_Higgs_CT_L4.at(4).at(1).at(1).at(4) = da2;
  Curvature_Higgs_CT_L4.at(4).at(1).at(4).at(1) = da2;
  Curvature_Higgs_CT_L4.at(4).at(2).at(2).at(4) = da2;
  Curvature_Higgs_CT_L4.at(4).at(2).at(4).at(2) = da2;
  Curvature_Higgs_CT_L4.at(4).at(3).at(3).at(4) = da2;
  Curvature_Higgs_CT_L4.at(4).at(3).at(4).at(3) = da2;
  Curvature_Higgs_CT_L4.at(4).at(4).at(0).at(0) = da2;
  Curvature_Higgs_CT_L4.at(4).at(4).at(1).at(1) = da2;
  Curvature_Higgs_CT_L4.at(4).at(4).at(2).at(2) = da2;
  Curvature_Higgs_CT_L4.at(4).at(4).at(3).at(3) = da2;
  Curvature_Higgs_CT_L4.at(4).at(4).at(4).at(4) = 6*db4;

}

/**
 * console output of all Parameters
 */
void Class_RxSM::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "lam  = " << lam << std::endl
     << "a_1  = " << a1 << std::endl
     << "a_2  = " << a2 << std::endl
     << "b_3  = " << b3 << std::endl
     << "b_4  = " << b4 << std::endl
     << "v_S  = " << vS << std::endl
     << "v_H  = " << vH << " (fixed to SM value)" << std::endl
     << "mu^2 = " << musq << " (via tadpole eqs.)" << std::endl
     << "b_2  = " << b2 << " (via tadpole eqs.)" << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
  ss << "dmu^2 = " << dmusq << std::endl
     << "dlam  = " << dlam << std::endl
     << "da1   = " << da1 << std::endl
     << "da2   = " << da2 << std::endl
     << "db2   = " << db2 << std::endl
     << "db3   = " << db3 << std::endl
     << "db4   = " << db4 << std::endl
     << "dT1   = " << dT1 << std::endl
     << "dT2   = " << dT2 << std::endl
     << "dT3   = " << dT3 << std::endl
     << "dT4   = " << dT4 << std::endl
     << "dT5   = " << dT5 << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrixEnsuredConvention[i][j];
    }
  }

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n"
     << "m_{G^+}^2 = " << HiggsMasses[pos_G1] << " GeV \n"
     << "m_{G^-}^2 = " << HiggsMasses[pos_G2] << " GeV \n"
     << "m_{G^0}^2 = " << HiggsMasses[pos_G0] << " GeV \n"
     << "m_h_SM    = " << std::sqrt(HiggsMasses[pos_h_SM]) << " GeV \n"
     << "m_h_H     = " << std::sqrt(HiggsMasses[pos_h_H]) << " GeV \n";

  ss << "The neutral mixing Matrix is given by :\n";
  ss << "h_SM = " << HiggsRot(pos_h_SM, 3) << " zeta_1 ";
  bool IsNegative = HiggsRot(pos_h_SM, 4) < 0;
  if (IsNegative)
    ss << "- ";
  else
    ss << "+ ";
  ss << std::abs(HiggsRot(pos_h_SM, 4)) << " zeta_S\n"
     << "h_H  = " << HiggsRot(pos_h_H, 3) << " zeta_1 ";
  IsNegative = HiggsRot(pos_h_H, 4) < 0;
  if (IsNegative)
    ss << "- ";
  else
    ss << "+ ";
  ss << std::abs(HiggsRot(pos_h_H, 4)) << " zeta_S" << std::endl;
  ss << "The mixing angle is: alpha = " << alpha << std::endl;

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_RxSM::calc_CT() const
{

  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  if (UnbrokenSingletPhase)
  {
    std::cout << "Currently not implemented." << std::endl;
  }
  else
  {
    // Free parameter t chosen such that all tadpole CTs vanish

    // // ren-1: db3 = db4 = 0
    // parCT.push_back((3.0/2.0)*HesseWeinberg(2,2) - 1.0/2.0*HesseWeinberg(3,3) - 1.0/2.0*HesseWeinberg(3,4)*vS/vH + HesseWeinberg(4,4)*std::pow(vS, 2)/std::pow(vH, 2) - NablaWeinberg(4)*vS/std::pow(vH, 2)); //dmusq
    // parCT.push_back((1.0/2.0)*HesseWeinberg(2,2)/std::pow(vH, 2) - 1.0/2.0*HesseWeinberg(3,3)/std::pow(vH, 2)); //dlam
    // parCT.push_back(4*HesseWeinberg(4,4)*vS/std::pow(vH, 2) - 4*NablaWeinberg(4)/std::pow(vH, 2)); //da1
    // parCT.push_back(-HesseWeinberg(3,4)/(vH*vS) - 2*HesseWeinberg(4,4)/std::pow(vH, 2) + 2*NablaWeinberg(4)/(std::pow(vH, 2)*vS)); //da2
    // parCT.push_back((1.0/2.0)*HesseWeinberg(3,4)*vH/vS - NablaWeinberg(4)/vS); //db2
    // parCT.push_back(0); //db3
    // parCT.push_back(0); //db4
    // parCT.push_back(-NablaWeinberg(0)); //dT1
    // parCT.push_back(-NablaWeinberg(1)); //dT2
    // parCT.push_back(-NablaWeinberg(2)); //dT3
    // parCT.push_back(HesseWeinberg(2,2)*vH - NablaWeinberg(3)); //dT4
    // parCT.push_back(0); //dT5

    // // ren-2: dmusq = db2 = 0
    // parCT.push_back(0); //dmusq
    // parCT.push_back((1.0/2.0)*HesseWeinberg(2,2)/std::pow(vH, 2) - 1.0/2.0*HesseWeinberg(3,3)/std::pow(vH, 2)); //dlam
    // parCT.push_back(-6*HesseWeinberg(2,2)/vS + 2*HesseWeinberg(3,3)/vS + 2*HesseWeinberg(3,4)/vH); //da1
    // parCT.push_back(3*HesseWeinberg(2,2)/std::pow(vS, 2) - HesseWeinberg(3,3)/std::pow(vS, 2) - 2*HesseWeinberg(3,4)/(vH*vS)); //da2
    // parCT.push_back(0); //db2
    // parCT.push_back((3.0/2.0)*HesseWeinberg(2,2)*std::pow(vH, 2)/std::pow(vS, 3) - 1.0/2.0*HesseWeinberg(3,3)*std::pow(vH, 2)/std::pow(vS, 3) + (1.0/2.0)*HesseWeinberg(3,4)*vH/std::pow(vS, 2) + HesseWeinberg(4,4)/vS - 3*NablaWeinberg(4)/std::pow(vS, 2)); //db3
    // parCT.push_back(-3.0/2.0*HesseWeinberg(2,2)*std::pow(vH, 2)/std::pow(vS, 4) + (1.0/2.0)*HesseWeinberg(3,3)*std::pow(vH, 2)/std::pow(vS, 4) - HesseWeinberg(4,4)/std::pow(vS, 2) + 2*NablaWeinberg(4)/std::pow(vS, 3)); //db4
    // parCT.push_back(-NablaWeinberg(0)); //dT1
    // parCT.push_back(-NablaWeinberg(1)); //dT2
    // parCT.push_back(-NablaWeinberg(2)); //dT3
    // parCT.push_back(HesseWeinberg(2,2)*vH - NablaWeinberg(3)); //dT4
    // parCT.push_back(0); //dT5

    // // ren-3: da2 = db3 = 0
    // parCT.push_back((3.0/2.0)*HesseWeinberg(2,2) - 1.0/2.0*HesseWeinberg(3,3) - HesseWeinberg(3,4)*vS/vH); //dmusq
    // parCT.push_back((1.0/2.0)*HesseWeinberg(2,2)/std::pow(vH, 2) - 1.0/2.0*HesseWeinberg(3,3)/std::pow(vH, 2)); //dlam
    // parCT.push_back(-2*HesseWeinberg(3,4)/vH); //da1
    // parCT.push_back(0); //da2
    // parCT.push_back((3.0/4.0)*HesseWeinberg(3,4)*vH/vS + (1.0/2.0)*HesseWeinberg(4,4) - 3.0/2.0*NablaWeinberg(4)/vS); //db2
    // parCT.push_back(0); //db3
    // parCT.push_back(-1.0/4.0*HesseWeinberg(3,4)*vH/std::pow(vS, 3) - 1.0/2.0*HesseWeinberg(4,4)/std::pow(vS, 2) + (1.0/2.0)*NablaWeinberg(4)/std::pow(vS, 3)); //db4
    // parCT.push_back(-NablaWeinberg(0)); //dT1
    // parCT.push_back(-NablaWeinberg(1)); //dT2
    // parCT.push_back(-NablaWeinberg(2)); //dT3
    // parCT.push_back(HesseWeinberg(2,2)*vH - NablaWeinberg(3)); //dT4
    // parCT.push_back(0); //dT5

    // ren-4: da2 = db4 = 0
    parCT.push_back((3.0/2.0)*HesseWeinberg(2,2) - 1.0/2.0*HesseWeinberg(3,3) - HesseWeinberg(3,4)*vS/vH); //dmusq
    parCT.push_back((1.0/2.0)*HesseWeinberg(2,2)/std::pow(vH, 2) - 1.0/2.0*HesseWeinberg(3,3)/std::pow(vH, 2)); //dlam
    parCT.push_back(-2*HesseWeinberg(3,4)/vH); //da1
    parCT.push_back(0); //da2
    parCT.push_back(HesseWeinberg(3,4)*vH/vS + HesseWeinberg(4,4) - 2*NablaWeinberg(4)/vS); //db2
    parCT.push_back(-1.0/2.0*HesseWeinberg(3,4)*vH/std::pow(vS, 2) - HesseWeinberg(4,4)/vS + NablaWeinberg(4)/std::pow(vS, 2)); //db3
    parCT.push_back(0); //db4
    parCT.push_back(-NablaWeinberg(0)); //dT1
    parCT.push_back(-NablaWeinberg(1)); //dT2
    parCT.push_back(-NablaWeinberg(2)); //dT3
    parCT.push_back(HesseWeinberg(2,2)*vH - NablaWeinberg(3)); //dT4
    parCT.push_back(0); //dT5

  }

  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_RxSM::AdjustRotationMatrix()
{
  const double ZeroThreshold = 1e-5;

  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (!CheckRotationMatrix()) // Check whether generically generated rotation
                              // matrix is proper rotation matrix
  {
    throw std::runtime_error("Error in rotation matrix.");
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  // initialize position indices (new initialization for each point in multiline
  // files)
  pos_G1 = -1, pos_G2 = -1, pos_G0 = -1, pos_h = -1, pos_H = -1;

  // Doublet: Phi1 = 1/Sqrt(2)*{rho1 + I*eta1, zeta1 + vH + I*psi1}
  // Singlet: S = zetaS + vS
  // higgsbasis = {rho1, eta1, psi1, zeta1, zetaS}

  // interaction basis
  // rho1, eta1, psi1, zeta1, zetaS
  int pos_rho1 = 0, pos_eta1 = 1, pos_psi1 = 2,
      pos_zeta1 = 3, pos_zetaS = 4;

  for (std::size_t i = 0; i < NHiggs;
       i++) // mass base index i corresponds to mass vector sorted in ascending
            // mass
  {
    // The Goldstones are all on the diagonal, there is no mixing
    if (std::abs(HiggsRot(i, pos_rho1)) > ZeroThreshold)
    {
      pos_G1 = i;
    }
    else if (std::abs(HiggsRot(i, pos_eta1)) > ZeroThreshold)
    {
      pos_G2 = i;
    }
    else if (std::abs(HiggsRot(i, pos_psi1)) > ZeroThreshold)
    {
      pos_G0 = i;
    }
    else if (std::abs(HiggsRot(i, pos_zeta1)) + std::abs(HiggsRot(i, pos_zetaS)) >
        ZeroThreshold) // use that mh < mH
    {
      if (pos_h == -1)
      {
        pos_h = i;
      }
      else
      {
        pos_H = i;
      }
    }
  }

  // check if all position indices are set
  if (pos_G1 == -1 or pos_G2 == -1 or pos_G0 == -1 or
      pos_h == -1 or pos_H == -1)
  {
    throw std::runtime_error("Error. Not all position indices are set.");
  }

  // check if all other elements of rotation matrix are zero
  bool zero_element = false;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      int ii = int(i);
      int jj = int(j);
      if (not((jj == pos_rho1 and ii == pos_G1) or
              (jj == pos_eta1 and ii == pos_G2) or
              (jj == pos_psi1 and ii == pos_G0) or
              (jj == pos_zeta1 and (ii == pos_h or ii == pos_H)) or
              (jj == pos_zetaS and (ii == pos_h or ii == pos_H))))
      {
        zero_element = true;
      }
      if (zero_element and std::abs(HiggsRot(i, j)) > ZeroThreshold)
      {
        throw std::runtime_error("Error. Invalid rotation matrix detected.");
      }
      zero_element = false;
    }
  }

  // Determine the additional indices for the SM-like
  // and exotic Higgses
  pos_h_SM = -1, pos_h_H = -1;

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  // Due to the masses being ordered, we will always have
  //  HiggsMasses[pos_h] <= HiggsMasses[pos_H]
  double diff1 = std::abs(std::sqrt(HiggsMasses[pos_h])
                          - SMConstants.C_MassSMHiggs);
  double diff2 = std::abs(std::sqrt(HiggsMasses[pos_H])
                          - SMConstants.C_MassSMHiggs);
  if (diff1 < diff2)
  {
    pos_h_SM = pos_h;
    pos_h_H = pos_H;
  }
  else
  {
    pos_h_H = pos_h;
    pos_h_SM = pos_H;
  }

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = HiggsRot.row(i);
  }

  // charged submatrix
  if (HiggsRotFixed(pos_G1, pos_rho1) < 0) // G1 rho1 (= +1)
  {
    HiggsRotFixed.row(pos_G1) *= -1;
  }
  if (HiggsRotFixed(pos_G2, pos_eta1) < 0) // G2 eta1 (= +1)
  {
    HiggsRotFixed.row(pos_G2) *= -1;
  }

  // check neutral, CP-odd submatrix
  if (HiggsRotFixed(pos_G0, pos_psi1) < 0) // G0 psi1 (= +1)
  {
    HiggsRotFixed.row(pos_G0) *= -1;
  }

  // // check neutral, CP-even submatrix
  if (HiggsRotFixed(pos_h, pos_zeta1) < 0) // h zeta1 (+ cos(alpha))
  {
    // if negative, rotate h
    HiggsRotFixed.row(pos_h) *= -1; // h
  }
  if (HiggsRotFixed(pos_H, pos_zetaS) < 0) // H zetaS (+ cos(alpha))
  {
    // if negative, rotate H
    HiggsRotFixed.row(pos_H) *= -1; // H
  }

  // Extract the fixed mixing angle
  alpha = std::asin(HiggsRotFixed(pos_h, pos_zetaS)); // h zetaS (+ sin(alpha))

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }

  return;
}

void Class_RxSM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (CalculatedTripleCopulings) return;
  CalculatedTripleCopulings = true;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrixEnsuredConvention[i][j];
    }
  }

  std::vector<double> HiggsOrder(NHiggs);
  HiggsOrder[0] = pos_G1;
  HiggsOrder[1] = pos_G2;
  HiggsOrder[2] = pos_G0;
  HiggsOrder[3] = pos_h_SM;
  HiggsOrder[4] = pos_h_H;

  MatrixXd HiggsRotSort(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_RxSM::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2.at(0).at(0) = -musq;
  Curvature_Higgs_L2.at(1).at(1) = -musq;
  Curvature_Higgs_L2.at(2).at(2) = -musq;
  Curvature_Higgs_L2.at(3).at(3) = -musq;
  Curvature_Higgs_L2.at(4).at(4) = b2;

  Curvature_Higgs_L3.at(0).at(0).at(4) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(0).at(4).at(0) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(1).at(1).at(4) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(1).at(4).at(1) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(2).at(2).at(4) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(2).at(4).at(2) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(3).at(3).at(4) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(3).at(4).at(3) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(4).at(0).at(0) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(4).at(1).at(1) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(4).at(2).at(2) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(4).at(3).at(3) = (1.0/2.0)*a1;
  Curvature_Higgs_L3.at(4).at(4).at(4) = 2*b3;

  Curvature_Higgs_L4.at(0).at(0).at(0).at(0) = 6*lam;
  Curvature_Higgs_L4.at(0).at(0).at(1).at(1) = 2*lam;
  Curvature_Higgs_L4.at(0).at(0).at(2).at(2) = 2*lam;
  Curvature_Higgs_L4.at(0).at(0).at(3).at(3) = 2*lam;
  Curvature_Higgs_L4.at(0).at(0).at(4).at(4) = a2;
  Curvature_Higgs_L4.at(0).at(1).at(0).at(1) = 2*lam;
  Curvature_Higgs_L4.at(0).at(1).at(1).at(0) = 2*lam;
  Curvature_Higgs_L4.at(0).at(2).at(0).at(2) = 2*lam;
  Curvature_Higgs_L4.at(0).at(2).at(2).at(0) = 2*lam;
  Curvature_Higgs_L4.at(0).at(3).at(0).at(3) = 2*lam;
  Curvature_Higgs_L4.at(0).at(3).at(3).at(0) = 2*lam;
  Curvature_Higgs_L4.at(0).at(4).at(0).at(4) = a2;
  Curvature_Higgs_L4.at(0).at(4).at(4).at(0) = a2;
  Curvature_Higgs_L4.at(1).at(0).at(0).at(1) = 2*lam;
  Curvature_Higgs_L4.at(1).at(0).at(1).at(0) = 2*lam;
  Curvature_Higgs_L4.at(1).at(1).at(0).at(0) = 2*lam;
  Curvature_Higgs_L4.at(1).at(1).at(1).at(1) = 6*lam;
  Curvature_Higgs_L4.at(1).at(1).at(2).at(2) = 2*lam;
  Curvature_Higgs_L4.at(1).at(1).at(3).at(3) = 2*lam;
  Curvature_Higgs_L4.at(1).at(1).at(4).at(4) = a2;
  Curvature_Higgs_L4.at(1).at(2).at(1).at(2) = 2*lam;
  Curvature_Higgs_L4.at(1).at(2).at(2).at(1) = 2*lam;
  Curvature_Higgs_L4.at(1).at(3).at(1).at(3) = 2*lam;
  Curvature_Higgs_L4.at(1).at(3).at(3).at(1) = 2*lam;
  Curvature_Higgs_L4.at(1).at(4).at(1).at(4) = a2;
  Curvature_Higgs_L4.at(1).at(4).at(4).at(1) = a2;
  Curvature_Higgs_L4.at(2).at(0).at(0).at(2) = 2*lam;
  Curvature_Higgs_L4.at(2).at(0).at(2).at(0) = 2*lam;
  Curvature_Higgs_L4.at(2).at(1).at(1).at(2) = 2*lam;
  Curvature_Higgs_L4.at(2).at(1).at(2).at(1) = 2*lam;
  Curvature_Higgs_L4.at(2).at(2).at(0).at(0) = 2*lam;
  Curvature_Higgs_L4.at(2).at(2).at(1).at(1) = 2*lam;
  Curvature_Higgs_L4.at(2).at(2).at(2).at(2) = 6*lam;
  Curvature_Higgs_L4.at(2).at(2).at(3).at(3) = 2*lam;
  Curvature_Higgs_L4.at(2).at(2).at(4).at(4) = a2;
  Curvature_Higgs_L4.at(2).at(3).at(2).at(3) = 2*lam;
  Curvature_Higgs_L4.at(2).at(3).at(3).at(2) = 2*lam;
  Curvature_Higgs_L4.at(2).at(4).at(2).at(4) = a2;
  Curvature_Higgs_L4.at(2).at(4).at(4).at(2) = a2;
  Curvature_Higgs_L4.at(3).at(0).at(0).at(3) = 2*lam;
  Curvature_Higgs_L4.at(3).at(0).at(3).at(0) = 2*lam;
  Curvature_Higgs_L4.at(3).at(1).at(1).at(3) = 2*lam;
  Curvature_Higgs_L4.at(3).at(1).at(3).at(1) = 2*lam;
  Curvature_Higgs_L4.at(3).at(2).at(2).at(3) = 2*lam;
  Curvature_Higgs_L4.at(3).at(2).at(3).at(2) = 2*lam;
  Curvature_Higgs_L4.at(3).at(3).at(0).at(0) = 2*lam;
  Curvature_Higgs_L4.at(3).at(3).at(1).at(1) = 2*lam;
  Curvature_Higgs_L4.at(3).at(3).at(2).at(2) = 2*lam;
  Curvature_Higgs_L4.at(3).at(3).at(3).at(3) = 6*lam;
  Curvature_Higgs_L4.at(3).at(3).at(4).at(4) = a2;
  Curvature_Higgs_L4.at(3).at(4).at(3).at(4) = a2;
  Curvature_Higgs_L4.at(3).at(4).at(4).at(3) = a2;
  Curvature_Higgs_L4.at(4).at(0).at(0).at(4) = a2;
  Curvature_Higgs_L4.at(4).at(0).at(4).at(0) = a2;
  Curvature_Higgs_L4.at(4).at(1).at(1).at(4) = a2;
  Curvature_Higgs_L4.at(4).at(1).at(4).at(1) = a2;
  Curvature_Higgs_L4.at(4).at(2).at(2).at(4) = a2;
  Curvature_Higgs_L4.at(4).at(2).at(4).at(2) = a2;
  Curvature_Higgs_L4.at(4).at(3).at(3).at(4) = a2;
  Curvature_Higgs_L4.at(4).at(3).at(4).at(3) = a2;
  Curvature_Higgs_L4.at(4).at(4).at(0).at(0) = a2;
  Curvature_Higgs_L4.at(4).at(4).at(1).at(1) = a2;
  Curvature_Higgs_L4.at(4).at(4).at(2).at(2) = a2;
  Curvature_Higgs_L4.at(4).at(4).at(3).at(3) = a2;
  Curvature_Higgs_L4.at(4).at(4).at(4).at(4) = 6*b4;

  Curvature_Gauge_G2H2.at(0).at(0).at(0).at(0) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(0).at(0).at(1).at(1) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(0).at(0).at(2).at(2) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(0).at(0).at(3).at(3) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(0).at(3).at(0).at(3) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(0).at(3).at(1).at(2) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(0).at(3).at(2).at(1) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(0).at(3).at(3).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(1).at(1).at(0).at(0) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(1).at(1).at(1).at(1) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(1).at(1).at(2).at(2) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(1).at(1).at(3).at(3) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(1).at(3).at(0).at(2) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(1).at(3).at(1).at(3) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(1).at(3).at(2).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(1).at(3).at(3).at(1) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(2).at(2).at(0).at(0) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(2).at(2).at(1).at(1) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(2).at(2).at(2).at(2) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(2).at(2).at(3).at(3) = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2.at(2).at(3).at(0).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(2).at(3).at(1).at(1) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(2).at(3).at(2).at(2) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(2).at(3).at(3).at(3) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(0).at(0).at(3) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(0).at(1).at(2) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(0).at(2).at(1) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(0).at(3).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(1).at(0).at(2) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(1).at(1).at(3) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(1).at(2).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(1).at(3).at(1) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(2).at(0).at(0) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(2).at(1).at(1) = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(2).at(2).at(2) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(2).at(3).at(3) = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2.at(3).at(3).at(0).at(0) = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2.at(3).at(3).at(1).at(1) = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2.at(3).at(3).at(2).at(2) = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2.at(3).at(3).at(3).at(3) = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);

  Curvature_Lepton_F2H1.at(0).at(3).at(0) = SMConstants.C_MassElectron/vH;
  Curvature_Lepton_F2H1.at(0).at(3).at(1) = SMConstants.C_MassElectron*II/vH;
  Curvature_Lepton_F2H1.at(1).at(4).at(0) = SMConstants.C_MassMu/vH;
  Curvature_Lepton_F2H1.at(1).at(4).at(1) = SMConstants.C_MassMu*II/vH;
  Curvature_Lepton_F2H1.at(2).at(5).at(0) = SMConstants.C_MassTau/vH;
  Curvature_Lepton_F2H1.at(2).at(5).at(1) = SMConstants.C_MassTau*II/vH;
  Curvature_Lepton_F2H1.at(3).at(0).at(0) = SMConstants.C_MassElectron/vH;
  Curvature_Lepton_F2H1.at(3).at(0).at(1) = SMConstants.C_MassElectron*II/vH;
  Curvature_Lepton_F2H1.at(3).at(6).at(2) = SMConstants.C_MassElectron*II/vH;
  Curvature_Lepton_F2H1.at(3).at(6).at(3) = SMConstants.C_MassElectron/vH;
  Curvature_Lepton_F2H1.at(4).at(1).at(0) = SMConstants.C_MassMu/vH;
  Curvature_Lepton_F2H1.at(4).at(1).at(1) = SMConstants.C_MassMu*II/vH;
  Curvature_Lepton_F2H1.at(4).at(7).at(2) = SMConstants.C_MassMu*II/vH;
  Curvature_Lepton_F2H1.at(4).at(7).at(3) = SMConstants.C_MassMu/vH;
  Curvature_Lepton_F2H1.at(5).at(2).at(0) = SMConstants.C_MassTau/vH;
  Curvature_Lepton_F2H1.at(5).at(2).at(1) = SMConstants.C_MassTau*II/vH;
  Curvature_Lepton_F2H1.at(5).at(8).at(2) = SMConstants.C_MassTau*II/vH;
  Curvature_Lepton_F2H1.at(5).at(8).at(3) = SMConstants.C_MassTau/vH;
  Curvature_Lepton_F2H1.at(6).at(3).at(2) = SMConstants.C_MassElectron*II/vH;
  Curvature_Lepton_F2H1.at(6).at(3).at(3) = SMConstants.C_MassElectron/vH;
  Curvature_Lepton_F2H1.at(7).at(4).at(2) = SMConstants.C_MassMu*II/vH;
  Curvature_Lepton_F2H1.at(7).at(4).at(3) = SMConstants.C_MassMu/vH;
  Curvature_Lepton_F2H1.at(8).at(5).at(2) = SMConstants.C_MassTau*II/vH;
  Curvature_Lepton_F2H1.at(8).at(5).at(3) = SMConstants.C_MassTau/vH;

  std::complex<double> Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb;
  Vud = SMConstants.C_Vud;
  Vus = SMConstants.C_Vus;
  Vub = SMConstants.C_Vub;
  Vcd = SMConstants.C_Vcd;
  Vcs = SMConstants.C_Vcs;
  Vcb = SMConstants.C_Vcb;
  Vtd = SMConstants.C_Vtd;
  Vts = SMConstants.C_Vts;
  Vtb = SMConstants.C_Vtb;

  Curvature_Quark_F2H1.at(0).at(6).at(2) = -SMConstants.C_MassUp*II/vH;
  Curvature_Quark_F2H1.at(0).at(6).at(3) = SMConstants.C_MassUp/vH;
  Curvature_Quark_F2H1.at(0).at(9).at(0) = -SMConstants.C_MassUp*conj(Vud)/vH;
  Curvature_Quark_F2H1.at(0).at(9).at(1) = SMConstants.C_MassUp*II*conj(Vud)/vH;
  Curvature_Quark_F2H1.at(0).at(10).at(0) = -SMConstants.C_MassUp*conj(Vus)/vH;
  Curvature_Quark_F2H1.at(0).at(10).at(1) = SMConstants.C_MassUp*II*conj(Vus)/vH;
  Curvature_Quark_F2H1.at(0).at(11).at(0) = -SMConstants.C_MassUp*conj(Vub)/vH;
  Curvature_Quark_F2H1.at(0).at(11).at(1) = SMConstants.C_MassUp*II*conj(Vub)/vH;
  Curvature_Quark_F2H1.at(1).at(7).at(2) = -SMConstants.C_MassCharm*II/vH;
  Curvature_Quark_F2H1.at(1).at(7).at(3) = SMConstants.C_MassCharm/vH;
  Curvature_Quark_F2H1.at(1).at(9).at(0) = -SMConstants.C_MassCharm*conj(Vcd)/vH;
  Curvature_Quark_F2H1.at(1).at(9).at(1) = SMConstants.C_MassCharm*II*conj(Vcd)/vH;
  Curvature_Quark_F2H1.at(1).at(10).at(0) = -SMConstants.C_MassCharm*conj(Vcs)/vH;
  Curvature_Quark_F2H1.at(1).at(10).at(1) = SMConstants.C_MassCharm*II*conj(Vcs)/vH;
  Curvature_Quark_F2H1.at(1).at(11).at(0) = -SMConstants.C_MassCharm*conj(Vcb)/vH;
  Curvature_Quark_F2H1.at(1).at(11).at(1) = SMConstants.C_MassCharm*II*conj(Vcb)/vH;
  Curvature_Quark_F2H1.at(2).at(8).at(2) = -SMConstants.C_MassTop*II/vH;
  Curvature_Quark_F2H1.at(2).at(8).at(3) = SMConstants.C_MassTop/vH;
  Curvature_Quark_F2H1.at(2).at(9).at(0) = -SMConstants.C_MassTop*conj(Vtd)/vH;
  Curvature_Quark_F2H1.at(2).at(9).at(1) = SMConstants.C_MassTop*II*conj(Vtd)/vH;
  Curvature_Quark_F2H1.at(2).at(10).at(0) = -SMConstants.C_MassTop*conj(Vts)/vH;
  Curvature_Quark_F2H1.at(2).at(10).at(1) = SMConstants.C_MassTop*II*conj(Vts)/vH;
  Curvature_Quark_F2H1.at(2).at(11).at(0) = -SMConstants.C_MassTop*conj(Vtb)/vH;
  Curvature_Quark_F2H1.at(2).at(11).at(1) = SMConstants.C_MassTop*II*conj(Vtb)/vH;
  Curvature_Quark_F2H1.at(3).at(6).at(0) = SMConstants.C_MassDown*Vud/vH;
  Curvature_Quark_F2H1.at(3).at(6).at(1) = SMConstants.C_MassDown*II*Vud/vH;
  Curvature_Quark_F2H1.at(3).at(7).at(0) = SMConstants.C_MassDown*Vcd/vH;
  Curvature_Quark_F2H1.at(3).at(7).at(1) = SMConstants.C_MassDown*II*Vcd/vH;
  Curvature_Quark_F2H1.at(3).at(8).at(0) = SMConstants.C_MassDown*Vtd/vH;
  Curvature_Quark_F2H1.at(3).at(8).at(1) = SMConstants.C_MassDown*II*Vtd/vH;
  Curvature_Quark_F2H1.at(3).at(9).at(2) = SMConstants.C_MassDown*II/vH;
  Curvature_Quark_F2H1.at(3).at(9).at(3) = SMConstants.C_MassDown/vH;
  Curvature_Quark_F2H1.at(4).at(6).at(0) = SMConstants.C_MassStrange*Vus/vH;
  Curvature_Quark_F2H1.at(4).at(6).at(1) = SMConstants.C_MassStrange*II*Vus/vH;
  Curvature_Quark_F2H1.at(4).at(7).at(0) = SMConstants.C_MassStrange*Vcs/vH;
  Curvature_Quark_F2H1.at(4).at(7).at(1) = SMConstants.C_MassStrange*II*Vcs/vH;
  Curvature_Quark_F2H1.at(4).at(8).at(0) = SMConstants.C_MassStrange*Vts/vH;
  Curvature_Quark_F2H1.at(4).at(8).at(1) = SMConstants.C_MassStrange*II*Vts/vH;
  Curvature_Quark_F2H1.at(4).at(10).at(2) = SMConstants.C_MassStrange*II/vH;
  Curvature_Quark_F2H1.at(4).at(10).at(3) = SMConstants.C_MassStrange/vH;
  Curvature_Quark_F2H1.at(5).at(6).at(0) = SMConstants.C_MassBottom*Vub/vH;
  Curvature_Quark_F2H1.at(5).at(6).at(1) = SMConstants.C_MassBottom*II*Vub/vH;
  Curvature_Quark_F2H1.at(5).at(7).at(0) = SMConstants.C_MassBottom*Vcb/vH;
  Curvature_Quark_F2H1.at(5).at(7).at(1) = SMConstants.C_MassBottom*II*Vcb/vH;
  Curvature_Quark_F2H1.at(5).at(8).at(0) = SMConstants.C_MassBottom*Vtb/vH;
  Curvature_Quark_F2H1.at(5).at(8).at(1) = SMConstants.C_MassBottom*II*Vtb/vH;
  Curvature_Quark_F2H1.at(5).at(11).at(2) = SMConstants.C_MassBottom*II/vH;
  Curvature_Quark_F2H1.at(5).at(11).at(3) = SMConstants.C_MassBottom/vH;
  Curvature_Quark_F2H1.at(6).at(0).at(2) = -SMConstants.C_MassUp*II/vH;
  Curvature_Quark_F2H1.at(6).at(0).at(3) = SMConstants.C_MassUp/vH;
  Curvature_Quark_F2H1.at(6).at(3).at(0) = SMConstants.C_MassDown*Vud/vH;
  Curvature_Quark_F2H1.at(6).at(3).at(1) = SMConstants.C_MassDown*II*Vud/vH;
  Curvature_Quark_F2H1.at(6).at(4).at(0) = SMConstants.C_MassStrange*Vus/vH;
  Curvature_Quark_F2H1.at(6).at(4).at(1) = SMConstants.C_MassStrange*II*Vus/vH;
  Curvature_Quark_F2H1.at(6).at(5).at(0) = SMConstants.C_MassBottom*Vub/vH;
  Curvature_Quark_F2H1.at(6).at(5).at(1) = SMConstants.C_MassBottom*II*Vub/vH;
  Curvature_Quark_F2H1.at(7).at(1).at(2) = -SMConstants.C_MassCharm*II/vH;
  Curvature_Quark_F2H1.at(7).at(1).at(3) = SMConstants.C_MassCharm/vH;
  Curvature_Quark_F2H1.at(7).at(3).at(0) = SMConstants.C_MassDown*Vcd/vH;
  Curvature_Quark_F2H1.at(7).at(3).at(1) = SMConstants.C_MassDown*II*Vcd/vH;
  Curvature_Quark_F2H1.at(7).at(4).at(0) = SMConstants.C_MassStrange*Vcs/vH;
  Curvature_Quark_F2H1.at(7).at(4).at(1) = SMConstants.C_MassStrange*II*Vcs/vH;
  Curvature_Quark_F2H1.at(7).at(5).at(0) = SMConstants.C_MassBottom*Vcb/vH;
  Curvature_Quark_F2H1.at(7).at(5).at(1) = SMConstants.C_MassBottom*II*Vcb/vH;
  Curvature_Quark_F2H1.at(8).at(2).at(2) = -SMConstants.C_MassTop*II/vH;
  Curvature_Quark_F2H1.at(8).at(2).at(3) = SMConstants.C_MassTop/vH;
  Curvature_Quark_F2H1.at(8).at(3).at(0) = SMConstants.C_MassDown*Vtd/vH;
  Curvature_Quark_F2H1.at(8).at(3).at(1) = SMConstants.C_MassDown*II*Vtd/vH;
  Curvature_Quark_F2H1.at(8).at(4).at(0) = SMConstants.C_MassStrange*Vts/vH;
  Curvature_Quark_F2H1.at(8).at(4).at(1) = SMConstants.C_MassStrange*II*Vts/vH;
  Curvature_Quark_F2H1.at(8).at(5).at(0) = SMConstants.C_MassBottom*Vtb/vH;
  Curvature_Quark_F2H1.at(8).at(5).at(1) = SMConstants.C_MassBottom*II*Vtb/vH;
  Curvature_Quark_F2H1.at(9).at(0).at(0) = -SMConstants.C_MassUp*conj(Vud)/vH;
  Curvature_Quark_F2H1.at(9).at(0).at(1) = SMConstants.C_MassUp*II*conj(Vud)/vH;
  Curvature_Quark_F2H1.at(9).at(1).at(0) = -SMConstants.C_MassCharm*conj(Vcd)/vH;
  Curvature_Quark_F2H1.at(9).at(1).at(1) = SMConstants.C_MassCharm*II*conj(Vcd)/vH;
  Curvature_Quark_F2H1.at(9).at(2).at(0) = -SMConstants.C_MassTop*conj(Vtd)/vH;
  Curvature_Quark_F2H1.at(9).at(2).at(1) = SMConstants.C_MassTop*II*conj(Vtd)/vH;
  Curvature_Quark_F2H1.at(9).at(3).at(2) = SMConstants.C_MassDown*II/vH;
  Curvature_Quark_F2H1.at(9).at(3).at(3) = SMConstants.C_MassDown/vH;
  Curvature_Quark_F2H1.at(10).at(0).at(0) = -SMConstants.C_MassUp*conj(Vus)/vH;
  Curvature_Quark_F2H1.at(10).at(0).at(1) = SMConstants.C_MassUp*II*conj(Vus)/vH;
  Curvature_Quark_F2H1.at(10).at(1).at(0) = -SMConstants.C_MassCharm*conj(Vcs)/vH;
  Curvature_Quark_F2H1.at(10).at(1).at(1) = SMConstants.C_MassCharm*II*conj(Vcs)/vH;
  Curvature_Quark_F2H1.at(10).at(2).at(0) = -SMConstants.C_MassTop*conj(Vts)/vH;
  Curvature_Quark_F2H1.at(10).at(2).at(1) = SMConstants.C_MassTop*II*conj(Vts)/vH;
  Curvature_Quark_F2H1.at(10).at(4).at(2) = SMConstants.C_MassStrange*II/vH;
  Curvature_Quark_F2H1.at(10).at(4).at(3) = SMConstants.C_MassStrange/vH;
  Curvature_Quark_F2H1.at(11).at(0).at(0) = -SMConstants.C_MassUp*conj(Vub)/vH;
  Curvature_Quark_F2H1.at(11).at(0).at(1) = SMConstants.C_MassUp*II*conj(Vub)/vH;
  Curvature_Quark_F2H1.at(11).at(1).at(0) = -SMConstants.C_MassCharm*conj(Vcb)/vH;
  Curvature_Quark_F2H1.at(11).at(1).at(1) = SMConstants.C_MassCharm*II*conj(Vcb)/vH;
  Curvature_Quark_F2H1.at(11).at(2).at(0) = -SMConstants.C_MassTop*conj(Vtb)/vH;
  Curvature_Quark_F2H1.at(11).at(2).at(1) = SMConstants.C_MassTop*II*conj(Vtb)/vH;
  Curvature_Quark_F2H1.at(11).at(5).at(2) = SMConstants.C_MassBottom*II/vH;
  Curvature_Quark_F2H1.at(11).at(5).at(3) = SMConstants.C_MassBottom/vH;

}

bool Class_RxSM::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_RxSM::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_RxSM::VTreeSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  return res;
}

double Class_RxSM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  return res;
}

void Class_RxSM::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
