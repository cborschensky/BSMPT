// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later


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
  Model =
      ModelID::ModelIDs::RXSM; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 3;               // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 9; // number of parameters in the tree-Level Lagrangian
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
  labels.push_back("dmu2");
  labels.push_back("dMS");
  labels.push_back("dKapS");
  labels.push_back("dKapSH");
  labels.push_back("dLam");
  labels.push_back("dLamS");
  labels.push_back("dLamSH");
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
  labels.push_back("v_c"); // Label for the critical vev
  labels.push_back(
      "v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
  // out += "Your VEV order"; // Now you have to put the label for your vevs
  labels.push_back("omegav");
  labels.push_back("omegavs");
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
  particles[3] = "h";
  particles[4] = "H";

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
  labels.push_back("omegav");
  labels.push_back("omegavs");
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

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 9; k++)
  {
    ss >> tmp;
    if (k == 1)
      par[0] = tmp; // v
    else if (k == 2)
      par[1] = tmp; // vs
    else if (k == 3)
      par[2] = tmp; // mu2
    else if (k == 4)
      par[3] = tmp; // MS
    else if (k == 5)
      par[4] = tmp; // KapS
    else if (k == 6)
      par[5] = tmp; // KapSH
    else if (k == 7)
      par[6] = tmp; // Lam
    else if (k == 8)
      par[7] = tmp; // LamS
    else if (k == 9)
      par[8] = tmp; // LamSH;
  }

  set_gen(par); // This you have to call so that everything will be set
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_RxSM::set_gen(const std::vector<double> &par)
{
  vh     = par[0]; // Class member is set accordingly to the input parameters
  vs    = par[1]; // Class member is set accordingly to the input parameters
  //mu2i   = par[2];
  //MSi    = par[3];
  KapS  = par[4];
  KapSH = par[5];
  Lam   = par[6];
  LamS  = par[7];
  LamSH = par[8];
  mu2   = -KapSH*vs - 1.0/2.0*Lam*std::pow(vh, 2) - 1.0/2.0*LamSH*std::pow(vs, 2);
  MS     = -KapS*vs - 1.0/2.0*KapSH*std::pow(vh, 2)/vs - 2*LamS*std::pow(vs, 2) - 1.0/2.0*LamSH*std::pow(vh, 2);
  
   
  //scale = sqrt(-KapSH*246); // Renormalisation scale is set to the SM VEV
  scale = 246.22; // Renormalisation scale is set to the SM VEV
  
  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = vh;
  vevTreeMin[1] = vs;
  
  
  vevTree       = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_RxSM::set_CT_Pot_Par(const std::vector<double> &par)
{

  dmu2      = par[0];
  dMS       = par[1];
  dKapS     = par[2];
  dKapSH    = par[3];
  dLam      = par[4];
  dLamS     = par[5];
  dLamSH    = par[6];
  dT1       = par[7];
  dT2       = par[8];
  dT3       = par[9];
  dT4       = par[10];
  dT5       = par[11]*0;
  
  
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L2[0][0] = dmu2;
  Curvature_Higgs_CT_L2[1][1] = dmu2;
  Curvature_Higgs_CT_L2[2][2] = dmu2;
  Curvature_Higgs_CT_L2[3][3] = dmu2;
  Curvature_Higgs_CT_L2[4][4] = dMS;
  Curvature_Higgs_CT_L3[0][0][4] = dKapSH;
  Curvature_Higgs_CT_L3[0][4][0] = dKapSH;
  Curvature_Higgs_CT_L3[1][1][4] = dKapSH;
  Curvature_Higgs_CT_L3[1][4][1] = dKapSH;
  Curvature_Higgs_CT_L3[2][2][4] = dKapSH;
  Curvature_Higgs_CT_L3[2][4][2] = dKapSH;
  Curvature_Higgs_CT_L3[3][3][4] = dKapSH;
  Curvature_Higgs_CT_L3[3][4][3] = dKapSH;
  Curvature_Higgs_CT_L3[4][0][0] = dKapSH;
  Curvature_Higgs_CT_L3[4][1][1] = dKapSH;
  Curvature_Higgs_CT_L3[4][2][2] = dKapSH;
  Curvature_Higgs_CT_L3[4][3][3] = dKapSH;
  Curvature_Higgs_CT_L3[4][4][4] = 2*dKapS;
  Curvature_Higgs_CT_L4[0][0][0][0] = 3*dLam;
  Curvature_Higgs_CT_L4[0][0][1][1] = dLam;
  Curvature_Higgs_CT_L4[0][0][2][2] = dLam;
  Curvature_Higgs_CT_L4[0][0][3][3] = dLam;
  Curvature_Higgs_CT_L4[0][0][4][4] = dLamSH;
  Curvature_Higgs_CT_L4[0][1][0][1] = dLam;
  Curvature_Higgs_CT_L4[0][1][1][0] = dLam;
  Curvature_Higgs_CT_L4[0][2][0][2] = dLam;
  Curvature_Higgs_CT_L4[0][2][2][0] = dLam;
  Curvature_Higgs_CT_L4[0][3][0][3] = dLam;
  Curvature_Higgs_CT_L4[0][3][3][0] = dLam;
  Curvature_Higgs_CT_L4[0][4][0][4] = dLamSH;
  Curvature_Higgs_CT_L4[0][4][4][0] = dLamSH;
  Curvature_Higgs_CT_L4[1][0][0][1] = dLam;
  Curvature_Higgs_CT_L4[1][0][1][0] = dLam;
  Curvature_Higgs_CT_L4[1][1][0][0] = dLam;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3*dLam;
  Curvature_Higgs_CT_L4[1][1][2][2] = dLam;
  Curvature_Higgs_CT_L4[1][1][3][3] = dLam;
  Curvature_Higgs_CT_L4[1][1][4][4] = dLamSH;
  Curvature_Higgs_CT_L4[1][2][1][2] = dLam;
  Curvature_Higgs_CT_L4[1][2][2][1] = dLam;
  Curvature_Higgs_CT_L4[1][3][1][3] = dLam;
  Curvature_Higgs_CT_L4[1][3][3][1] = dLam;
  Curvature_Higgs_CT_L4[1][4][1][4] = dLamSH;
  Curvature_Higgs_CT_L4[1][4][4][1] = dLamSH;
  Curvature_Higgs_CT_L4[2][0][0][2] = dLam;
  Curvature_Higgs_CT_L4[2][0][2][0] = dLam;
  Curvature_Higgs_CT_L4[2][1][1][2] = dLam;
  Curvature_Higgs_CT_L4[2][1][2][1] = dLam;
  Curvature_Higgs_CT_L4[2][2][0][0] = dLam;
  Curvature_Higgs_CT_L4[2][2][1][1] = dLam;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3*dLam;
  Curvature_Higgs_CT_L4[2][2][3][3] = dLam;
  Curvature_Higgs_CT_L4[2][2][4][4] = dLamSH;
  Curvature_Higgs_CT_L4[2][3][2][3] = dLam;
  Curvature_Higgs_CT_L4[2][3][3][2] = dLam;
  Curvature_Higgs_CT_L4[2][4][2][4] = dLamSH;
  Curvature_Higgs_CT_L4[2][4][4][2] = dLamSH;
  Curvature_Higgs_CT_L4[3][0][0][3] = dLam;
  Curvature_Higgs_CT_L4[3][0][3][0] = dLam;
  Curvature_Higgs_CT_L4[3][1][1][3] = dLam;
  Curvature_Higgs_CT_L4[3][1][3][1] = dLam;
  Curvature_Higgs_CT_L4[3][2][2][3] = dLam;
  Curvature_Higgs_CT_L4[3][2][3][2] = dLam;
  Curvature_Higgs_CT_L4[3][3][0][0] = dLam;
  Curvature_Higgs_CT_L4[3][3][1][1] = dLam;
  Curvature_Higgs_CT_L4[3][3][2][2] = dLam;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3*dLam;
  Curvature_Higgs_CT_L4[3][3][4][4] = dLamSH;
  Curvature_Higgs_CT_L4[3][4][3][4] = dLamSH;
  Curvature_Higgs_CT_L4[3][4][4][3] = dLamSH;
  Curvature_Higgs_CT_L4[4][0][0][4] = dLamSH;
  Curvature_Higgs_CT_L4[4][0][4][0] = dLamSH;
  Curvature_Higgs_CT_L4[4][1][1][4] = dLamSH;
  Curvature_Higgs_CT_L4[4][1][4][1] = dLamSH;
  Curvature_Higgs_CT_L4[4][2][2][4] = dLamSH;
  Curvature_Higgs_CT_L4[4][2][4][2] = dLamSH;
  Curvature_Higgs_CT_L4[4][3][3][4] = dLamSH;
  Curvature_Higgs_CT_L4[4][3][4][3] = dLamSH;
  Curvature_Higgs_CT_L4[4][4][0][0] = dLamSH;
  Curvature_Higgs_CT_L4[4][4][1][1] = dLamSH;
  Curvature_Higgs_CT_L4[4][4][2][2] = dLamSH;
  Curvature_Higgs_CT_L4[4][4][3][3] = dLamSH;
  Curvature_Higgs_CT_L4[4][4][4][4] = 12*dLamS;

  
}

/**
 * console output of all Parameters
 */
void Class_RxSM::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "mu2 = " << mu2 << std::endl
     << "MS = " << MS << std::endl
     << "KapS = " << KapS << std::endl
     << "KapSH = " << KapSH << std::endl
     << "Lam = " << Lam << std::endl
     << "LamS = " << LamS << std::endl
     << "LamSH = " << LamSH << std::endl
     << "vh = " << vh << std::endl 
     << "vs = " << vs << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
  ss << "dmu2 = " << dmu2 << std::endl
     << "dMS = " << dMS << std::endl
     << "dKapS = " << dKapS << std::endl
     << "dKapSH = " << dKapSH << std::endl
     << "dLam = " << dLam << std::endl
     << "dLamS = " << dLamS << std::endl
     << "dLamSH = " << dLamSH << std::endl
     << "dT1 = " << dT1 << std::endl
     << "dT2 = " << dT2 << std::endl
     << "dT3 = " << dT3 << std::endl
     << "dT4 = " << dT4 << std::endl 
     << "dT5 = " << dT5 << std::endl;


  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_RxSM::calc_CT() const
{

  std::vector<double> parCT;
  parCT.resize(nParCT);

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


  // Here you have to use your formulae for the counterterm scheme
  //(-1.0/2.0*HesseWeinberg(2,2)*(vh*vh)/(vs*vs*vs) + (1.0/2.0)*HesseWeinberg(2,4)*vh/(vs*vs) + (3.0/2.0)*HesseWeinberg(3,3)*(vh*vh)/(vs*vs*vs) + HesseWeinberg(4,4)/vs - 3*NablaWeinberg(4)/(vs*vs) - 3*dT5/(vs*vs))
  //(HesseWeinberg(2,2)/vs + HesseWeinberg(2,4)/vh - 3*HesseWeinberg(3,3)/vs)
  //
  double ZeroCut = 1e-5;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(NablaWeinberg[i]) <= ZeroCut) NablaWeinberg[i] = 0;
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(HesseWeinberg(i, j)) <= ZeroCut) HesseWeinberg(i, j) = 0;
    }
  }
  
  
  
  parCT[0] = (0); //dmu2
  parCT[1] = (0); //dMS
  parCT[2] = ((3.0/2.0)*HesseWeinberg(2,2)*std::pow(vh, 2)/std::pow(vs, 3) - 1.0/2.0*HesseWeinberg(3,3)*std::pow(vh, 2)/std::pow(vs, 3) + (1.0/2.0)*HesseWeinberg(3,4)*vh/std::pow(vs, 2) + HesseWeinberg(4,4)/vs - 3*NablaWeinberg(4)/std::pow(vs, 2)); //dKapS
  parCT[3] = (-3*HesseWeinberg(2,2)/vs + HesseWeinberg(3,3)/vs + HesseWeinberg(3,4)/vh); //dKapSH
  parCT[4] = (HesseWeinberg(2,2)/std::pow(vh, 2) - HesseWeinberg(3,3)/std::pow(vh, 2)); //dLam
  parCT[5] = (-3.0/4.0*HesseWeinberg(2,2)*std::pow(vh, 2)/std::pow(vs, 4) + (1.0/4.0)*HesseWeinberg(3,3)*std::pow(vh, 2)/std::pow(vs, 4) - 1.0/2.0*HesseWeinberg(4,4)/std::pow(vs, 2) + NablaWeinberg(4)/std::pow(vs, 3)); //dLamS
  parCT[6] = (3*HesseWeinberg(2,2)/std::pow(vs, 2) - HesseWeinberg(3,3)/std::pow(vs, 2) - 2*HesseWeinberg(3,4)/(vh*vs)); //dLamSH
  
  parCT[7] = (-NablaWeinberg(0)); //dT1
  parCT[8] = (-NablaWeinberg(1)); //dT2
  parCT[9] = (-NablaWeinberg(2)); //dT3
  parCT[10] = (HesseWeinberg(2,2)*vh - NablaWeinberg(3)); //dT4
  parCT[11] = (0); //dT5  
  
  //parCT[0] = (0); //dmu2
  //parCT[1] = (0); //dMS
  //parCT[2] = ((3.0/2.0)*HesseWeinberg(2,2)*std::pow(vh, 2)/std::pow(vs, 3) - 1.0/2.0*HesseWeinberg(3,3)*std::pow(vh, 2)/std::pow(vs, 3) + (1.0/2.0)*HesseWeinberg(3,4)*vh/std::pow(vs, 2) + HesseWeinberg(4,4)/vs - 3*NablaWeinberg(4)/std::pow(vs, 2) - 3*dT5/std::pow(vs, 2)); //dKapS
  //parCT[3] = (-3*HesseWeinberg(2,2)/vs + HesseWeinberg(3,3)/vs + HesseWeinberg(3,4)/vh); //dKapSH
  //parCT[4] = (HesseWeinberg(2,2)/std::pow(vh, 2) - HesseWeinberg(3,3)/std::pow(vh, 2)); //dLam
  //parCT[5] = (-3.0/4.0*HesseWeinberg(2,2)*std::pow(vh, 2)/std::pow(vs, 4) + (1.0/4.0)*HesseWeinberg(3,3)*std::pow(vh, 2)/std::pow(vs, 4) - 1.0/2.0*HesseWeinberg(4,4)/std::pow(vs, 2) + NablaWeinberg(4)/std::pow(vs, 3) + dT5/std::pow(vs, 3)); //dLamS
  //parCT[6] = (3*HesseWeinberg(2,2)/std::pow(vs, 2) - HesseWeinberg(3,3)/std::pow(vs, 2) - 2*HesseWeinberg(3,4)/(vh*vs)); //dLamSH
  
  //parCT[7] = (-NablaWeinberg(0)); //dT1
  //parCT[8] = (-NablaWeinberg(1)); //dT2
  //parCT[9] = (-NablaWeinberg(2)); //dT3
  //parCT[10] = (HesseWeinberg(2,2)*vh - NablaWeinberg(3)); //dT4
  //parCT[11] = (dT5); //dT5  
  for (std::size_t i = 0; i < nParCT; i++)
  {
    if (std::abs(parCT[i]) <= ZeroCut) parCT[i] = 0;
  }


  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_RxSM::AdjustRotationMatrix()
{
}

void Class_RxSM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
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

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
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
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = mu2;
  Curvature_Higgs_L2[1][1] = mu2;
  Curvature_Higgs_L2[2][2] = mu2;
  Curvature_Higgs_L2[3][3] = mu2;
  Curvature_Higgs_L2[4][4] = MS;
  Curvature_Higgs_L3[0][0][4] = KapSH;
  Curvature_Higgs_L3[0][4][0] = KapSH;
  Curvature_Higgs_L3[1][1][4] = KapSH;
  Curvature_Higgs_L3[1][4][1] = KapSH;
  Curvature_Higgs_L3[2][2][4] = KapSH;
  Curvature_Higgs_L3[2][4][2] = KapSH;
  Curvature_Higgs_L3[3][3][4] = KapSH;
  Curvature_Higgs_L3[3][4][3] = KapSH;
  Curvature_Higgs_L3[4][0][0] = KapSH;
  Curvature_Higgs_L3[4][1][1] = KapSH;
  Curvature_Higgs_L3[4][2][2] = KapSH;
  Curvature_Higgs_L3[4][3][3] = KapSH;
  Curvature_Higgs_L3[4][4][4] = 0.2e1*KapS;
  Curvature_Higgs_L4[0][0][0][0] = 0.3e1*Lam;
  Curvature_Higgs_L4[0][0][1][1] = Lam;
  Curvature_Higgs_L4[0][0][2][2] = Lam;
  Curvature_Higgs_L4[0][0][3][3] = Lam;
  Curvature_Higgs_L4[0][0][4][4] = LamSH;
  Curvature_Higgs_L4[0][1][0][1] = Lam;
  Curvature_Higgs_L4[0][1][1][0] = Lam;
  Curvature_Higgs_L4[0][2][0][2] = Lam;
  Curvature_Higgs_L4[0][2][2][0] = Lam;
  Curvature_Higgs_L4[0][3][0][3] = Lam;
  Curvature_Higgs_L4[0][3][3][0] = Lam;
  Curvature_Higgs_L4[0][4][0][4] = LamSH;
  Curvature_Higgs_L4[0][4][4][0] = LamSH;
  Curvature_Higgs_L4[1][0][0][1] = Lam;
  Curvature_Higgs_L4[1][0][1][0] = Lam;
  Curvature_Higgs_L4[1][1][0][0] = Lam;
  Curvature_Higgs_L4[1][1][1][1] = 0.3e1*Lam;
  Curvature_Higgs_L4[1][1][2][2] = Lam;
  Curvature_Higgs_L4[1][1][3][3] = Lam;
  Curvature_Higgs_L4[1][1][4][4] = LamSH;
  Curvature_Higgs_L4[1][2][1][2] = Lam;
  Curvature_Higgs_L4[1][2][2][1] = Lam;
  Curvature_Higgs_L4[1][3][1][3] = Lam;
  Curvature_Higgs_L4[1][3][3][1] = Lam;
  Curvature_Higgs_L4[1][4][1][4] = LamSH;
  Curvature_Higgs_L4[1][4][4][1] = LamSH;
  Curvature_Higgs_L4[2][0][0][2] = Lam;
  Curvature_Higgs_L4[2][0][2][0] = Lam;
  Curvature_Higgs_L4[2][1][1][2] = Lam;
  Curvature_Higgs_L4[2][1][2][1] = Lam;
  Curvature_Higgs_L4[2][2][0][0] = Lam;
  Curvature_Higgs_L4[2][2][1][1] = Lam;
  Curvature_Higgs_L4[2][2][2][2] = 0.3e1*Lam;
  Curvature_Higgs_L4[2][2][3][3] = Lam;
  Curvature_Higgs_L4[2][2][4][4] = LamSH;
  Curvature_Higgs_L4[2][3][2][3] = Lam;
  Curvature_Higgs_L4[2][3][3][2] = Lam;
  Curvature_Higgs_L4[2][4][2][4] = LamSH;
  Curvature_Higgs_L4[2][4][4][2] = LamSH;
  Curvature_Higgs_L4[3][0][0][3] = Lam;
  Curvature_Higgs_L4[3][0][3][0] = Lam;
  Curvature_Higgs_L4[3][1][1][3] = Lam;
  Curvature_Higgs_L4[3][1][3][1] = Lam;
  Curvature_Higgs_L4[3][2][2][3] = Lam;
  Curvature_Higgs_L4[3][2][3][2] = Lam;
  Curvature_Higgs_L4[3][3][0][0] = Lam;
  Curvature_Higgs_L4[3][3][1][1] = Lam;
  Curvature_Higgs_L4[3][3][2][2] = Lam;
  Curvature_Higgs_L4[3][3][3][3] = 0.3e1*Lam;
  Curvature_Higgs_L4[3][3][4][4] = LamSH;
  Curvature_Higgs_L4[3][4][3][4] = LamSH;
  Curvature_Higgs_L4[3][4][4][3] = LamSH;
  Curvature_Higgs_L4[4][0][0][4] = LamSH;
  Curvature_Higgs_L4[4][0][4][0] = LamSH;
  Curvature_Higgs_L4[4][1][1][4] = LamSH;
  Curvature_Higgs_L4[4][1][4][1] = LamSH;
  Curvature_Higgs_L4[4][2][2][4] = LamSH;
  Curvature_Higgs_L4[4][2][4][2] = LamSH;
  Curvature_Higgs_L4[4][3][3][4] = LamSH;
  Curvature_Higgs_L4[4][3][4][3] = LamSH;
  Curvature_Higgs_L4[4][4][0][0] = LamSH;
  Curvature_Higgs_L4[4][4][1][1] = LamSH;
  Curvature_Higgs_L4[4][4][2][2] = LamSH;
  Curvature_Higgs_L4[4][4][3][3] = LamSH;
  Curvature_Higgs_L4[4][4][4][4] = 0.12e2*LamS;

  Curvature_Gauge_G2H2[0][0][0][0] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[0][0][1][1] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[0][0][2][2] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[0][0][3][3] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[0][3][0][3] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][1][2] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][2][1] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][3][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][1][0][0] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[1][1][1][1] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[1][1][2][2] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[1][1][3][3] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[1][3][0][2] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][3][1][3] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][3][3][1] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][3][2][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][2][0][0] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[2][2][1][1] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[2][2][2][2] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[2][2][3][3] = (1.0/2.0)*std::pow(SMConstants.C_g, 2);
  Curvature_Gauge_G2H2[2][3][0][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][3][1][1] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][3][2][2] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][3][3][3] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][0][3] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][1][2] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][3][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][2][1] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][0][2] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][1][3] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][3][1] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][2][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][0][0] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][1][1] = (1.0/2.0)*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][2][2] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][3][3] = -1.0/2.0*SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][3][0][0] = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2[3][3][1][1] = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2[3][3][2][2] = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);
  Curvature_Gauge_G2H2[3][3][3][3] = (1.0/2.0)*std::pow(SMConstants.C_gs, 2);


   std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  std::complex<double> VC11, VC12, VC13, VC21, VC22, VC23, VC31, VC32, VC33;
  VC11 = std::conj(SMConstants.C_Vud);
  VC12 = std::conj(SMConstants.C_Vus);
  VC13 = std::conj(SMConstants.C_Vub);
  VC21 = std::conj(SMConstants.C_Vcd);
  VC22 = std::conj(SMConstants.C_Vcs);
  VC23 = std::conj(SMConstants.C_Vcb);
  VC31 = std::conj(SMConstants.C_Vtd);
  VC32 = std::conj(SMConstants.C_Vts);
  VC33 = std::conj(SMConstants.C_Vtb);

  Curvature_Lepton_F2H1[0][1][2] = II / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = 0.1e1 / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = II / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = 0.1e1 / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = II / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = 0.1e1 / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = II / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = 0.1e1 / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = II / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = 0.1e1 / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = II / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = 0.1e1 / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / vh * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / vh * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / vh * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / vh * SMConstants.C_MassTau;

  Curvature_Quark_F2H1[0][6][2] = -II / vh * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][6][3] = 0.1e1 / vh * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][9][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1] = II / vh * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / vh * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / vh * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2]  = -II / vh * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][7][3]  = 0.1e1 / vh * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][9][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1] = II / vh * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] =
      II / vh * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] =
      II / vh * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2] = -II / vh * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][8][3] = 0.1e1 / vh * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][9][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1] = II / vh * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / vh * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / vh * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0]  = 0.1e1 / vh * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][6][1]  = II / vh * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][7][0]  = V21 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][1]  = II * V21 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = 0.1e1 / vh * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]  = II / vh * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[3][9][2]  = II / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][3]  = 0.1e1 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][6][0]  = 0.1e1 / vh * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][6][1]  = II / vh * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][7][0]  = V22 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][8][0]  = 0.1e1 / vh * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][8][1]  = II / vh * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][10][2] = II / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = 0.1e1 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[5][6][0]  = V13 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][6][1]  = II / vh * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[5][7][0]  = V23 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][7][1]  = II / vh * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[5][8][0]  = V33 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][8][1]  = II / vh * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[5][11][2] = II / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = 0.1e1 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][0][2]  = -II / vh * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][0][3]  = 0.1e1 / vh * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][3][0]  = 0.1e1 / vh * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[6][3][1]  = II / vh * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[6][4][0]  = 0.1e1 / vh * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][4][1]  = II / vh * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][5][0]  = V13 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][5][1]  = II / vh * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[7][1][2]  = -II / vh * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][1][3]  = 0.1e1 / vh * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][3][0]  = V21 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][3][1]  = II * V21 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][4][0]  = V22 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][5][0]  = V23 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[7][5][1]  = II / vh * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[8][2][2]  = -II / vh * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][2][3]  = 0.1e1 / vh * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][3][0]  = 0.1e1 / vh * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]  = II / vh * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][4][0]  = 0.1e1 / vh * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][4][1]  = II / vh * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][5][0]  = V33 / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][5][1]  = II / vh * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[9][0][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1] = II / vh * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1] = II / vh * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1] = II / vh * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2] = II / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[9][3][3] = 0.1e1 / vh * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[10][0][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / vh * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] =
      II / vh * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / vh * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = II / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = 0.1e1 / vh * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] =
      -0.1e1 / vh * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / vh * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] =
      -0.1e1 / vh * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] =
      II / vh * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] =
      -0.1e1 / vh * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / vh * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = II / vh * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = 0.1e1 / vh * SMConstants.C_MassBottom;

  SetCurvatureDone = true;
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

  return 0;
}

double Class_RxSM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  return 0;
}

void Class_RxSM::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
