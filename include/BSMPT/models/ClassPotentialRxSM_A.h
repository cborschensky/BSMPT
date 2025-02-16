// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Model file for the SM + real singlet with the tree-level potential
 * Alain Verduras Schaeidt
 */

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_Template class
 * Template for implementing a new model
 */
class Class_RxSM : public Class_Potential_Origin
{
public:
  Class_RxSM(const ISMConstants &smConstants);
  virtual ~Class_RxSM();

  // Add here your parameters for the Lagrangian as well as for the counterterm
  // potential Add here your variables in which you will save the Debye
  // correction factors

  double mu2, MS, Lam, LamS, LamSH, KapS, KapSH;

  double dmu2, dMS, dLam, dLamS, dLamSH, dKapS, dKapSH, dT1, dT2, dT3, dT4, dT5;
  
  double vh, vs;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;
  
  /**
   * @brief set_gen
   * @param par[0] = v
   * @param par[1] = vs
   * @param par[2] = mu2
   * @param par[3] = MS
   * @param par[4] = Kaps
   * @param par[5] = KapSH
   * @param par[6] = Lam
   * @param par[7] = LamS
   * @param par[8] = LamSH
   */

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void AdjustRotationMatrix() override;
  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
