// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the electroweak phase transition for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for unique_ptr
#include <stdlib.h> // for EXIT_FAILURE, atoi
#include <string>   // for getline, string
#include <utility>  // for pair
#include <vector>   // for vector

// CB: added v
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/ClassPotentialR2HDM.h>
// CB: added ^

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int Line{};
  std::string InputFile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{true};

  CLIOptions(const BSMPT::parser &argparser);
  bool good() const;
};

BSMPT::parser prepare_parser();

std::vector<std::string> convert_input(int argc, char *argv[]);

int main(int argc, char *argv[])
try
{
  const auto SMConstants = GetSMConstants();
  auto argparser         = prepare_parser();
  argparser.add_input(convert_input(argc, argv));
  const CLIOptions args(argparser);
  if (not args.good())
  {
    return EXIT_FAILURE;
  }

  int linecounter = 1;
  std::ifstream infile(args.InputFile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default, "Input file not found ");
    return EXIT_FAILURE;
  }

  Logger::Write(LoggingLevel::ProgDetailed, "Found file");

  std::string linestr;
  // std::unique_ptr<Class_Potential_Origin> modelPointer =
  //     ModelID::FChoose(args.Model, SMConstants);
  std::shared_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model, SMConstants);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  while (getline(infile, linestr))
  {
    if (linecounter > args.Line) break;

    if (linecounter == 1)
    {

      modelPointer->setUseIndexCol(linestr);
    }
    if (linecounter == args.Line and linecounter != 1)
    {
      Logger::Write(LoggingLevel::ProgDetailed, "Found line");
      modelPointer->initModel(linestr);
      modelPointer->write();

      // CB: v
      double eps = 0.01;
      std::vector<double> vevTree = modelPointer->MinimizeOrderVEV(modelPointer->get_vevTreeMin());

      // auto VTree0 = [&](std::vector<double> x) {
      //   return modelPointer->VTree(x, 0, false);
      // };
      // auto VEff0 = [&](std::vector<double> x) {
      //   return modelPointer->VEff(x, 0, 0, 0);
      // };
      // auto VEff1 = [&](std::vector<double> x) {
      //   return modelPointer->VEff(x, 0, 0, 1);
      // };
      // std::vector<std::vector<double>> hessenum0 = HessianNumerical(vevTree, VTree0, eps);
      // std::vector<std::vector<double>> hessenum10 = HessianNumerical(vevTree, VEff0, eps);
      // std::vector<std::vector<double>> hessenum1 = HessianNumerical(vevTree, VEff1, eps);


      // std::function<double(std::vector<double>)> V0  = std::bind(&BSMPT::Class_Potential_Origin::VTree, modelPointer, std::placeholders::_1, 0, false);
      // std::function<double(std::vector<double>)> V10 = std::bind(&BSMPT::Class_Potential_Origin::VEff, modelPointer, std::placeholders::_1, 0, 0, 0);
      // std::function<double(std::vector<double>)> V1  = std::bind(&BSMPT::Class_Potential_Origin::VEff, modelPointer, std::placeholders::_1, 0, 0, 1);
      std::function<double(std::vector<double>)> V0  = std::bind(&BSMPT::Class_Potential_Origin::VTree, std::ref(modelPointer), std::placeholders::_1, 0, false);
      std::function<double(std::vector<double>)> V10 = std::bind(&BSMPT::Class_Potential_Origin::VEff, std::ref(modelPointer), std::placeholders::_1, 0, 0, 0);
      std::function<double(std::vector<double>)> V1  = std::bind(&BSMPT::Class_Potential_Origin::VEff, std::ref(modelPointer), std::placeholders::_1, 0, 0, 1);
      
      std::vector<std::vector<double>> hessenum0 = HessianNumerical(vevTree, V0, eps);
      std::vector<std::vector<double>> hessenum10 = HessianNumerical(vevTree, V10, eps);
      std::vector<std::vector<double>> hessenum1 = HessianNumerical(vevTree, V1, eps);


      auto MassMatrix = modelPointer->HiggsMassMatrix(vevTree, 0, 0);

      std::size_t NHiggs = modelPointer->get_NHiggs();

      std::cout << "HiggsMassMatrix:" << std::endl;
      for (std::size_t i = 0; i < NHiggs; ++i) {
        for (std::size_t j = 0; j < NHiggs; ++j) {
          std::cout << MassMatrix(i, j) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;

/*
      std::cout << "Hesse VTree:" << std::endl;
      for (std::size_t i = 0; i < NHiggs; ++i) {
        for (std::size_t j = 0; j < NHiggs; ++j) {
          std::cout << hessenum0[i][j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;

      std::cout << "Hesse VEff 0:" << std::endl;
      for (std::size_t i = 0; i < NHiggs; ++i) {
        for (std::size_t j = 0; j < NHiggs; ++j) {
          std::cout << hessenum10[i][j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;

      std::cout << "Hesse VEff 1:" << std::endl;
      for (std::size_t i = 0; i < NHiggs; ++i) {
        for (std::size_t j = 0; j < NHiggs; ++j) {
          std::cout << hessenum1[i][j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
*/

      // for (std::size_t i = 0; i < NHiggs; ++i) {
      //   for (std::size_t j = 0; j < NHiggs; ++j) {
      //     MassMatrix(i, j) = hessenum1[i][j];
      //   }
      // }

      SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
      auto EV = es.eigenvalues();
      for (std::size_t i = 0; i < NHiggs; ++i) {
        if (EV[i] < 0) {
          EV[i] = EV[i]/std::abs(EV[i])*std::sqrt(std::abs(EV[i]));
        } else {
          EV[i] = std::sqrt(EV[i]);
        }
      }

      std::cout << "EV=" << EV << std::endl;

      // Test at tree-level
      modelPointer->SetUseTreeLevel(true);
      std::shared_ptr<MinimumTracer> mintracer(new MinimumTracer(modelPointer, 1 + 2 + 4, true));
      auto glob_min_0 = mintracer->ConvertToVEVDim(mintracer->GetGlobalMinimum(0));

      // Test at 1-loop effective
      modelPointer->SetUseTreeLevel(false);
      mintracer.reset(new MinimumTracer(modelPointer, 1 + 2 + 4, true));
      auto glob_min = mintracer->ConvertToVEVDim(mintracer->GetGlobalMinimum(0));

      std::cout << "glob_min_0(tree-level)=" << glob_min_0 << std::endl;
      std::cout << std::sqrt(glob_min_0[0]*glob_min_0[0] + glob_min_0[1]*glob_min_0[1] + glob_min_0[2]*glob_min_0[2] + glob_min_0[3]*glob_min_0[3]) << std::endl;

      std::cout << "glob_min(1-loop)=" << glob_min << std::endl;
      std::cout << std::sqrt(glob_min[0]*glob_min[0] + glob_min[1]*glob_min[1] + glob_min[2]*glob_min[2] + glob_min[3]*glob_min[3]) << std::endl;

      std::vector<double> vevGlob = modelPointer->MinimizeOrderVEV(glob_min);

      std::vector<std::vector<double>> hessenum1Glob = HessianNumerical(vevGlob, V1, eps);
      for (std::size_t i = 0; i < NHiggs; ++i) {
        for (std::size_t j = 0; j < NHiggs; ++j) {
          MassMatrix(i, j) = hessenum1Glob[i][j];
        }
      }

      // std::cout << "hessenum1Glob:" << std::endl;
      // for (std::size_t i = 0; i < NHiggs; ++i) {
      //   for (std::size_t j = 0; j < NHiggs; ++j) {
      //     std::cout << MassMatrix(i, j) << " ";
      //   }
      //   std::cout << std::endl;
      // }
      // std::cout << std::endl;

      std::cout << std::endl;

      SelfAdjointEigenSolver<MatrixXd> esGlob(MassMatrix, EigenvaluesOnly);
      EV = esGlob.eigenvalues();
      for (std::size_t i = 0; i < NHiggs; ++i) {
        if (EV[i] < 0) {
          EV[i] = EV[i]/std::abs(EV[i])*std::sqrt(std::abs(EV[i]));
        } else {
          EV[i] = std::sqrt(EV[i]);
        }
      }

      std::cout << "EV=" << EV << std::endl;

      int ewsr_status = mintracer->IsThereEWSymmetryRestoration();
      auto status_ewsr = mintracer->GetStatusEWSR(ewsr_status);
      std::cout << "EWSR status=" << ewsr_status << ", " << status_ewsr << std::endl;
      // ewsr_status:
      //   3: EWSR
      //   2: no EWSR (i.e. potential is BFB, but minimum is not at the origin, or what does this mean?)
      //   1, 0: higher orders needed/not conclusive
      //  -1: not bounded from below at high T -> no EWSR

      std::cout << std::endl;
      // CB: ^

      std::vector<double> dummy;
      modelPointer->Debugging(dummy, dummy);
      ModelTests::CheckImplementation(*modelPointer, args.WhichMinimizer);
    }
    linecounter++;
    if (infile.eof()) break;
  }
  return EXIT_SUCCESS;
}
catch (int)
{
  return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  return EXIT_FAILURE;
}

bool CLIOptions::good() const
{
  if (UseGSL and not Minimizer::UseGSLDefault)
  {
    throw std::runtime_error(
        "You set --useGSL=true but GSL was not found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error(
        "You set --useCMAES=true but CMAES was not found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error(
        "You set --useNLopt=true but NLopt was not found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }
  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  if (Line < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
    return false;
  }
  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  argparser.check_required_parameters();
  Model     = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile = argparser.get_value("input");
  Line      = argparser.get_value<int>("line");

  try
  {
    UseGSL = argparser.get_value<bool>("useGSL");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseCMAES = argparser.get_value<bool>("useCMAES");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseNLopt = argparser.get_value<bool>("useNLopt");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseMultithreading = argparser.get_value<bool>("useMultithreading");
  }
  catch (BSMPT::parserException &)
  {
  }

  WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser;
  argparser.add_argument("model", "The model you want to investigate.", true);
  argparser.add_argument("input", "The input file in tsv format.", true);
  argparser.add_argument(
      "line",
      "The line in the input file with the parameter point used to "
      "check the model.",
      true);

  std::stringstream ss;
  ss << "Test performs a serious of tests on the given model. "
        "Intended for testing new models."
     << std::endl
     << "It is called either by " << std::endl
     << "./Test model input Line" << std::endl
     << "or with the following arguments" << std::endl;
  argparser.set_help_header(ss.str());

  argparser.enable_minimizer_options();

  return argparser;
}

std::vector<std::string> convert_input(int argc, char *argv[])
{
  std::vector<std::string> arguments;
  if (argc == 1) return arguments;
  auto first_arg = std::string(argv[1]);

  bool UsePrefix =
      StringStartsWith(first_arg, "--") or StringStartsWith(first_arg, "-");

  if (UsePrefix)
  {
    for (int i{1}; i < argc; ++i)
    {
      arguments.emplace_back(argv[i]);
    }
  }
  else
  {
    if (argc >= 2)
    {
      arguments.emplace_back("--model=" + std::string(argv[1]));
    }
    if (argc >= 3)
    {
      arguments.emplace_back("--input=" + std::string(argv[2]));
    }
    if (argc >= 4)
    {
      arguments.emplace_back("--line=" + std::string(argv[3]));
    }
  }
  return arguments;
}
