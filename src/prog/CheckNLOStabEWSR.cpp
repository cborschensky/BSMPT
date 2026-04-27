// Copyright (C) 2026 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program checks NLO stability and electroweak symmetry restoration
 *
 */

#include "BSMPT/minimum_tracer/minimum_tracer.h"   // MinimumTracer
#include <BSMPT/models/IncludeAllModels.h>
#include "BSMPT/transition_tracer/transition_tracer.h"   // TransitionTracer
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <utility>  // for pair
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int firstline{0}, lastline{0};
  std::string inputfile, outputfile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{false};
  int CheckEWSymmetryRestoration{1};
  int CheckNLOStability{1};

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

  std::ifstream infile(args.inputfile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Input file " + args.inputfile + " not found ");
    return EXIT_FAILURE;
  }

  Logger::Write(LoggingLevel::ProgDetailed, "Found file");

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model, SMConstants);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  std::string linestr, linestr_store;
  int linecounter   = 1;

  // output contents storage
  std::stringstream output_contents;
  std::vector<std::string> legend;

  legend.push_back("status_nlo_stability");
  legend.push_back("status_ewsr");
  legend.push_back("runtime");

  // write to output file
  std::ofstream outfile(args.outputfile);
  if (!outfile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Can not create file " + args.outputfile);
    return EXIT_FAILURE;
  }

  while (getline(infile, linestr))
  {
    if (linecounter == 1) linestr_store = linestr;

    if (linecounter > args.lastline) break;

    if (linecounter >= args.firstline and linecounter <= args.lastline)
    {
      output_contents.precision(
          std::numeric_limits<double>::max_digits10);

      Logger::Write(LoggingLevel::ProgDetailed,
                    "Currently at line " + std::to_string(linecounter));

      modelPointer->setUseIndexCol(linestr_store);

      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);

      if (args.firstline == args.lastline)
      {
        modelPointer->write();
      }

      auto start = std::chrono::high_resolution_clock::now();

      std::shared_ptr<MinimumTracer> mintracer(new MinimumTracer(
          modelPointer, args.WhichMinimizer, args.UseMultithreading));

      status_codes status;

      // NLO stability check
      if (args.CheckNLOStability)
      {
        auto glob_min = mintracer->ConvertToVEVDim(mintracer->GetGlobalMinimum(0));
        Logger::Write(LoggingLevel::TransitionDetailed,
                      "Global minimum at T = 0 found at " +
                          vec_to_string(glob_min));
        status.status_nlo_stability =
            mintracer->GetStatusNLOVEV(modelPointer->CheckNLOVEV(glob_min));
      }
      else
      {
        Logger::Write(LoggingLevel::TransitionDetailed,
                      "Check for NLO stability is disabled.");
        status.status_nlo_stability = StatusNLOStability::Off;
      }

      if (status.status_nlo_stability ==
              StatusNLOStability::Success or
          status.status_nlo_stability ==
              StatusNLOStability::Off)
      {
        // Electroweak Symmetry Restoration check
        if (args.CheckEWSymmetryRestoration > 0)
        {
          double ewsr_status = mintracer->IsThereEWSymmetryRestoration();
          status.status_ewsr = mintracer->GetStatusEWSR(ewsr_status);

          // If no minimum was found at high temperature
          if (args.CheckEWSymmetryRestoration == 2 && ewsr_status < 2)
          {

            Logger::Write(
                LoggingLevel::TransitionDetailed,
                "EW symmetry restoration check failed. Point will be filtered out");
          }
          // If EW was not restored
          if (args.CheckEWSymmetryRestoration == 3 && ewsr_status < 3)
          {
            Logger::Write(
                LoggingLevel::TransitionDetailed,
                "EW symmetry restoration check failed. Point will be filtered out");
          }
        }
        else
        {
          Logger::Write(LoggingLevel::TransitionDetailed,
                        "Check for EW symmetry restoration is disabled.");
          status.status_ewsr = BSMPT::StatusEWSR::Off;
        }
      }


      auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::high_resolution_clock::now() - start)
                      .count() /
                  1000.;

      BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
                           "\nTook\t" + std::to_string(time) + " seconds.\n");

      output_contents
          << linestr << sep << parameters.second << sep
          << status.status_nlo_stability << sep
          << status.status_ewsr << sep << time;

      std::stringstream full_legend;
      full_legend << linestr_store << sep << modelPointer->addLegendCT() << sep
                  << legend;
      outfile << full_legend.str() << std::endl;
      outfile << output_contents.str() << std::endl;
    }

    linecounter++;
    if (infile.eof()) break;
  }
  outfile.close();
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
    throw std::runtime_error("You set --UseGSL=true but GSL was not "
                             "found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error("You set --UseCMAES=true but CMAES was not "
                             "found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error("You set --UseNLopt=true but NLopt was not "
                             "found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }

  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(LoggingLevel::Default,
                  "Your Model parameter does not match with the "
                  "implemented Models.");
    ShowInputError();
    return false;
  }
  if (firstline == 0 or lastline == 0)
  {
    Logger::Write(LoggingLevel::Default, "firstline or lastline not set.");
    return false;
  }
  if (firstline < 0 or lastline < 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid input for first- or lastline.");
    return false;
  }
  if (firstline > lastline)
  {
    Logger::Write(LoggingLevel::Default, "lastline is smaller then firstline.");
    return false;
  }
  if (CheckEWSymmetryRestoration > 2 or CheckEWSymmetryRestoration < 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid choice for CheckEWSymmetryRestoration.");
    return false;
  }
  if (CheckNLOStability < 0 or CheckNLOStability > 1)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid choice for CheckNLOStability.");
    return false;
  }

  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  std::stringstream ss;
  argparser.check_required_parameters();

  // required arguments
  Model      = BSMPT::ModelID::getModel(argparser.get_value("model"));
  inputfile  = argparser.get_value("input");
  outputfile = argparser.get_value("output");
  firstline  = argparser.get_value<int>("firstline");
  lastline   = argparser.get_value<int>("lastline");

  std::string GSLhelp   = Minimizer::UseGSLDefault ? "true" : "false";
  std::string CMAEShelp = Minimizer::UseLibCMAESDefault ? "true" : "false";
  std::string NLoptHelp = Minimizer::UseNLoptDefault ? "true" : "false";
  try
  {
    UseGSL = (argparser.get_value("usegsl") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usegsl not set, using default value: " << GSLhelp << "\n";
  }

  try
  {
    UseCMAES = (argparser.get_value("usecmaes") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usecmaes not set, using default value: " << CMAEShelp << "\n";
  }

  try
  {
    UseNLopt = (argparser.get_value("usenlopt") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usenlopt not set, using default value: " << NLoptHelp << "\n";
  }

  try
  {
    UseMultithreading = (argparser.get_value("usemultithreading") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usemultithreading not set, using default value: false\n";
  }

  // CheckNLOStability
  try
  {
    auto nlo_string = argparser.get_value("checknlo");
    if (nlo_string == "on")
    {
      CheckNLOStability = 1;
    }
    else if (nlo_string == "off")
    {
      CheckNLOStability = 0;
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--checknlo not set, using default value: on\n";
  }

  // CheckEWSymmetryRestoration
  try
  {
    auto ewsr_string = argparser.get_value("checkewsr");
    if (ewsr_string == "on")
    {
      CheckEWSymmetryRestoration = 1;
    }
    else if (ewsr_string == "off")
    {
      CheckEWSymmetryRestoration = 0;
    }
    else if (ewsr_string == "keep_bfb")
    {
      CheckEWSymmetryRestoration = 2;
    }
    else if (ewsr_string == "keep_ewsr")
    {
      CheckEWSymmetryRestoration = 3;
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--checkewsr not set, using default value: on\n";
  }

  WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);

  Logger::Write(LoggingLevel::ProgDetailed, ss.str());
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser(true);
  argparser.add_argument("help", "shows this menu", false);
  argparser.add_argument("model", "[*] model name", true);
  argparser.add_argument("input", "[*] input file (in tsv format)", true);
  argparser.add_argument("output", "[*] output file (in tsv format)", true);
  argparser.add_argument(
      "firstline", "[*] line number of first line in input file", true);
  argparser.add_subtext("    (expects line 1 to be a legend)");
  argparser.add_argument(
      "lastline", "[*] line number of last line in input file", true);
  argparser.add_argument("checknlo", "check for NLO stability", "on", false);
  argparser.add_subtext("on: only keep NLO stable points");
  argparser.add_subtext("off: check disabled");
  argparser.add_argument(
      "checkewsr", "check for EWSR at high temperature", "on", false);
  argparser.add_subtext("on: perform check and add info");
  argparser.add_subtext("keep_bfb: only keep BFB points");
  argparser.add_subtext("keep_ewsr: only keep EWSR points");
  argparser.add_subtext("off: check disabled");

  std::string GSLhelp   = Minimizer::UseGSLDefault ? "true" : "false";
  std::string CMAEShelp = Minimizer::UseLibCMAESDefault ? "true" : "false";
  std::string NLoptHelp = Minimizer::UseNLoptDefault ? "true" : "false";

  argparser.add_argument(
      "usegsl", "use GSL library for minimization", GSLhelp, false);
  argparser.add_argument(
      "usecmaes", "use CMAES library  for minimization", CMAEShelp, false);
  argparser.add_argument(
      "usenlopt", "use NLopt library for minimization", NLoptHelp, false);
  argparser.add_argument("usemultithreading",
                         "enable multi-threading for minimizers",
                         "false",
                         false);
  argparser.add_argument(
      "json", "use a json file instead of cli parameters", false);

  std::stringstream ss;
  ss << "CalcTemps calculates characteristic temperatures for phase "
        "transitions\nit is called "
        "by\n\n\t./bin/CalcTemps model input output firstline "
        "lastline\n\nor "
        "with arguments\n\n\t./bin/CalcTemps [arguments]\n\nwith the "
        "following arguments, ([*] are required arguments, others "
        "are optional):\n";
  argparser.set_help_header(ss.str());

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
      arguments.emplace_back("--output=" + std::string(argv[3]));
    }
    if (argc >= 5)
    {
      arguments.emplace_back("--firstline=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--lastline=" + std::string(argv[5]));
    }
  }
  return arguments;
}
