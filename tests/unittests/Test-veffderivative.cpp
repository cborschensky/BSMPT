#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using Approx = Catch::Approx;

#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>

TEST_CASE("Check of analytical T-derivative of Veff (CxSM example point)",
          "[veffderivative]")
{
  const std::vector<double> example_point_CXSM{/* v = */ 245.34120667410863,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* msq = */ -15650,
                                               /* lambda = */ 0.52,
                                               /* delta2 = */ 0.55,
                                               /* b2 = */ -8859,
                                               /* d2 = */ 0.5,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SMConstants);
  modelPointer->initModel(example_point_CXSM);

  std::vector<double> T_values{20., 50., 100., 200., 500., 1000.};
  std::vector<std::vector<double>> phi_values{{0., 0., 0.},
                                              {100., 0., 0.},
                                              {0., 100., 0.},
                                              {0., 0., 100.},
                                              {100., 200., 0.},
                                              {50., 100., 200.}};

  auto dTV_num = [=](std::vector<double> p, double T, double h)
  {
    double V1 = modelPointer->VEff(modelPointer->MinimizeOrderVEV(p), T - h, 0);
    double V2 = modelPointer->VEff(modelPointer->MinimizeOrderVEV(p), T + h, 0);
    return (V2 - V1) / 2. / h;
  };

  for (double T_val : T_values)
  {
    for (auto &phi_val : phi_values)
    {
      double x_num = dTV_num(phi_val, T_val, 1e-3);
      double x_an  = modelPointer->VEff(
          modelPointer->MinimizeOrderVEV(phi_val), T_val, -1);
      REQUIRE(x_num == Approx(x_an).epsilon(5e-3));
    }
  }
}

TEST_CASE("Check of analytical T-derivative of Veff (R2HDM example point)",
          "[veffderivative]")
{
  const std::vector<double> example_point_R2HDM{
      /* lambda_1 = */ 6.9309437685026,
      /* lambda_2 = */ 0.26305141403285998,
      /* lambda_3 = */ 1.2865950045595,
      /* lambda_4 = */ 4.7721306931875001,
      /* lambda_5 = */ 4.7275722046239004,
      /* m_{12}^2 = */ 18933.440789693999,
      /* tan(beta) = */ 16.577896825227999,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  std::vector<double> T_values{20., 50., 100., 200., 500., 1000.};
  std::vector<std::vector<double>> phi_values{{0., 0., 0., 0.},
                                              {100., 0., 0., 0.},
                                              {0., 100., 0., 0.},
                                              {0., 0., 100., 0.},
                                              {0., 50., 200., 0.},
                                              {50., 100., 200., 10.}};

  auto dTV_num = [=](std::vector<double> p, double T, double h)
  {
    double V1 = modelPointer->VEff(modelPointer->MinimizeOrderVEV(p), T - h, 0);
    double V2 = modelPointer->VEff(modelPointer->MinimizeOrderVEV(p), T + h, 0);
    return (V2 - V1) / 2. / h;
  };

  for (double T_val : T_values)
  {
    for (auto &phi_val : phi_values)
    {
      double x_num = dTV_num(phi_val, T_val, 1e-3);
      double x_an  = modelPointer->VEff(
          modelPointer->MinimizeOrderVEV(phi_val), T_val, -1);
      REQUIRE(x_num == Approx(x_an).epsilon(5e-3));
    }
  }
}
