// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file transition history evaluator
 */

#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{

TransitionTracer::TransitionTracer(user_input &input)
{
  num_vev = input.modelPointer->get_nVEV();

  std::shared_ptr<MinimumTracer> mintracer(new MinimumTracer(
      input.modelPointer, input.which_minimizer, input.use_multithreading));

  // initialize legend
  output_store.legend = mintracer->GetLegend(0, input.gw_calculation);

  // NLO stability check
  if (input.nlo_check)
  {
    auto glob_min = mintracer->ConvertToVEVDim(mintracer->GetGlobalMinimum(0));
    Logger::Write(LoggingLevel::TransitionDetailed,
                  "Global minimum at T = 0 found at " +
                      vec_to_string(glob_min));
    output_store.status.status_nlo_stability =
        mintracer->GetStatusNLOVEV(input.modelPointer->CheckNLOVEV(glob_min));
  }
  else
  {
    Logger::Write(LoggingLevel::TransitionDetailed,
                  "Check for NLO stability is disabled.");
    output_store.status.status_nlo_stability = BSMPT::StatusNLOStability::Off;
  }

  if (output_store.status.status_nlo_stability ==
          BSMPT::StatusNLOStability::Success or
      output_store.status.status_nlo_stability ==
          BSMPT::StatusNLOStability::Off)
  {
    // Electroweak Symmetry Restoration check
    bool point_filtered_out_due_to_snr = false;
    if (input.ewsr_check > 0)
    {
      double ewsr_status = mintracer->IsThereEWSymmetryRestoration();
      output_store.status.status_ewsr = mintracer->GetStatusEWSR(ewsr_status);

      // If no minimum was found at high temperature
      if (input.ewsr_check == 2 && ewsr_status < 2)
      {

        Logger::Write(
            LoggingLevel::TransitionDetailed,
            "EW symmetry restoration check failed. Point will be filtered out");

        point_filtered_out_due_to_snr = true;
      }
      // If EW was not restored
      if (input.ewsr_check == 3 && ewsr_status < 3)
      {
        Logger::Write(
            LoggingLevel::TransitionDetailed,
            "EW symmetry restoration check failed. Point will be filtered out");

        point_filtered_out_due_to_snr = true;
      }
    }
    else
    {
      Logger::Write(LoggingLevel::TransitionDetailed,
                    "Check for EW symmetry restoration is disabled.");
      output_store.status.status_ewsr = StatusEWSR::Off;
    }

    if (not point_filtered_out_due_to_snr)
    {
      Logger::Write(
          LoggingLevel::TransitionDetailed,
          "Track phases in between T_low = " + std::to_string(input.T_low) +
              " GeV and T_high = " + std::to_string(input.T_high) + " GeV");

      Vacuum vac(input.T_low,
                 input.T_high,
                 mintracer,
                 input.modelPointer,
                 input.multistepmode,
                 input.num_points);

      vec_coex = vac.CoexPhasesList;

      output_store.num_coex_phase_pairs = vec_coex.size();

      Logger::Write(LoggingLevel::TransitionDetailed,
                    "\nIdentified " +
                        std::to_string(output_store.num_coex_phase_pairs) +
                        " coexisiting phase pair(s) in total.");

      output_store.status.status_tracing    = vac.status_vacuum;
      output_store.status.status_coex_pairs = vac.status_coex_pairs;

      if ((output_store.status.status_tracing == StatusTracing::Success) &&
          (output_store.status.status_coex_pairs == StatusCoexPair::Success))
      {
        output_store.legend = mintracer->GetLegend(
            output_store.num_coex_phase_pairs, input.gw_calculation);

        for (auto pair : vec_coex)
        {
          transition_data new_transition_data;
          gw_data new_gw_data;

          Logger::Write(LoggingLevel::TransitionDetailed,
                        "Pair " + std::to_string(pair.coex_pair_id) +
                            " (phase " + std::to_string(pair.false_phase.id) +
                            " -> phase " + std::to_string(pair.true_phase.id) +
                            ") with Tc = " + std::to_string(pair.crit_temp) +
                            " (" + StatusCritToString.at(pair.crit_status) +
                            ")");

          output_store.status.status_crit.push_back(pair.crit_status);
          if ((pair.crit_status == BSMPT::StatusCrit::Success) ||
              (pair.crit_status == BSMPT::StatusCrit::TrueLower))
          {

            Logger::Write(LoggingLevel::TransitionDetailed,
                          "Calculate bounce solution, for more output, use "
                          "--logginglevel::bouncedetailed=true.");

            new_transition_data.crit_temp = pair.crit_temp;
            // CB changes v

            // new_transition_data.crit_true_vev =
                // pair.true_phase.Get(pair.crit_temp).point;
            // new_transition_data.crit_false_vev =
                // pair.false_phase.Get(pair.crit_temp).point;

            Minimum true_phase_crit = pair.true_phase.Get(pair.crit_temp);
            Minimum false_phase_crit = pair.false_phase.Get(pair.crit_temp);
            new_transition_data.crit_true_vev = true_phase_crit.point;
            new_transition_data.crit_false_vev = false_phase_crit.point;
            // Calculate the difference of Veff between the false and the true minimum;
            // in the case of the critical temperature, this should be zero or close to zero
            new_transition_data.crit_deltaVif = false_phase_crit.potential - true_phase_crit.potential;
            // Calculate (m(T)/T)^2
            auto [crit_m2T2_scalar, crit_m2T2_gauge] = CheckMassRatio_sep(
                input, new_transition_data.crit_false_vev, pair.crit_temp);
            new_transition_data.crit_m2T2_scalar = crit_m2T2_scalar;
            new_transition_data.crit_m2T2_gauge = crit_m2T2_gauge;
            // CB changes ^

            // (void)CheckMassRatio(
                // input, new_transition_data.crit_false_vev, pair.crit_temp);


            BounceSolution bounce(input.modelPointer,
                                  mintracer,
                                  pair,
                                  input.vwall,
                                  input.epsturb,
                                  input.maxpathintegrations,
                                  input.number_of_initial_scan_temperatures,
                                  input.PNLO_scaling);

            ListBounceSolution.push_back(bounce);

            output_store.status.status_bounce_sol.push_back(
                bounce.status_bounce_sol);

            if (bounce.status_bounce_sol == StatusGW::Success)
            {
              bounce.CalculateNucleationTempApprox();

              output_store.status.status_nucl_approx.push_back(
                  bounce.status_nucl_approx);
              if (bounce.status_nucl_approx ==
                  BSMPT::StatusTemperature::Success)
              {
                new_transition_data.nucl_approx_temp =
                    bounce.GetNucleationTempApprox();
                // CB changes v

                // new_transition_data.nucl_approx_true_vev =
                //     pair.true_phase
                //         .Get(new_transition_data.nucl_approx_temp.value_or(
                //             EmptyValue))
                //         .point;
                // new_transition_data.nucl_approx_false_vev =
                //     pair.false_phase
                //         .Get(new_transition_data.nucl_approx_temp.value_or(
                //             EmptyValue))
                //         .point;

                Minimum true_phase_nucl_approx = pair.true_phase.Get(new_transition_data.nucl_approx_temp.value_or(EmptyValue));
                Minimum false_phase_nucl_approx = pair.false_phase.Get(new_transition_data.nucl_approx_temp.value_or(EmptyValue));
                new_transition_data.nucl_approx_true_vev = true_phase_nucl_approx.point;
                new_transition_data.nucl_approx_false_vev = false_phase_nucl_approx.point;
                // Calculate the difference of Veff between the false and the true minimum
                new_transition_data.nucl_approx_deltaVif = false_phase_nucl_approx.potential - true_phase_nucl_approx.potential;
                // Calculate (m(T)/T)^2
                auto [nucl_approx_m2T2_scalar, nucl_approx_m2T2_gauge] = CheckMassRatio_sep(input,
                                     new_transition_data.nucl_approx_false_vev,
                                     bounce.GetNucleationTempApprox());
                new_transition_data.nucl_approx_m2T2_scalar = nucl_approx_m2T2_scalar;
                new_transition_data.nucl_approx_m2T2_gauge = nucl_approx_m2T2_gauge;
                // CB changes ^

                // (void)CheckMassRatio(input,
                //                      new_transition_data.nucl_approx_false_vev,
                //                      bounce.GetNucleationTempApprox());
              }
              else
              {
                new_transition_data.nucl_approx_true_vev =
                    std::vector<double>(num_vev, EmptyValue);
                new_transition_data.nucl_approx_false_vev =
                    std::vector<double>(num_vev, EmptyValue);
              }

              bounce.CalculateNucleationTemp();

              output_store.status.status_nucl.push_back(bounce.status_nucl);
              if (bounce.status_nucl == BSMPT::StatusTemperature::Success)
              {
                new_transition_data.nucl_temp = bounce.GetNucleationTemp();
                // CB changes v

                // new_transition_data.nucl_true_vev =
                //     pair.true_phase
                //         .Get(new_transition_data.nucl_temp.value_or(EmptyValue))
                //         .point;
                // new_transition_data.nucl_false_vev =
                //     pair.false_phase
                //         .Get(new_transition_data.nucl_temp.value_or(EmptyValue))
                //         .point;

                Minimum true_phase_nucl = pair.true_phase.Get(new_transition_data.nucl_temp.value_or(EmptyValue));
                Minimum false_phase_nucl = pair.false_phase.Get(new_transition_data.nucl_temp.value_or(EmptyValue));
                new_transition_data.nucl_true_vev = true_phase_nucl.point;
                new_transition_data.nucl_false_vev = false_phase_nucl.point;
                // Calculate the difference of Veff between the false and the true minimum
                new_transition_data.nucl_deltaVif = false_phase_nucl.potential - true_phase_nucl.potential;
                // Calculate (m(T)/T)^2
                auto [nucl_m2T2_scalar, nucl_m2T2_gauge] = CheckMassRatio_sep(input,
                                     new_transition_data.nucl_false_vev,
                                     bounce.GetNucleationTemp());
                new_transition_data.nucl_m2T2_scalar = nucl_m2T2_scalar;
                new_transition_data.nucl_m2T2_gauge = nucl_m2T2_gauge;
                // CB changes ^

               // (void)CheckMassRatio(input,
               //                       new_transition_data.nucl_false_vev,
               //                       bounce.GetNucleationTemp());
              }
              else
              {
                new_transition_data.nucl_true_vev =
                    std::vector<double>(num_vev, EmptyValue);
                new_transition_data.nucl_false_vev =
                    std::vector<double>(num_vev, EmptyValue);
              }

              bounce.CalculatePercolationTemp();

              output_store.status.status_perc.push_back(bounce.status_perc);
              if (bounce.status_perc == BSMPT::StatusTemperature::Success)
              {
                new_transition_data.perc_temp = bounce.GetPercolationTemp();
                // CB changes v

                // new_transition_data.perc_true_vev =
                //     pair.true_phase
                //         .Get(new_transition_data.perc_temp.value_or(EmptyValue))
                //         .point;
                // new_transition_data.perc_false_vev =
                //     pair.false_phase
                //         .Get(new_transition_data.perc_temp.value_or(EmptyValue))
                //         .point;

                Minimum true_phase_perc = pair.true_phase.Get(new_transition_data.perc_temp.value_or(EmptyValue));
                Minimum false_phase_perc = pair.false_phase.Get(new_transition_data.perc_temp.value_or(EmptyValue));
                new_transition_data.perc_true_vev = true_phase_perc.point;
                new_transition_data.perc_false_vev = false_phase_perc.point;
                // Calculate the difference of Veff between the false and the true minimum
                new_transition_data.perc_deltaVif = false_phase_perc.potential - true_phase_perc.potential;
                // Calculate (m(T)/T)^2
                auto [perc_m2T2_scalar, perc_m2T2_gauge] = CheckMassRatio_sep(input,
                                     new_transition_data.perc_false_vev,
                                     bounce.GetPercolationTemp());
                new_transition_data.perc_m2T2_scalar = perc_m2T2_scalar;
                new_transition_data.perc_m2T2_gauge = perc_m2T2_gauge;
                // CB changes ^

                // (void)CheckMassRatio(input,
                //                      new_transition_data.perc_false_vev,
                //                      bounce.GetPercolationTemp());
              }
              else
              {
                new_transition_data.perc_true_vev =
                    std::vector<double>(num_vev, EmptyValue);
                new_transition_data.perc_false_vev =
                    std::vector<double>(num_vev, EmptyValue);
              }

              bounce.CalculateCompletionTemp();

              output_store.status.status_compl.push_back(bounce.status_compl);
              if (bounce.status_compl == BSMPT::StatusTemperature::Success)
              {
                new_transition_data.compl_temp = bounce.GetCompletionTemp();
                // CB changes v

                // new_transition_data.compl_true_vev =
                //     pair.true_phase
                //         .Get(
                //             new_transition_data.compl_temp.value_or(EmptyValue))
                //         .point;
                // new_transition_data.compl_false_vev =
                //     pair.false_phase
                //         .Get(
                //             new_transition_data.compl_temp.value_or(EmptyValue))
                //         .point;

                Minimum true_phase_compl = pair.true_phase.Get(new_transition_data.compl_temp.value_or(EmptyValue));
                Minimum false_phase_compl = pair.false_phase.Get(new_transition_data.compl_temp.value_or(EmptyValue));
                new_transition_data.compl_true_vev = true_phase_compl.point;
                new_transition_data.compl_false_vev = false_phase_compl.point;
                // Calculate the difference of Veff between the false and the true minimum
                new_transition_data.compl_deltaVif = false_phase_compl.potential - true_phase_compl.potential;
                // Calculate (m(T)/T)^2
                auto [compl_m2T2_scalar, compl_m2T2_gauge] = CheckMassRatio_sep(input,
                                     new_transition_data.compl_false_vev,
                                     bounce.GetCompletionTemp());
                new_transition_data.compl_m2T2_scalar = compl_m2T2_scalar;
                new_transition_data.compl_m2T2_gauge = compl_m2T2_gauge;
                // CB changes ^

                // (void)CheckMassRatio(input,
                //                      new_transition_data.compl_false_vev,
                //                      bounce.GetCompletionTemp());
              }
              else
              {
                new_transition_data.compl_true_vev =
                    std::vector<double>(num_vev, EmptyValue);
                new_transition_data.compl_false_vev =
                    std::vector<double>(num_vev, EmptyValue);
              }

              BSMPT::StatusTemperature trans_status =
                  BSMPT::StatusTemperature::NotSet;
              if (input.which_transition_temp ==
                  TransitionTemperature::ApproxNucleation)
              {
                trans_status = bounce.status_nucl_approx;
              }
              else if (input.which_transition_temp ==
                       TransitionTemperature::Nucleation)
              {
                trans_status = bounce.status_nucl;
              }
              else if (input.which_transition_temp ==
                       TransitionTemperature::Percolation)
              {
                trans_status = bounce.status_perc;
              }
              else if (input.which_transition_temp ==
                       TransitionTemperature::Completion)
              {
                trans_status = bounce.status_compl;
              }

              if (trans_status == BSMPT::StatusTemperature::Success &&
                  input.gw_calculation)
              {
                Logger::Write(LoggingLevel::TransitionDetailed,
                              "Start GW parameters calculation.");

                GravitationalWave gw(bounce, input.which_transition_temp);

                new_gw_data.status_gw  = gw.data.status;
                new_gw_data.trans_temp = gw.data.transitionTemp;
                new_gw_data.reh_temp   = gw.data.reheatingTemp;

                new_gw_data.alpha       = gw.data.PTStrength;
                new_gw_data.beta_over_H = gw.data.betaH;
                new_gw_data.vwall       = gw.data.vw;

                if (new_gw_data.status_gw != StatusGW::Failure)
                {
                  gw.CalcPeakCollision();
                  new_gw_data.fb_col = gw.data.CollisionParameter.f_b.value();
                  new_gw_data.omegab_col =
                      gw.data.CollisionParameter.Omega_b.value();

                  gw.CalcPeakSoundWave();
                  new_gw_data.f1_sw = gw.data.SoundWaveParameter.f_1.value();
                  new_gw_data.f2_sw = gw.data.SoundWaveParameter.f_2.value();
                  new_gw_data.omega_2_sw =
                      gw.data.SoundWaveParameter.Omega_2.value();

                  gw.CalcPeakTurbulence();
                  new_gw_data.f1_turb = gw.data.TurbulanceParameter.f_1.value();
                  new_gw_data.f2_turb = gw.data.TurbulanceParameter.f_2.value();
                  new_gw_data.omega_2_turb =
                      gw.data.TurbulanceParameter.Omega_2.value();

                  // SNR of Collision
                  gw.data.collisionON = true;
                  gw.data.swON        = false;
                  gw.data.turbON      = false;
                  new_gw_data.SNR_col = gw.GetSNR(1e-6, 10);

                  // SNR of Sound Waves
                  gw.data.collisionON = false;
                  gw.data.swON        = true;
                  gw.data.turbON      = false;
                  new_gw_data.SNR_sw  = gw.GetSNR(1e-6, 10);

                  // SNR of Turbulence
                  gw.data.collisionON  = false;
                  gw.data.swON         = false;
                  gw.data.turbON       = true;
                  new_gw_data.SNR_turb = gw.GetSNR(1e-6, 10);

                  // SNR of all contributions
                  gw.data.collisionON = true;
                  gw.data.swON        = true;
                  gw.data.turbON      = true;
                  new_gw_data.SNR     = gw.GetSNR(1e-6, 10);

                  new_gw_data.kappa_col    = gw.data.kappa_col;
                  new_gw_data.kappa_sw     = gw.data.kappa_sw;
                  new_gw_data.Epsilon_Turb = gw.data.Epsilon_Turb;
                  new_gw_data.cs_f         = gw.data.Csound_false;
                  new_gw_data.cs_t         = gw.data.Csound_true;

                  new_gw_data.status_gw = gw.data.status;
                }
              }
              else if (input.gw_calculation &&
                       trans_status != BSMPT::StatusTemperature::Success)
              {
                Logger::Write(LoggingLevel::TransitionDetailed,
                              "Requested transition temperature could not be "
                              "calculated.");
              }
            }
            else
            {
              output_store.status.status_nucl_approx.push_back(
                  BSMPT::StatusTemperature::NaN);
              output_store.status.status_nucl.push_back(
                  BSMPT::StatusTemperature::NaN);
              output_store.status.status_perc.push_back(
                  BSMPT::StatusTemperature::NaN);
              output_store.status.status_compl.push_back(
                  BSMPT::StatusTemperature::NaN);

              new_transition_data.nucl_approx_true_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.nucl_approx_false_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.nucl_true_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.nucl_false_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.perc_true_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.perc_false_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.compl_true_vev =
                  std::vector<double>(num_vev, EmptyValue);
              new_transition_data.compl_false_vev =
                  std::vector<double>(num_vev, EmptyValue);
            }
          }
          else
          {
            new_transition_data.crit_true_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.crit_false_vev =
                std::vector<double>(num_vev, EmptyValue);

            output_store.status.status_bounce_sol.push_back(StatusGW::NotSet);

            output_store.status.status_nucl_approx.push_back(
                BSMPT::StatusTemperature::NaN);
            output_store.status.status_nucl.push_back(
                BSMPT::StatusTemperature::NaN);
            output_store.status.status_perc.push_back(
                BSMPT::StatusTemperature::NaN);
            output_store.status.status_compl.push_back(
                BSMPT::StatusTemperature::NaN);

            new_transition_data.nucl_approx_true_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.nucl_approx_false_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.nucl_true_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.nucl_false_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.perc_true_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.perc_false_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.compl_true_vev =
                std::vector<double>(num_vev, EmptyValue);
            new_transition_data.compl_false_vev =
                std::vector<double>(num_vev, EmptyValue);
          }

          new_transition_data.crit_temp = pair.crit_temp;
          output_store.vec_trans_data.push_back(new_transition_data);
          output_store.vec_gw_data.push_back(new_gw_data);
        }

        // transition history evaluator
        std::vector<int> transition_history;
        std::vector<int> pair_history;
        bool final_true_phase_reached = false;
        int tmp_phase_id              = 0; // initial false phase
        double tmp_compl_temp         = -1;

        int tmp_next_phase_id;
        double pair_compl_temp;
        int tmp_pair_id;

        while (not final_true_phase_reached)
        {
          // store current false phase
          transition_history.push_back(tmp_phase_id);
          tmp_next_phase_id = -1;
          tmp_compl_temp    = -1;

          for (auto pair : vec_coex)
          {
            // get pair with matching false phase id
            if (pair.false_phase.id == tmp_phase_id)
            {
              pair_compl_temp =
                  output_store.vec_trans_data.at(pair.coex_pair_id)
                      .compl_temp.value_or(EmptyValue);

              if (std::isnan(pair_compl_temp)) // completion temperature not
                                               // reached in pair
              {
                continue;
              }
              else
              {
                // update next phase id with true phase with highest
                // completion temperature
                if (tmp_compl_temp == -1)
                {
                  tmp_compl_temp    = pair_compl_temp;
                  tmp_next_phase_id = pair.true_phase.id;
                  tmp_pair_id       = pair.coex_pair_id;
                }
                else
                {
                  if (pair_compl_temp > tmp_compl_temp)
                  {
                    tmp_compl_temp    = pair_compl_temp;
                    tmp_next_phase_id = pair.true_phase.id;
                    tmp_pair_id       = pair.coex_pair_id;
                  }
                }
              }
            }
          }

          if (tmp_next_phase_id == -1)
          {
            final_true_phase_reached = true;
          }
          else
          {
            tmp_phase_id      = tmp_next_phase_id;
            tmp_next_phase_id = -1;
            pair_history.push_back(tmp_pair_id);
            tmp_pair_id = -1;
          }
        }


        // CB added v
        // Check if the last reached phase (as in transition_history.back()) corresponds to
        // the correct EW phase at T = 0 with the minimum as given in the model definition
        // (from modelPointer->get_vevTreeMin())
        auto p = vac.PhasesList[transition_history.back()];
        if (p.T_low > 0. ||
            !almost_the_same(input.modelPointer->get_vevTreeMin(),
                             p.Get(0.).point, false, 0.01, 1e-5))
        {
          output_store.status.status_last_phase_ew = StatusLastPhaseEW::Failure;
        }
        else
        {
          output_store.status.status_last_phase_ew = StatusLastPhaseEW::Success;
        }
        // CB added ^


        output_store.transition_history =
            std::to_string(transition_history.at(0));
        if (transition_history.size() > 1)
        {
          for (std::size_t i = 1; i < transition_history.size(); i++)
          {
            output_store.transition_history +=
                "-(" + std::to_string(pair_history.at(i - 1)) + ")->" +
                std::to_string(transition_history.at(i));
          }
        }
      }
    }
  }
  else
  {
    Logger::Write(LoggingLevel::TransitionDetailed, "Point is not NLO stable.");
  }
  return;
}

TransitionTracer::~TransitionTracer()
{
}

// CB: changed to output the mass/temp ratio for scalar and gauge bosons separately
double TransitionTracer::CheckMassRatio(const user_input &input,
                                        const std::vector<double> &vec,
                                        const double &temp) const
{
  auto [m2T2_scalar, m2T2_gauge] = CheckMassRatio_sep(input, vec, temp);
  return std::max(m2T2_scalar, m2T2_gauge);
}

std::pair<double, double> TransitionTracer::CheckMassRatio_sep(const user_input &input,
                                        const std::vector<double> &vec,
                                        const double &temp) const
{
  std::stringstream ss;
  std::vector<double> massOverTempSq, massOverTempSqScalar, massOverTempSqGauge;
  massOverTempSqScalar = input.modelPointer->HiggsMassesSquared(
                       input.modelPointer->MinimizeOrderVEV(vec), temp) /
                   std::pow(temp, 2);
  massOverTempSqGauge = input.modelPointer->GaugeMassesSquared(
                            input.modelPointer->MinimizeOrderVEV(vec), temp) /
                        std::pow(temp, 2);

  massOverTempSq = massOverTempSqScalar;
  massOverTempSq.insert(massOverTempSq.end(),
                        massOverTempSqGauge.begin(),
                        massOverTempSqGauge.end());

  int color = 0;
  for (auto el : massOverTempSq)
  {
    if (el > 0.25) // m/T > 0.5
    {
      color = 1;
      if (el > 1) // m/T > 1.0
      {
        color = 2;
        break;
      }
    }
  }

  if (color == 0)
  {
    ss << "\n\033[1;92mm^2(vev_false, T = " << std::to_string(temp)
       << ") / T^2 = " << massOverTempSq << "\033[0m\n";
  }
  else if (color == 1)
  {
    ss << "\n\033[1;93mm^2(vev_false, T = " << std::to_string(temp)
       << ") / T^2 = " << massOverTempSq << "\033[0m\n";
  }
  else
  {
    ss << "\n\033[1;91mm^2(vev_false, T = " << std::to_string(temp)
       << ") / T^2 = " << massOverTempSq << "\033[0m\n";
  }

  Logger::Write(LoggingLevel::TransitionDetailed, ss.str());
  return std::make_pair(*std::max_element(massOverTempSqScalar.begin(), massOverTempSqScalar.end()),
                        *std::max_element(massOverTempSqGauge.begin(), massOverTempSqGauge.end()));
}

} // namespace BSMPT
