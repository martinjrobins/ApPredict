/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifdef CHASTE_CVODE

#ifndef _TESTPOPULATIONOFMODELS_HPP_
#define _TESTPOPULATIONOFMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>

#include "AbstractCvodeCell.hpp"
#include "ApPredictMethods.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SetupModel.hpp"

/**
 * Tests for the "population of models" feature: one or more oxmeta parameters are given a set
 * of pre-sampled values in a file (--parameter-samples-file), and each concentration is
 * simulated for every sampled cell to produce a distribution of results.
 */
class TestPopulationOfModels : public CxxTest::TestSuite
{
private:
    /** Resolve the herg (IKr) conductance parameter name for a model, coping with the
     *  '..._scaling_factor' variant. */
    std::string ResolveHergName(boost::shared_ptr<AbstractCvodeCell> pModel)
    {
        std::string name = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
        if (pModel->HasParameter(name))
        {
            return name;
        }
        return name + "_scaling_factor";
    }

    /** Write a samples file scattering the herg conductance around its default, and return the
     *  path (and the parameter name written in the header). */
    std::string WriteHergSamplesFile(OutputFileHandler& rHandler, const std::string& rFileName)
    {
        SetupModel setup(1.0, 2u); // ten Tusscher 2006 epi
        boost::shared_ptr<AbstractCvodeCell> p_model = setup.GetModel();
        std::string herg_name = ResolveHergName(p_model);
        double default_value = p_model->GetParameter(herg_name);

        const double factors[5] = { 0.7, 0.85, 1.0, 1.15, 1.3 };
        out_stream p_file = rHandler.OpenOutputFile(rFileName);
        (*p_file) << herg_name << "\n";
        for (unsigned i = 0; i < 5; i++)
        {
            (*p_file) << default_value * factors[i] << "\n";
        }
        p_file->close();
        return rHandler.GetOutputDirectoryFullPath() + rFileName;
    }

public:
    /**
     * With no drug, every cell's APD90 is concentration-independent, so the per-sample change
     * in APD90 must be exactly zero at every concentration EVEN THOUGH the raw APD90 spreads
     * across the population. This is the key check that deltas are paired per-sample 
     * (a scalar 'raw percentile minus a single control' would wrongly report a spread of deltas).
     */
    void TestNoDrugDeltasAreZeroButRawApd90Spreads()
    {
        OutputFileHandler handler("TestPopulationOfModels_NoDrug");
        std::string samples_path = WriteHergSamplesFile(handler, "herg_samples.txt");

        std::stringstream args;
        args << "--model 2 --plasma-concs 0 10 --pacing-max-time 0.1 --credible-intervals "
             << "--parameter-samples-file " << samples_path;
        CommandLineArgumentsMocker wrapper(args.str());

        ApPredictMethods methods;
        methods.SetOutputDirectory("TestPopulationOfModels_NoDrug_output/");
        methods.Run();

        std::vector<double> concs = methods.GetConcentrations();
        const unsigned num_concs = concs.size();
        TS_ASSERT_LESS_THAN_EQUALS(2u, num_concs);

        // Central line: one APD90 per concentration.
        std::vector<double> apd90s = methods.GetApd90s();
        TS_ASSERT_EQUALS(apd90s.size(), num_concs);

        // Raw APD90 credible regions: population spread, so lower < upper at each concentration.
        std::vector<std::vector<double> > raw = methods.GetApd90CredibleRegions();
        TS_ASSERT_EQUALS(raw.size(), num_concs);
        for (unsigned c = 0; c < raw.size(); c++)
        {
            TS_ASSERT_EQUALS(raw[c].size(), 2u); // default 95% interval => 2 percentiles
            TS_ASSERT_LESS_THAN(raw[c][0], raw[c][1]);
        }

        // Delta APD90 credible regions: per-sample, so exactly zero everywhere with no drug
        std::vector<std::vector<double> > deltas = methods.GetDeltaApd90CredibleRegions();
        TS_ASSERT_EQUALS(deltas.size(), num_concs);
        for (unsigned c = 0; c < deltas.size(); c++)
        {
            TS_ASSERT_EQUALS(deltas[c].size(), 2u);
            TS_ASSERT_DELTA(deltas[c][0], 0.0, 1e-6);
            TS_ASSERT_DELTA(deltas[c][1], 0.0, 1e-6);
        }
    }

    /**
     * A drug (median dose-response) may be applied on top of the population. Blocking hERG
     * prolongs the AP, so the per-sample delta APD90 should become clearly positive at the
     * non-zero concentration, while the control delta stays zero.
     */
    void TestDrugComposesWithPopulation()
    {
        OutputFileHandler handler("TestPopulationOfModels_Drug");
        std::string samples_path = WriteHergSamplesFile(handler, "herg_samples.txt");

        std::stringstream args;
        args << "--model 2 --plasma-concs 0 30 --pacing-max-time 1 --credible-intervals "
             << "--pic50-herg 5 "
             << "--parameter-samples-file " << samples_path;
        CommandLineArgumentsMocker wrapper(args.str());

        ApPredictMethods methods;
        methods.SetOutputDirectory("TestPopulationOfModels_Drug_output/");
        methods.Run();

        std::vector<double> concs = methods.GetConcentrations();
        const unsigned num_concs = concs.size();
        std::vector<std::vector<double> > deltas = methods.GetDeltaApd90CredibleRegions();
        TS_ASSERT_EQUALS(deltas.size(), num_concs);

        // Control (0 uM) delta is still exactly zero.
        TS_ASSERT_DELTA(deltas[0][0], 0.0, 1e-6);
        TS_ASSERT_DELTA(deltas[0][1], 0.0, 1e-6);

        // At the highest concentration, hERG block prolongs the AP: the delta APD90 credible
        // region should be clearly positive.
        const unsigned last = num_concs - 1u;
        TS_ASSERT_LESS_THAN(0.5, deltas[last][0]);
        TS_ASSERT_LESS_THAN(0.0, deltas[last][1]);
    }

    /**
     * Guard-rails: drug uncertainty cannot be combined with a population of models, and the
     * named parameters must exist on the model, and a control concentration is required.
     */
    void TestGuardRails()
    {
        OutputFileHandler handler("TestPopulationOfModels_Guards");
        std::string samples_path = WriteHergSamplesFile(handler, "herg_samples.txt");

        // Brute force + samples file -> throws.
        {
            std::stringstream args;
            args << "--model 2 --plasma-concs 0 10 --pacing-max-time 0.1 --credible-intervals "
                 << "--brute-force 10 --parameter-samples-file " << samples_path;
            CommandLineArgumentsMocker wrapper(args.str());
            ApPredictMethods methods;
            methods.SetOutputDirectory("TestPopulationOfModels_Guards_output/");
            TS_ASSERT_THROWS_CONTAINS(methods.Run(), "cannot be combined with --brute-force");
        }

        // Drug spread + samples file -> throws.
        {
            std::stringstream args;
            args << "--model 2 --plasma-concs 0 10 --pacing-max-time 0.1 --credible-intervals "
                 << "--pic50-herg 5 --pic50-spread-herg 0.15 --parameter-samples-file " << samples_path;
            CommandLineArgumentsMocker wrapper(args.str());
            ApPredictMethods methods;
            methods.SetOutputDirectory("TestPopulationOfModels_Guards_output/");
            TS_ASSERT_THROWS_CONTAINS(methods.Run(), "cannot be combined with drug 'spread'");
        }

        // A parameter the model does not have -> throws.
        {
            out_stream p_file = handler.OpenOutputFile("bad_param.txt");
            (*p_file) << "not_a_real_oxmeta_parameter\n1.0\n2.0\n";
            p_file->close();
            std::string bad_path = handler.GetOutputDirectoryFullPath() + "bad_param.txt";

            std::stringstream args;
            args << "--model 2 --plasma-concs 0 10 --pacing-max-time 0.1 "
                 << "--parameter-samples-file " << bad_path;
            CommandLineArgumentsMocker wrapper(args.str());
            ApPredictMethods methods;
            methods.SetOutputDirectory("TestPopulationOfModels_Guards_output/");
            TS_ASSERT_THROWS_CONTAINS(methods.Run(), "does not have 'not_a_real_oxmeta_parameter'");
        }
    }
};

#endif //_TESTPOPULATIONOFMODELS_HPP_

#endif //CHASTE_CVODE
