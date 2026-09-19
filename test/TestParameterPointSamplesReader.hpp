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

#ifndef _TESTPARAMETERPOINTSAMPLESREADER_HPP_
#define _TESTPARAMETERPOINTSAMPLESREADER_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "ParameterPointSamplesReader.hpp"

/**
 * Tests for the ParameterPointSamplesReader, which reads a "population of models"
 * parameter samples file (header of oxmeta names, then one row per sample).
 */
class TestParameterPointSamplesReader : public CxxTest::TestSuite
{
private:
    /**
     * Helper to write a file with the given contents and return a FileFinder to it.
     */
    FileFinder WriteFile(OutputFileHandler& rHandler, const std::string& rFileName, const std::string& rContents)
    {
        out_stream p_file = rHandler.OpenOutputFile(rFileName);
        (*p_file) << rContents;
        p_file->close();
        return FileFinder(rHandler.GetOutputDirectoryFullPath() + rFileName, RelativeTo::Absolute);
    }

public:
    void TestReadingAGoodFile()
    {
        OutputFileHandler handler("TestParameterPointSamplesReader");

        // Two parameters, three samples. Use a mixture of spaces, tabs and commas.
        std::string contents = "membrane_fast_sodium_current_conductance\tmembrane_rapid_delayed_rectifier_potassium_current_conductance\n"
                               "1.0 2.0\n"
                               "1.1,2.2\n"
                               "0.9\t1.8\n";
        FileFinder file = WriteFile(handler, "good.txt", contents);

        ParameterPointSamplesReader reader(file);

        TS_ASSERT_EQUALS(reader.GetNumParameters(), 2u);
        TS_ASSERT_EQUALS(reader.GetNumSamples(), 3u);

        const std::vector<std::string>& r_names = reader.rGetParameterNames();
        TS_ASSERT_EQUALS(r_names.size(), 2u);
        TS_ASSERT_EQUALS(r_names[0], "membrane_fast_sodium_current_conductance");
        TS_ASSERT_EQUALS(r_names[1], "membrane_rapid_delayed_rectifier_potassium_current_conductance");

        // Values indexed [parameter][sample]
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(0)[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(0)[1], 1.1, 1e-12);
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(0)[2], 0.9, 1e-12);
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(1)[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(1)[1], 2.2, 1e-12);
        TS_ASSERT_DELTA(reader.rGetSamplesForParameter(1)[2], 1.8, 1e-12);

        // Row accessor indexed [parameter]
        std::vector<double> row1 = reader.GetSampleRow(1);
        TS_ASSERT_EQUALS(row1.size(), 2u);
        TS_ASSERT_DELTA(row1[0], 1.1, 1e-12);
        TS_ASSERT_DELTA(row1[1], 2.2, 1e-12);
    }

    void TestRaggedRowThrows()
    {
        OutputFileHandler handler("TestParameterPointSamplesReader", false);

        // Second row has too few values.
        std::string contents = "paramA paramB\n"
                               "1.0 2.0\n"
                               "1.1\n";
        FileFinder file = WriteFile(handler, "ragged_short.txt", contents);
        TS_ASSERT_THROWS_CONTAINS(ParameterPointSamplesReader reader(file),
                                  "Could not read a finite value for every parameter");

        // A row with too many values (the base class catches trailing unread items).
        std::string contents2 = "paramA paramB\n"
                                "1.0 2.0 3.0\n";
        FileFinder file2 = WriteFile(handler, "ragged_long.txt", contents2);
        TS_ASSERT_THROWS_CONTAINS(ParameterPointSamplesReader reader(file2),
                                  "unread items");
    }

    void TestHeaderOnlyThrows()
    {
        OutputFileHandler handler("TestParameterPointSamplesReader", false);
        std::string contents = "paramA paramB\n";
        FileFinder file = WriteFile(handler, "header_only.txt", contents);
        TS_ASSERT_THROWS_CONTAINS(ParameterPointSamplesReader reader(file),
                                  "No sample rows found");
    }

    void TestNonNumericThrows()
    {
        OutputFileHandler handler("TestParameterPointSamplesReader", false);
        std::string contents = "paramA paramB\n"
                               "1.0 not_a_number\n";
        FileFinder file = WriteFile(handler, "non_numeric.txt", contents);
        TS_ASSERT_THROWS_CONTAINS(ParameterPointSamplesReader reader(file),
                                  "Could not read a finite value for every parameter");
    }
};

#endif //_TESTPARAMETERPOINTSAMPLESREADER_HPP_
