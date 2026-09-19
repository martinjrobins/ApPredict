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

#ifndef PARAMETERPOINTSAMPLESREADER_HPP_
#define PARAMETERPOINTSAMPLESREADER_HPP_

#include <cmath>

#include "AbstractDataStructure.hpp"

/**
 * Helper class to read in a "population of models" parameter samples file.
 *
 * The file format is:
 *  * a header line listing the oxmeta parameter names (space, tab or comma separated),
 *  * followed by one row per sample, each row giving one (absolute) value per named
 *    parameter in the same column order as the header.
 *
 * Row i across all parameters defines one virtual cell (paired rows), so every row must
 * have exactly as many values as there are parameter names in the header.
 *
 * Values are stored internally as [parameter][sample] to match the shape used elsewhere
 * in ApPredict for sampled quantities (e.g. mSampledIc50s).
 */
class ParameterPointSamplesReader : public AbstractDataStructure
{
private:
    /** The oxmeta parameter names, in the order they appear in the header line. */
    std::vector<std::string> mParameterNames;

    /**
     * The sampled values.
     * The outer index is over parameters (matching #mParameterNames),
     * the inner index is over samples (one entry per data row in the file).
     */
    std::vector<std::vector<double> > mSamples;

    /** The number of samples (data rows) read from the file. */
    unsigned mNumSamples;

protected:
    /**
     * Read the header line of parameter names.
     *
     * @param rLine  the header line in stringstream format.
     * @return true (this class always has a header line).
     */
    bool LoadHeaderLine(std::stringstream& rLine)
    {
        std::string name;
        while (rLine >> name)
        {
            mParameterNames.push_back(name);
        }

        if (mParameterNames.empty())
        {
            EXCEPTION("No parameter names found on the header line of the parameter samples file.");
        }

        // Size the outer vector now we know how many parameters there are.
        mSamples.resize(mParameterNames.size());
        return true;
    }

    /**
     * Read a single sample (one row of the file).
     *
     * @param rLine  a data line in stringstream format.
     */
    void LoadALine(std::stringstream& rLine)
    {
        for (unsigned i = 0; i < mParameterNames.size(); i++)
        {
            double value;
            rLine >> value;
            if (rLine.fail() || !std::isfinite(value))
            {
                EXCEPTION("Could not read a finite value for every parameter on sample row "
                          << (mNumSamples + 1u) << " of the parameter samples file (expected "
                          << mParameterNames.size() << " numeric values per row).");
            }
            mSamples[i].push_back(value);
        }
        mNumSamples++;
    }

public:
    /**
     * Constructor.
     *
     * @param rFileFinder  a FileFinder pointing at the parameter samples file.
     */
    ParameterPointSamplesReader(FileFinder& rFileFinder)
            : AbstractDataStructure(),
              mNumSamples(0u)
    {
        LoadDataFromFile(rFileFinder.GetAbsolutePath(), 1u);

        if (mNumSamples == 0u)
        {
            EXCEPTION("No sample rows found in the parameter samples file '"
                      << rFileFinder.GetAbsolutePath() << "' (only a header line was present).");
        }
    }

    /**
     * @return the oxmeta parameter names in file (column) order.
     */
    const std::vector<std::string>& rGetParameterNames() const
    {
        return mParameterNames;
    }

    /**
     * @return the number of parameters (columns).
     */
    unsigned GetNumParameters() const
    {
        return mParameterNames.size();
    }

    /**
     * @return the number of samples (data rows).
     */
    unsigned GetNumSamples() const
    {
        return mNumSamples;
    }

    /**
     * @param paramIndex  the index of the parameter (matching rGetParameterNames()).
     * @return the sampled values for this parameter, indexed by sample.
     */
    const std::vector<double>& rGetSamplesForParameter(unsigned paramIndex) const
    {
        return mSamples[paramIndex];
    }

    /**
     * @param sampleIndex  the index of the sample (data row).
     * @return the values for this sample, indexed by parameter (matching rGetParameterNames()).
     */
    std::vector<double> GetSampleRow(unsigned sampleIndex) const
    {
        std::vector<double> row(mParameterNames.size());
        for (unsigned i = 0; i < mParameterNames.size(); i++)
        {
            row[i] = mSamples[i][sampleIndex];
        }
        return row;
    }
};

#endif // PARAMETERPOINTSAMPLESREADER_HPP_
