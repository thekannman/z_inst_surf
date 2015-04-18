//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

#include <fstream>
#include <armadillo>
#include <string>

using namespace arma;

#ifndef _Z_HISTOGRAM_HPP_
#define _Z_HISTOGRAM_HPP_

class Hist
{
    private:
        unsigned size, total;
        urowvec array;
        double min, max, interval;
        std::string filename_;
        std::ofstream file_;

    public:
        inline void set_filename(const std::string& filename)
        {
            filename_ = filename;
        }
        void print(const unsigned numBins, const double minBin, const double binSize,
            const double norm);
        inline void print(const std::string& filename, const unsigned num_bins,
            const double min_bin, const double bin_size, const double norm)
        {
            set_filename(filename);
            print(num_bins, min_bin, bin_size, norm);
        }
};

#endif
