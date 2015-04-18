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

#include <armadillo>

using namespace arma;

#ifndef _Z_TCF_HPP_
#define _Z_TCF_HPP_

template <class T>
class TCF
{
    private:
        unsigned length_, interval_, zeros_, number_, speclength_;
        unsigned corr_int_, corr_, num_corr_;
        Row<T> correlation_function_;

    public:
        TCF(unsigned length, unsigned interval = 1, unsigned zeros = 0):
            length_(length), interval_(interval), zeros_(zeros)
        {
            correlation_function_.zeros(length_+zeros_);
        }
        void correlate(const mat& vec1, const mat& vec2);
        void correlate(const mat& vec1, const mat& vec2,
            const unsigned mod);
        void correlate_one_direction(const mat& vec1, const mat& vec2);
        void correlate_one_direction(const mat& vec1, const mat& vec2,
            const unsigned mod);
        void correlate_1D(const rowvec& vec1, const rowvec& vec2);
        void correlate_1D(const rowvec& vec1, const rowvec& vec2,
            const unsigned mod);
        void correlate_1D_one_direction(const rowvec& vec1, const rowvec& vec2);
        void correlate_1D_one_direction(const rowvec& vec1, const rowvec& vec2,
            const unsigned mod);
        inline void correlate_one_direction(const rowvec& vec)
        {
            correlate_one_direction(vec, vec);
        }
        inline void correlate_one_direction(const rowvec& vec,
            const unsigned mod)
        {
            correlate_one_direction(vec, vec, mod);
        }
        inline void correlate_1D_one_direction(const rowvec& vec)
        {
            correlate_1D_one_direction(vec, vec);
        }
        inline void correlate_1D_one_direction(const rowvec& vec,
            const unsigned mod)
        {
            correlate_1D_one_direction(vec, vec, mod);
        }
        unsigned length() { return length_; }
        inline T tcf(int i) { return correlation_function_(i); }
};
#endif
