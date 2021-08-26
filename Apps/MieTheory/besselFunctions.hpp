/**
 * Copyright 2019 Sina Haeri
 *
 * Authors: See AUTHORS
 *
 * Contact: [s.haeri@ed.ac.uk and/or si.haeri@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @brief This is a part of MieLAM software package. besselFunctions class
 * implements robust and stable algorithms for calculation of Bessel functions
 * of complex argument and integer degree (with arbitrary precision but this
 * capability is removed here). 
 *
 * @author Sina Haeri (s.haeri@ed.ac.uk)
 * @author Soroosh Haeri (soroosh.haeri@gmail.com)
 */

#ifndef BESSEL_FUNCTIONS_HPP
#define BESSEL_FUNCTIONS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/multiprecision/eigen.hpp>
#include "boostAssertMessage.hpp"

#include <iostream>
#include <complex.h>
#include <exception>

class besselFunctions
{
private:

    typedef Eigen::Array<complexd, Eigen::Dynamic, 1> array1DComplexType;
    const double one = double(1.0);
    const double zero = double(0.0);
    const double two = double(2.0);
    const double mtwo = double(-2.0);
    const double half = double(0.5);
    const complexd c_i = {zero, one};
    const complexd c_one = {one, zero};

    int _maxOrder;
    complexd _argument;
    array1DComplexType _sphBesselj, _sphBessely, _sphBesselk;
    array1DComplexType _riccatiBesselRatio;
    array1DComplexType _riccatiBesselPsi, _riccatiBesselChi,
                       _riccatiBesselXi, _dRiccatiBesselPsi,
                       _dRiccatiBesselXi;

    double _minReal = std::numeric_limits<double>::min();

    static inline int _isOdd(int x) { return x & 1; }

    void _calculateAllSphBessel(int orderFrom)
    {
        using boost::math::epsilon_difference;
        int maxRequired = _maxOrder + 1;
        int nelem = maxRequired - orderFrom + 1;
        int ii,nn;

        complexd oneOArg = one / _argument;
        complexd sphBesselj0 = sin(_argument) * oneOArg;
        complexd firstTerm;
        if ( abs(sphBesselj0) < std::numeric_limits<double>::epsilon() )
        {
            //Assume it is a root 
            //Set j1 as the first term of the continued fractions
            firstTerm = - cos(_argument) * oneOArg;
            if(!orderFrom) _sphBesselj(0) = 0.;
            nn = 2;
        } else
        {
            firstTerm = sphBesselj0;
            if(!orderFrom) _sphBesselj(0) = sphBesselj0;
            nn = 1;
        }
        
        
        for(ii = std::max(1,orderFrom); ii <= maxRequired; ++ii)
        {
            _sphBesselj(ii) = firstTerm * _sphBesseljCFractions(ii, nn, _argument);
        }

        for(ii = orderFrom; ii <= maxRequired; ++ii)
        {
            _riccatiBesselRatio(ii) = double(-ii) / _argument +
                _riccatiBesselRatioCFractions(ii, _argument);
        }

        if(orderFrom)
        {

            _sphBessely.segment(orderFrom, nelem) =
                            ( _sphBesselj.segment(orderFrom, nelem) *
                              _sphBessely.segment(orderFrom-1, nelem) -
                              oneOArg * oneOArg ) /
                            _sphBesselj.segment(orderFrom-1, nelem);
        } else
        {
            complexd _sphBesselym1 = _sphBesselj(0);
            complexd _sphBesseljm1 = _sphBesselj(0) / _argument - _sphBesselj(1);
            _sphBessely(0) = ( _sphBesselj(0) * _sphBesselym1 -
                                oneOArg * oneOArg ) / _sphBesseljm1;

            _sphBessely.segment(1, nelem-1) =
                            ( _sphBesselj.segment(1, nelem-1) *
                              _sphBessely.segment(0, nelem-1) -
                              oneOArg * oneOArg ) /
                            _sphBesselj.segment(0, nelem-1);
        }

        //Defined assuming exp(-iwt) time dependance
        if(orderFrom > 1)
        {
            for(ii = orderFrom; ii <= maxRequired; ++ii)
                _sphBesselk(ii) = (two * ii - one) * _sphBesselk(ii-1) *
                                                    oneOArg - _sphBesselk(ii-2);
        }
        else
        {
            _sphBesselk(0) = sphBesselj0 - c_i * cos(_argument) * oneOArg;
            if(maxRequired >= 1)
            {
                _sphBesselk(1) = _sphBesselk(0) * ( oneOArg - c_i );
                for(ii = 2; ii <= maxRequired; ++ii)
                    _sphBesselk(ii) = (two * ii - one) * _sphBesselk(ii-1) *
                                                oneOArg - _sphBesselk(ii-2);
            }
        }
 
        //_sphBesselk.segment(orderFrom, nelem) =
        //                _sphBesselj.segment(orderFrom, nelem) +
        //                c_i * _sphBessely.segment(orderFrom, nelem);

    }

    complexd _sphBesseljCFractions(int n, int nn, complexd z)
    {
        using boost::math::float_distance;
        complexd numProduct, denProduct, nextTermNum, nextTermDen;
        int ii,jj;

        if(nn==2) return one;

        //Calculate the first n terms
        nextTermNum = _bnm(nn, 1, z);
        numProduct = nextTermNum;
        for( ii =2; ii <= n+2-nn; ++ii )
        {
            nextTermNum = _bnm(nn, ii, z) - one/nextTermNum;
            numProduct *= nextTermNum;
        }

        ii = n+2-nn;
        jj = 2;
        nextTermDen = _bnm(n, jj, z);
        denProduct = nextTermDen;
        do {

            ii++;
            nextTermNum = _bnm(nn, ii, z) - one/nextTermNum;
            numProduct *= nextTermNum;

            jj++;
            nextTermDen = _bnm(n, jj, z) - one/nextTermDen;
            denProduct *= nextTermDen;
            
        } while(fabs(float_distance(nextTermNum.imag(), nextTermDen.imag())) > zero ||
                fabs(float_distance(nextTermNum.real(), nextTermDen.real())) > zero   );

        return denProduct / numProduct;
    }

    complexd _riccatiBesselRatioCFractions(int n, complexd z)
    {
        using boost::math::float_distance;
        complexd numProduct, denProduct, nextTermNum, nextTermDen;
        double nu;
        int ii;

        nu = n + half;
        ii = 1;
        nextTermNum = _anun(ii, nu, z);
        numProduct = nextTermNum;

        ii = 2;
        nextTermDen = _anun(ii, nu, z);
        denProduct = nextTermDen;

        for(;;)
        {
            nextTermNum = _anun(ii, nu, z) + one/nextTermNum;
            numProduct *= nextTermNum;

            ii++;
            if( abs(float_distance(abs(nextTermNum.imag()),
                abs(nextTermDen.imag()))) < one &&
                abs(float_distance(abs(nextTermNum.real()),
                abs(nextTermDen.real()))) < one ) break;

            nextTermDen = _anun(ii, nu, z) + one/nextTermDen;
            denProduct *= nextTermDen;

        }

        return numProduct / denProduct;
    }

    inline complexd _bnm(int n, int m, complexd z)
    {
        return two * (n + m - half) / z;
    }

    inline complexd _anun(int n, double nu, complexd z)
    {
        if(_isOdd(n+1))
        {
            return mtwo * (nu + n - one) / z;
        }
        else
        {
            return two * (nu + n - one) / z;
        }

    }

public:

    besselFunctions() :
        _maxOrder(0), _argument(one),
        _sphBesselj(array1DComplexType::Constant(2, 1, zero)),
        _sphBessely(array1DComplexType::Constant(2, 1, zero)),
        _sphBesselk(array1DComplexType::Constant(2, 1, zero)),
        _riccatiBesselRatio(array1DComplexType::Constant(2, 1, zero)),
        _riccatiBesselPsi(array1DComplexType::Constant(1, 1, zero)),
        _riccatiBesselChi(array1DComplexType::Constant(1, 1, zero)),
        _riccatiBesselXi(array1DComplexType::Constant(1, 1, zero)),
        _dRiccatiBesselPsi(array1DComplexType::Constant(1, 1, zero)),
        _dRiccatiBesselXi(array1DComplexType::Constant(1, 1, zero))
    {
        _calculateAllSphBessel(0);
    }

    besselFunctions(int n, complexd z)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");
        BOOST_ASSERT_MSG(abs(z.imag()) >= _minReal || abs(z.real()) >= _minReal,
                "The besselFunctions class is only defined for abs(z) > 0");

        _maxOrder = n;
        _argument = z;
        _sphBesselj = array1DComplexType::Constant(n+2, 1, zero);
        _sphBessely = array1DComplexType::Constant(n+2, 1, zero);
        _sphBesselk = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselRatio = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselPsi = array1DComplexType::Constant(n, 1, zero);
        _riccatiBesselChi = array1DComplexType::Constant(n, 1, zero);
        _riccatiBesselXi = array1DComplexType::Constant(n, 1, zero);
        _dRiccatiBesselPsi = array1DComplexType::Constant(n, 1, zero);
        _dRiccatiBesselXi = array1DComplexType::Constant(n, 1, zero);

        _calculateAllSphBessel(0);
    }

    besselFunctions(int n) :
    _argument(one)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");

        _maxOrder = n;
        _sphBesselj = array1DComplexType::Constant(n+2, 1, zero);
        _sphBessely = array1DComplexType::Constant(n+2, 1, zero);
        _sphBesselk = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselRatio = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselPsi = array1DComplexType::Constant(n, 1, zero);
        _riccatiBesselChi = array1DComplexType::Constant(n, 1, zero);
        _riccatiBesselXi = array1DComplexType::Constant(n, 1, zero);
        _dRiccatiBesselPsi = array1DComplexType::Constant(n, 1, zero);
        _dRiccatiBesselXi = array1DComplexType::Constant(n, 1, zero);
        _calculateAllSphBessel(0);
    }

    besselFunctions(const besselFunctions& other) :
        _maxOrder( other._maxOrder ),
        _argument( other._argument ),
        _sphBesselj( other._sphBesselj ),
        _sphBessely( other._sphBessely ),
        _sphBesselk( other._sphBesselk ),
        _riccatiBesselRatio( other._riccatiBesselRatio ),
        _riccatiBesselPsi(other._riccatiBesselPsi),
        _riccatiBesselChi(other._riccatiBesselChi),
        _riccatiBesselXi(other._riccatiBesselXi),
        _dRiccatiBesselPsi(other._dRiccatiBesselPsi),
        _dRiccatiBesselXi(other._dRiccatiBesselXi)
    {}

    besselFunctions& operator=(const besselFunctions& other)
    {
        if (this != &other)
        {
            _argument = other._argument;
            _maxOrder = other._maxOrder;
            _sphBesselj = other._sphBesselj;
            _sphBessely = other._sphBessely;
            _sphBesselk = other._sphBesselk;
            _riccatiBesselRatio = other._riccatiBesselRatio;
            _riccatiBesselPsi = other._riccatiBesselPsi;
            _riccatiBesselChi = other._riccatiBesselChi;
            _riccatiBesselXi  = other._riccatiBesselXi;
            _dRiccatiBesselPsi= other._dRiccatiBesselPsi;
            _dRiccatiBesselXi = other._dRiccatiBesselXi;
        }
        return *this;
    }

    void setArgument(complexd z)
    {
        BOOST_ASSERT_MSG(abs(z.imag()) >= _minReal || abs(z.real()) >= _minReal,
                "The besselFunctions class is only defined for abs(z) > 0");
        _argument = z;
        _calculateAllSphBessel(0);

    }

    void setMaxOrder(int n)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");

        if( n > _maxOrder )
        {
            int oldOrder = _maxOrder;
            _maxOrder = n;
            _sphBesselj.conservativeResize(_maxOrder+2);
            _sphBessely.conservativeResize(_maxOrder+2);
            _sphBesselk.conservativeResize(_maxOrder+2);
            _riccatiBesselRatio.conservativeResize(_maxOrder+2);
            _riccatiBesselPsi.conservativeResize(_maxOrder);
            _riccatiBesselChi.conservativeResize(_maxOrder);
            _riccatiBesselXi.conservativeResize(_maxOrder);
            _dRiccatiBesselPsi.conservativeResize(_maxOrder);
            _dRiccatiBesselXi.conservativeResize(_maxOrder);
            _calculateAllSphBessel(oldOrder+1);

        } else if( n < _maxOrder )
        {
            _maxOrder = n;
            _sphBesselj.conservativeResize(_maxOrder+2);
            _sphBessely.conservativeResize(_maxOrder+2);
            _sphBesselk.conservativeResize(_maxOrder+2);
            _riccatiBesselRatio.conservativeResize(_maxOrder+2);
            _riccatiBesselPsi.conservativeResize(_maxOrder);
            _riccatiBesselChi.conservativeResize(_maxOrder);
            _riccatiBesselXi.conservativeResize(_maxOrder);
            _dRiccatiBesselPsi.conservativeResize(_maxOrder);
            _dRiccatiBesselXi.conservativeResize(_maxOrder);
        }
    }

    /*const array1DComplexType& getSphBesselj() const
    {
        return _sphBesselj.segment(1,_maxOrder);
    }

    const array1DComplexType& getSphBessely() const
    {
        return _sphBessely.segment(1,_maxOrder);
    }

    const array1DComplexType& getSphBesselk() const
    {
        return _sphBesselk.segment(1,_maxOrder);
    }

    array1DComplexType getRiccatiBesselRatio() const 
    {
        return _riccatiBesselRatio.segment(1,_maxOrder);
    }*/

    void computeRiccatiBesselPsi(int startFrom)
    {
        BOOST_ASSERT_MSG(startFrom == 1 || !startFrom,
            "The computeRiccatiBesselPsi function only takes 0 and 1");
        _riccatiBesselPsi = _argument * 
                            _sphBesselj.segment( startFrom , _maxOrder );
        std::cout.precision(17);
        std::cout << "In computeRiccatiBesselPsi, _argument = " << _argument << std::endl;
        std::cout << _sphBesselj << std::endl;
        std::cout << _riccatiBesselPsi << std::endl;
    }

    array1DComplexType getRiccatiBesselPsi() const { return _riccatiBesselPsi; }

    void computeRiccatiBesselChi()
    {
        _riccatiBesselChi = _argument * _sphBessely.segment(1,_maxOrder);
    }

    const array1DComplexType& getRiccatiBesselChi() const { return _riccatiBesselChi; }

    void computeRiccatiBesselXi(int startFrom) 
    {
        BOOST_ASSERT_MSG(startFrom == 1 || !startFrom,
            "The computeRiccatiBesselXi function only takes 0 and 1");
        _riccatiBesselXi = _argument * 
                           _sphBesselk.segment( startFrom , _maxOrder );
    }

    const array1DComplexType& getRiccatiBesselXi() const { return _riccatiBesselXi; }

    void computeDRiccatiBesselPsi() 
    {
        _dRiccatiBesselPsi = half * 
                            ( _argument * ( _sphBesselj.segment(0,_maxOrder) -
                                            _sphBesselj.segment(2,_maxOrder) ) +
                              _sphBesselj.segment(1,_maxOrder) );
    }

    const array1DComplexType& getDRiccatiBesselPsi() const { return _dRiccatiBesselPsi; }

    void computeDRiccatiBesselXi()
    {
        _dRiccatiBesselXi = half * 
                            ( _argument * ( _sphBesselk.segment(0,_maxOrder) -
                                            _sphBesselk.segment(2,_maxOrder) ) +
                              _sphBesselk.segment(1,_maxOrder) );
    }

    const array1DComplexType& getDRiccatiBesselXi() const { return _dRiccatiBesselXi; }
};
#endif //BESSEL_FUNCTIONS_HPP
