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
 * of complex argument and integer degree with arbitrary precision. 
 *
 * @author Sina Haeri (s.haeri@ed.ac.uk)
 * @author Soroosh Haeri (soroosh.haeri@gmail.com)
 */

#ifndef BESSEL_FUNCTIONS_HPP
#define BESSEL_FUNCTIONS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/multiprecision/eigen.hpp>
#include "boostAssertMessage.hpp"

#include <iostream>
#include <exception>

template<typename intType, typename floatType, typename complexType>
class besselFunctions
{
private:

    typedef Eigen::Array<complexType, Eigen::Dynamic, 1> array1DComplexType;
    const complexType c_i {"0.0e+00", "1.0e+00"};
    const complexType c_one {"1.0e+00", "0.0e+00"};
    const floatType one = floatType("1.0");
    const floatType zero = floatType("0.0");
    const floatType two = floatType("2.0");
    const floatType mtwo = floatType("-2.0");
    const floatType half = floatType("0.5");

    intType _maxOrder;
    complexType _argument;
    array1DComplexType _sphBesselj, _sphBessely, _sphBesselk;
    array1DComplexType _riccatiBesselRatio;
    floatType _minReal = std::numeric_limits<floatType>::min();

    void _calculateAllSphBessel(intType orderFrom)
    {
        intType maxRequired = _maxOrder + 1;
        intType nelem = maxRequired - orderFrom + 1;
        intType ii;

        complexType oneOArg = one / _argument;
        complexType sphBesselj0 = sin(_argument) * oneOArg;

        for(ii = orderFrom; ii <= maxRequired; ++ii)
        {
            _sphBesselj(ii) = sphBesselj0 * _sphBesseljCFractions(ii, _argument);
        }

        for(ii = orderFrom; ii <= maxRequired; ++ii)
        {
            _riccatiBesselRatio(ii) = floatType(-ii) / _argument +
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
            complexType _sphBesselym1 = _sphBesselj(0);
            complexType _sphBesseljm1 = _sphBesselj(0) / _argument - _sphBesselj(1);
            _sphBessely(0) = ( _sphBesselj(0) * _sphBesselym1 -
                                oneOArg * oneOArg ) / _sphBesseljm1;

            _sphBessely.segment(1, nelem-1) =
                            ( _sphBesselj.segment(1, nelem-1) *
                              _sphBessely.segment(0, nelem-1) -
                              oneOArg * oneOArg ) /
                            _sphBesselj.segment(0, nelem-1);
        }

        if(orderFrom > 1)
        {
            for(ii = orderFrom; ii <= maxRequired; ++ii)
                _sphBesselk(ii) = (two * ii - one) * _sphBesselk(ii-1) *
                                                    oneOArg - _sphBesselk(ii-2);
        }
        else
        {
            _sphBesselk(0) = sphBesselj0 + c_i * cos(_argument) * oneOArg;
            if(maxRequired >= 1)
            {
                _sphBesselk(1) = _sphBesselk(0) * ( oneOArg + c_i );
                for(ii = 2; ii <= maxRequired; ++ii)
                    _sphBesselk(ii) = (two * ii - one) * _sphBesselk(ii-1) *
                                                oneOArg - _sphBesselk(ii-2);
            }
        }
    }

    complexType _sphBesseljCFractions(intType n, complexType z)
    {
        using boost::math::float_distance;
        complexType numProduct, denProduct, nextTermNum, nextTermDen;
        intType ii,jj;

        if(!n) return one;

        //Calculate the first n terms
        nextTermNum = _bnm(1, 1, z);
        numProduct = nextTermNum;
        for( ii =2; ii <= n+1; ++ii )
        {
            nextTermNum = _bnm(1, ii, z) - one/nextTermNum;
            numProduct *= nextTermNum;
        }

        ii = n+1;
        jj = 2;
        nextTermDen = _bnm(n, jj, z);
        denProduct = nextTermDen;
        do {

            ii++;
            nextTermNum = _bnm(1, ii, z) - one/nextTermNum;
            numProduct *= nextTermNum;

            jj++;
            nextTermDen = _bnm(n, jj, z) - one/nextTermDen;
            denProduct *= nextTermDen;

        } while(abs(float_distance(nextTermNum.imag(), nextTermDen.imag())) > zero ||
                abs(float_distance(nextTermNum.real(), nextTermDen.real())) > zero);
        return denProduct / numProduct;
    }

    complexType _riccatiBesselRatioCFractions(intType n, complexType z)
    {
        using boost::math::float_distance;
        complexType numProduct, denProduct, nextTermNum, nextTermDen;
        floatType nu;
        intType ii;

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

    inline complexType _bnm(intType n, intType m, complexType z)
    {
        return two * (n + m - half) / z;
    }

    inline complexType _anun(intType n, floatType nu, complexType z)
    {
        if(functions::isOdd(n+1))
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
        _argument(one), _maxOrder(0),
        _sphBesselj(array1DComplexType::Constant(2, 1, zero)),
        _sphBessely(array1DComplexType::Constant(2, 1, zero)),
        _sphBesselk(array1DComplexType::Constant(2, 1, zero)),
        _riccatiBesselRatio(array1DComplexType::Constant(2, 1, zero))
    {
        _calculateAllSphBessel(0);
    }

    besselFunctions(intType n, complexType z)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");
        BOOST_ASSERT_MSG(abs(z.imag()) >= _minReal || abs(z.real()) >= _minReal,
                "The besselFunctions class is only defined for abs(z) > 0");

        _argument = z;
        _maxOrder = n;
        _sphBesselj = array1DComplexType::Constant(n+2, 1, zero);
        _sphBessely = array1DComplexType::Constant(n+2, 1, zero);
        _sphBesselk = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselRatio = array1DComplexType::Constant(n+2, 1, zero);

        _calculateAllSphBessel(0);
    }

    besselFunctions(intType n) :
    _argument(one)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");

        _maxOrder = n;
        _sphBesselj = array1DComplexType::Constant(n+2, 1, zero);
        _sphBessely = array1DComplexType::Constant(n+2, 1, zero);
        _sphBesselk = array1DComplexType::Constant(n+2, 1, zero);
        _riccatiBesselRatio = array1DComplexType::Constant(n+2, 1, zero);

        _calculateAllSphBessel(0);
    }

    besselFunctions(const besselFunctions& other) :
        _argument( other._argument ),
        _maxOrder( other._maxOrder ),
        _sphBesselj( other._sphBesselj ),
        _sphBessely( other._sphBessely ),
        _sphBesselk( other._sphBesselk ),
        _riccatiBesselRatio( other._riccatiBesselRatio )
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
        }
        return *this;
    }

    void setArgument(complexType z)
    {
        BOOST_ASSERT_MSG(abs(z.imag()) >= _minReal || abs(z.real()) >= _minReal,
                "The besselFunctions class is only defined for abs(z) > 0");
        _argument = z;
        _calculateAllSphBessel(0);

    }

    void setMaxOrder(intType n)
    {
        BOOST_ASSERT_MSG(n >= 0,
            "The besselFunctions class is only defined for n >= 0.");

        if( n > _maxOrder )
        {
            intType oldOrder = _maxOrder;
            _maxOrder = n;
            _sphBesselj.conservativeResize(_maxOrder+2);
            _sphBessely.conservativeResize(_maxOrder+2);
            _sphBesselk.conservativeResize(_maxOrder+2);
            _riccatiBesselRatio.conservativeResize(_maxOrder+2);

            _calculateAllSphBessel(oldOrder+1);

        } else if( n < _maxOrder )
        {
            _maxOrder = n;
            _sphBesselj.conservativeResize(_maxOrder+2);
            _sphBessely.conservativeResize(_maxOrder+2);
            _sphBesselk.conservativeResize(_maxOrder+2);
            _riccatiBesselRatio.conservativeResize(_maxOrder+2);

        }
    }

    array1DComplexType getSphBesselj()
    {
        return _sphBesselj.segment(1,_maxOrder);
    }

    array1DComplexType getSphBessely()
    {
        return _sphBessely.segment(1,_maxOrder);
    }

    array1DComplexType getSphBesselk()
    {
        return _sphBesselk.segment(1,_maxOrder);
    }

    array1DComplexType getRiccatiBesselRatio()
    {
        return _riccatiBesselRatio.segment(1,_maxOrder);
    }

    array1DComplexType computeRiccatiBesselPsi(intType startFrom)
    {
        BOOST_ASSERT_MSG(startFrom == 1 || !startFrom,
            "The computeRiccatiBesselPsi function only takes 0 and 1");
        return _argument * _sphBesselj.segment( startFrom , _maxOrder );
    }

    array1DComplexType computeRiccatiBesselChi()
    {
        return _argument * _sphBessely.segment(1,_maxOrder);
    }

    array1DComplexType computeRiccatiBesselXi(intType startFrom)
    {
        BOOST_ASSERT_MSG(startFrom == 1 || !startFrom,
            "The computeRiccatiBesselXi function only takes 0 and 1");
        return _argument * _sphBesselk.segment( startFrom , _maxOrder );
    }

    array1DComplexType computeDRiccatiBesselPsi()
    {
        return half * ( _argument * ( _sphBesselj.segment(0,_maxOrder) -
                                      _sphBesselj.segment(2,_maxOrder) ) +
                                    _sphBesselj.segment(1,_maxOrder) );
    }

    array1DComplexType computeDRiccatiBesselXi()
    {

        return half * ( _argument * ( _sphBesselk.segment(0,_maxOrder) -
                                      _sphBesselk.segment(2,_maxOrder) ) +
                                    _sphBesselk.segment(1,_maxOrder) );
    }
};
#endif //BESSEL_FUNCTIONS_HPP
