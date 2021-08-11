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
 * @brief This is a simplified version of the Particle class and is 
 * taken from MieLAM software package. The scatterer class implements 
 * robust and stable algorithms for calculation 
 * of scattering coefficients of internal and scattered fields 
 * (arbitrary precision capablilites and modified versions of an...dn
 * coefficients are removed). 
 *
 * @author Sina Haeri (s.haeri@ed.ac.uk)
 * @author Soroosh Haeri (soroosh.haeri@gmail.com)
 */

#ifndef SCATTERER_HPP
#define SCATTERER_HPP

#include <boost/multiprecision/eigen.hpp>
#include <boost/math/constants/constants.hpp>
#include "boostAssertMessage.hpp"
#include "besselFunctions.hpp"

#include <iostream>
#include <exception>

class scatterer
{
private:

    const double one = double(1.0);
    const double zero = double(0.0);
    const double two = double(2.0);
    const complexd c_i = {zero, one};
    const complexd c_one = {one, zero};

    typedef Eigen::Array<complexd, Eigen::Dynamic, 1> array1DComplexType;
    typedef Eigen::Array<double, Eigen::Dynamic, 1> array1DFloatType;

    double _radius;
    complexd _relRefractiveIndex;
    double _beamWaveLength;
    double _envPermeability;
    double _particlePermeability;
    int _maxOrder;
    complexd _alpha;
    besselFunctions *_besselObjAlpha = NULL;
    complexd _beta;
    besselFunctions *_besselObjBeta = NULL;
    array1DComplexType _linArray;

    array1DComplexType _anDirect()
    {
        array1DComplexType dPsiBeta = _particlePermeability *
                                      _besselObjBeta->computeDRiccatiBesselPsi();
        complexd tmpComp = _envPermeability * _relRefractiveIndex;
        array1DComplexType psiBeta = tmpComp *
                                     _besselObjBeta->computeRiccatiBesselPsi(1);

        return ( dPsiBeta * _besselObjAlpha->computeRiccatiBesselPsi(1) -
                 psiBeta * _besselObjAlpha->computeDRiccatiBesselPsi()
               ) /
               (
                 dPsiBeta * _besselObjAlpha->computeRiccatiBesselXi(1) -
                 psiBeta * _besselObjAlpha->computeDRiccatiBesselXi()
               );
    }

    array1DComplexType _bnDirect()
    {
        complexd tmpComp = _envPermeability * _relRefractiveIndex;
        array1DComplexType dPsiBeta = tmpComp *
                                      _besselObjBeta->computeDRiccatiBesselPsi();
        array1DComplexType psiBeta = _particlePermeability *
                                     _besselObjBeta->computeRiccatiBesselPsi(1);

        return ( dPsiBeta * _besselObjAlpha->computeRiccatiBesselPsi(1) -
                 psiBeta * _besselObjAlpha->computeDRiccatiBesselPsi()
               ) /
               (
                 dPsiBeta * _besselObjAlpha->computeRiccatiBesselXi(1) -
                 psiBeta * _besselObjAlpha->computeDRiccatiBesselXi()
               );
    }

    array1DComplexType _cnDirect()
    {
        array1DComplexType xiAlpha = _besselObjAlpha->computeRiccatiBesselXi(1);
        array1DComplexType dXiAlpha= _besselObjAlpha->computeDRiccatiBesselXi();

        return ( _particlePermeability * _relRefractiveIndex *
                 ( xiAlpha * _besselObjAlpha->computeDRiccatiBesselPsi() -
                   dXiAlpha * _besselObjAlpha->computeRiccatiBesselPsi(1)  )
               ) /
               (
                 _particlePermeability * xiAlpha *
                 _besselObjBeta->computeDRiccatiBesselPsi() -
                 (_envPermeability * _relRefractiveIndex) * dXiAlpha *
                 _besselObjBeta->computeRiccatiBesselPsi(1)
               );
    }

    array1DComplexType _dnDirect()
    {
        array1DComplexType xiAlpha = _besselObjAlpha->computeRiccatiBesselXi(1);
        array1DComplexType dXiAlpha= _besselObjAlpha->computeDRiccatiBesselXi();

        return ( _particlePermeability * _relRefractiveIndex * _relRefractiveIndex *
                 ( xiAlpha * _besselObjAlpha->computeDRiccatiBesselPsi() -
                   dXiAlpha * _besselObjAlpha->computeRiccatiBesselPsi(1)  )
               ) /
               (
                 ( _envPermeability * _relRefractiveIndex ) * xiAlpha *
                 _besselObjBeta->computeDRiccatiBesselPsi() -
                 _particlePermeability * dXiAlpha *
                 _besselObjBeta->computeRiccatiBesselPsi(1)
               );
    }

public:
    scatterer():
        _radius(one),
        _relRefractiveIndex(one,zero),
        _beamWaveLength(one),
        _envPermeability(one),
        _particlePermeability(one),
        _maxOrder(1)
    {
        using boost::math::constants::pi;

        _alpha = two * pi<double>() * _radius / _beamWaveLength;
        _besselObjAlpha = new besselFunctions( _maxOrder, _alpha );

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta = new besselFunctions( _maxOrder, _beta );

        _linArray = array1DComplexType::LinSpaced( _maxOrder, 1, _maxOrder);

    }

    scatterer(int order, double radius, double lambda,
              double mup, double mue, double envRefractiveInd,
              complexd partRefractiveInd):
        _radius(radius),
        _relRefractiveIndex(partRefractiveInd/envRefractiveInd),
        _beamWaveLength(lambda/envRefractiveInd),
        _envPermeability(mup),
        _particlePermeability(mue)
    {
        using boost::math::constants::pi;

        BOOST_ASSERT_MSG(order >= 1, "Class Scatterer order must be greater than 0.");

        _maxOrder = order;

        _alpha = two * pi<double>() * _radius / _beamWaveLength;
        _besselObjAlpha = new besselFunctions( _maxOrder, _alpha );

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta = new besselFunctions( _maxOrder, _beta );

        _linArray = array1DComplexType::LinSpaced( _maxOrder, 1, _maxOrder);

    }

    scatterer(const scatterer& other) :
    _radius(other._radius),
    _relRefractiveIndex(other._relRefractiveIndex),
    _beamWaveLength(other._beamWaveLength),
    _envPermeability(other._envPermeability),
    _particlePermeability(other._particlePermeability),
    _maxOrder(other._maxOrder),
    _alpha(other._alpha),
    _beta(other._beta),
    _linArray(other._linArray)
    {
        _besselObjAlpha = new besselFunctions(*(other._besselObjAlpha));
        _besselObjBeta = new besselFunctions(*(other._besselObjBeta));
    }

    scatterer& operator=(const scatterer& other)
    {
        if (this != &other)
        {
            _radius = other._radius;
            _relRefractiveIndex = other._relRefractiveIndex;
            _beamWaveLength = other._beamWaveLength;
            _envPermeability = other._envPermeability;
            _particlePermeability = other._particlePermeability;
            _maxOrder = other._maxOrder;
            _alpha = other._alpha;
            _beta = other._beta;
            _linArray = other._linArray;
            if(_besselObjAlpha != NULL)
            {
                delete _besselObjAlpha;
                _besselObjAlpha = NULL;
            }

            if(_besselObjBeta != NULL)
            {
                delete _besselObjBeta;
                _besselObjBeta = NULL;
            }

            _besselObjAlpha = new besselFunctions(*(other._besselObjAlpha));
            _besselObjBeta = new besselFunctions(*(other._besselObjBeta));
        }

        return *this;
    }

    ~scatterer()
    {
        delete _besselObjBeta;
        delete _besselObjAlpha;
    }

    void setRadius(double rad)
    {
        using boost::math::constants::pi;
        _radius = rad;

        _alpha = two * pi<double>() * _radius / _beamWaveLength;
        _besselObjAlpha->setArgument(_alpha);

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);

    }

    void setMaxOrder(int order)
    {
        _maxOrder = order;
        _besselObjAlpha->setMaxOrder(_maxOrder);
        _besselObjBeta->setMaxOrder(_maxOrder);
        _linArray = array1DComplexType::LinSpaced( _maxOrder, 1, _maxOrder);
    }

    int getMaxOrder()
    { return _maxOrder; }

    void setWaveLength(double lambda)
    {
        using boost::math::constants::pi;
        _beamWaveLength = lambda;

        _alpha = two * pi<double>() * _radius / _beamWaveLength;
        _besselObjAlpha->setArgument(_alpha);

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);

    }

    void setRefractiveInd(complexd relRefInd)
    {
        using boost::math::constants::pi;
        _relRefractiveIndex = relRefInd;

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);
    }

    array1DComplexType computeAn()
    {
        return _anDirect();
    }

    array1DComplexType computeBn()
    {
        return _bnDirect();
    }

    array1DComplexType computeCn()
    {
        return _cnDirect();
    }

    array1DComplexType computeDn()
    {

        return _dnDirect();

    }

};

#endif //SCATTERER_HPP
