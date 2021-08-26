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

    array1DComplexType _An, _Bn, _Cn, _Dn;

    void _computeBesselsAlpha()
    {
        _besselObjAlpha->computeDRiccatiBesselPsi();
        _besselObjAlpha->computeRiccatiBesselPsi(1);
        _besselObjAlpha->computeDRiccatiBesselXi();
        _besselObjAlpha->computeRiccatiBesselXi(1);
    }

    void _computeBesselsBeta()
    {
        _besselObjBeta->computeDRiccatiBesselPsi();
        _besselObjBeta->computeRiccatiBesselPsi(1);
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
        _computeBesselsAlpha();

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta = new besselFunctions( _maxOrder, _beta );
        _computeBesselsBeta();
        
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
        _computeBesselsAlpha();

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta = new besselFunctions( _maxOrder, _beta );
        _computeBesselsBeta();

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
        _computeBesselsAlpha();

        _besselObjBeta = new besselFunctions(*(other._besselObjBeta));
        _computeBesselsBeta();     
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
            _computeBesselsAlpha();

            _besselObjBeta = new besselFunctions(*(other._besselObjBeta));
            _computeBesselsBeta();

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
        _computeBesselsAlpha();
        
        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);
        _computeBesselsBeta(); 
    }

    void setMaxOrder(int order)
    {
        _maxOrder = order;
        _besselObjAlpha->setMaxOrder(_maxOrder);
        _computeBesselsAlpha();

        _besselObjBeta->setMaxOrder(_maxOrder);
        _computeBesselsBeta(); 

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
        _computeBesselsAlpha();

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);
        _computeBesselsBeta();

    }

    void setRefractiveInd(complexd relRefInd)
    {
        using boost::math::constants::pi;
        _relRefractiveIndex = relRefInd;

        _beta = _alpha * _relRefractiveIndex;
        _besselObjBeta->setArgument(_beta);
        _computeBesselsBeta();
    }

    complexd getRelRefInd() { return _relRefractiveIndex;}
    complexd getSizeParam() {return _alpha;}
    complexd getBeta() {return _beta;}

    void computeAn()
    {
        array1DComplexType dPsiBeta = _particlePermeability *
                                      _besselObjBeta->getDRiccatiBesselPsi();
        complexd tmpComp = _envPermeability * _relRefractiveIndex;
        array1DComplexType psiBeta = tmpComp *
                                     _besselObjBeta->getRiccatiBesselPsi();

        //std::cout << "_alpha: " << _alpha << std::endl; 
        _An =  ( dPsiBeta * _besselObjAlpha->getRiccatiBesselPsi() -
                 psiBeta * _besselObjAlpha->getDRiccatiBesselPsi()
               ) /
               (
                 dPsiBeta * _besselObjAlpha->getRiccatiBesselXi() -
                 psiBeta * _besselObjAlpha->getDRiccatiBesselXi()
               );
        /*std::cout << "DRiccatiBesselPsi(beta)" << std::endl;
        std::cout << _besselObjBeta->getDRiccatiBesselPsi() << std::endl << std::endl;

        std::cout << std::endl;
        std::cout << "_particlePermeability" << _particlePermeability << std::endl;

        std::cout << std::endl;
        std::cout << "_particlePermeability * _besselObjBeta->getDRiccatiBesselPsi(): " << 
        dPsiBeta.size() << ";  " << dPsiBeta << std::endl;

        std::cout << std::endl;
        std::cout << "RiccatiBesselPsi(beta)" << std::endl;
        std::cout << _besselObjBeta->getRiccatiBesselPsi() << std::endl << std::endl;
        
        std::cout << std::endl;    
        std::cout << "_envPermeability * _relRefractiveIndex" << tmpComp << std::endl;

        std::cout << std::endl;
        std::cout << "getRiccatiBesselPsi(alpha)" << std::endl;
        std::cout << _besselObjAlpha->getRiccatiBesselPsi() << std::endl << std::endl;

        std::cout << std::endl;
        std::cout << "getDRiccatiBesselPsi(alpha)" << std::endl;
        std::cout << _besselObjAlpha->getDRiccatiBesselPsi() << std::endl << std::endl;

        std::cout << std::endl;
        std::cout << "getRiccatiBesselXi(alpha)" << std::endl;
        std::cout << _besselObjAlpha->getRiccatiBesselXi() << std::endl << std::endl;

        std::cout << std::endl;
        std::cout << "getDRiccatiBesselPsi(alpha)" << std::endl;
        std::cout << _besselObjAlpha->getDRiccatiBesselXi() << std::endl << std::endl;

        std::cout << "An from within the class:" << std::endl;
        std::cout << _An << std::endl << std::endl;*/
    }

    void computeBn()
    {
        complexd tmpComp = _envPermeability * _relRefractiveIndex;
        array1DComplexType dPsiBeta = tmpComp *
                                      _besselObjBeta->getDRiccatiBesselPsi();
        array1DComplexType psiBeta = _particlePermeability *
                                     _besselObjBeta->getRiccatiBesselPsi();

        _Bn =  ( dPsiBeta * _besselObjAlpha->getRiccatiBesselPsi() -
                 psiBeta * _besselObjAlpha->getDRiccatiBesselPsi()
               ) /
               (
                 dPsiBeta * _besselObjAlpha->getRiccatiBesselXi() -
                 psiBeta * _besselObjAlpha->getDRiccatiBesselXi()
               );
    }

    void computeCn()
    {
        array1DComplexType xiAlpha = _besselObjAlpha->getRiccatiBesselXi();
        array1DComplexType dXiAlpha= _besselObjAlpha->getDRiccatiBesselXi();

        _Cn = ( _particlePermeability * _relRefractiveIndex *
                 ( xiAlpha * _besselObjAlpha->getDRiccatiBesselPsi() -
                   dXiAlpha * _besselObjAlpha->getRiccatiBesselPsi()  )
               ) /
               (
                 _particlePermeability * xiAlpha *
                 _besselObjBeta->getDRiccatiBesselPsi() -
                 (_envPermeability * _relRefractiveIndex) * dXiAlpha *
                 _besselObjBeta->getRiccatiBesselPsi()
               );
    }

    void computeDn()
    {
        array1DComplexType xiAlpha = _besselObjAlpha->getRiccatiBesselXi();
        array1DComplexType dXiAlpha= _besselObjAlpha->getDRiccatiBesselXi();

        _Dn = ( _particlePermeability * _relRefractiveIndex * _relRefractiveIndex *
                 ( xiAlpha * _besselObjAlpha->getDRiccatiBesselPsi() -
                   dXiAlpha * _besselObjAlpha->getRiccatiBesselPsi()  )
               ) /
               (
                 ( _envPermeability * _relRefractiveIndex ) * xiAlpha *
                 _besselObjBeta->getDRiccatiBesselPsi() -
                 _particlePermeability * dXiAlpha *
                 _besselObjBeta->getRiccatiBesselPsi()
               );
    }


    const array1DComplexType& getAn() {return _An;}
    const array1DComplexType& getBn() {return _Bn;}
    const array1DComplexType& getCn() {return _Cn;}
    const array1DComplexType& getDn() {return _Dn;}

    void printBetaRicBessel()
    {
        std::cout << "DRiccatiBesselPsi(beta)" << std::endl << std::endl; 
        std::cout << _besselObjBeta->getDRiccatiBesselPsi() << std::endl;

        std::cout << "RiccatiBesselPsi(beta)" << std::endl << std::endl; 
        std::cout << _besselObjBeta->getRiccatiBesselPsi() << std::endl;

        std::cout << "RiccatiBesselChi(beta)" << std::endl << std::endl;
        _besselObjBeta->computeRiccatiBesselChi(); 
        std::cout << _besselObjBeta->getRiccatiBesselChi() << std::endl;

    }

    void printAlphaRicBessel()
    {
        std::cout << "DRiccatiBesselXi(alpha)" << std::endl; 
        std::cout << _besselObjAlpha->getDRiccatiBesselXi() << std::endl;

        std::cout << "RiccatiBesselXi(alpha)" << std::endl; 
        std::cout << _besselObjAlpha->getRiccatiBesselXi() << std::endl;

    }

};

#endif //SCATTERER_HPP
