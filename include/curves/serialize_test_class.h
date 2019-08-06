/**
* \file polynomial.h
* \brief Definition of a cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the polynomial struct.
* It allows the creation and evaluation of natural
* smooth splines of arbitrary dimension and order
*/


#ifndef _STRUCT_TEST_SER
#define _STRUCT_TEST_SER

#include "MathDefs.h"

#include "curve_abc.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <fstream>
#include <string>
#include <stdexcept>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"

namespace curves
{
struct serialize_test_class : public serialization::Serializable< serialize_test_class >
{
    public:

    serialize_test_class(int a) :a_(a) {}

    int a_;

    void serialize_to_file(std::string path)
    {
        std::ofstream ofile(path.c_str());
        boost::archive::text_oarchive oTextArchive(ofile);
        oTextArchive << (*this);    // sérialisation de d
    }

    void deserialize_from_file(std::string path)
    {
        std::ifstream ifile(path.c_str());
        boost::archive::text_iarchive iTextArchive(ifile);
        iTextArchive >> (*this);     // désérialisation dans d
    }

    // Serialization of the class
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(a_);
    }

}; //class polynomial
} // namespace curves
#endif //_STRUCT_POLYNOMIAL

