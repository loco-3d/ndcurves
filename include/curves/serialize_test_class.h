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

template<typename Time, std::size_t Dim=3>
struct TestStructParent {
    TestStructParent(){}
    virtual ~TestStructParent(){}
};

template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
         typename Point= Eigen::Matrix<Numeric, Eigen::Dynamic, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct TestStruct : public TestStructParent<Time, Dim> {//: public curve_abc<Time, Numeric, Safe, Point>{
public:
    TestStruct(){}
    TestStruct(int i, std::string c, float f1, float f2, Eigen::Matrix<double, Eigen::Dynamic, 1> mat){
        a=i;
        s=c;
        v.push_back(f1);
        v.push_back(f2);
        matrix_test = mat;
    }

    ~TestStruct()
    {
        // NOTHING
    }

    Point operator()(const time_t t) const
    {
        return Point(0);
    }

    ///  \brief Evaluation of the derivative of order N of spline at time t.
    ///  \param t : the time when to evaluate the spline.
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derivative spline at time t.
    Point derivate(const time_t t, const std::size_t order) const
    {
        return Point(0);
    }

    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$ lower bound of time range.
    Numeric min() const {return 0;}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$ upper bound of time range.
    Numeric max() const {return 0;}

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
         ar & a & s & v & matrix_test;
    }
    int a;
    std::string s;
    std::vector<float> v;
    Eigen::Matrix<double, Eigen::Dynamic, 1> matrix_test;
};

struct serialize_test_class
{
    public:
    typedef Eigen::Matrix<double,1,1> point_one;
    typedef TestStruct<double, double, 1, true, point_one> test_t;
    test_t t1, t2;

    serialize_test_class()
    {
        Eigen::Matrix<double, 3, 1> mat1(1,2,3);
        Eigen::Matrix<double, 3, 1> mat2(4,4,4);
        t1 = test_t(2, std::string("essai 1"), 2.0, 3.5, mat1);
        t2 = test_t(0, std::string("essai 2"), 0.0, 0.5, mat2);
    }

    void serialize_to_file(std::string path)
    {
        std::cout<<"serialize"<<std::endl;
        std::ofstream ofile(path.c_str());
        boost::archive::text_oarchive oTextArchive(ofile);
        oTextArchive << t1;    // serialization
    }

    void deserialize_from_file(std::string path)
    {
        std::cout<<"deserialize"<<std::endl;
        std::ifstream ifile(path.c_str());
        boost::archive::text_iarchive iTextArchive(ifile);
        iTextArchive >> t2;     // deserialization
    }

}; //class polynomial
} // namespace curves
#endif //_STRUCT_POLYNOMIAL

