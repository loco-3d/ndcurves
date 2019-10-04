// Copyright (c) 2015-2018, CNRS
// Authors: Justin Carpentier <jcarpent@laas.fr>

#ifndef __curves_serialization_fwd_hpp__
#define __curves_api_serialization_fwd_hpp__

#include <boost/serialization/nvp.hpp>

#define BOOST_SERIALIZATION_MAKE_NVP(member) boost::serialization::make_nvp(##member, member)

#endif  // ifndef __multicontact_api_serialization_fwd_hpp__
