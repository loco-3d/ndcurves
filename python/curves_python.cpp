#include "python_variables.h"
#include "archive_python_binding.h"
#include "optimization_python.h"

#include <vector>



namespace curves
{
  using namespace boost::python;

  /* base wrap of curve_abc : must implement all pure virtual method of curve_abc */
  struct CurveWrapper : curve_abc_t, wrapper<curve_abc_t>
  {
      point_t operator()(const real) { return this->get_override("operator()")();}
      point_t derivate(const real, const std::size_t) { return this->get_override("derivate")();}
      std::size_t dim() { return this->get_override("dim")();}
      real min() { return this->get_override("min")();}
      real max() { return this->get_override("max")();}

  };
  /* end base wrap of curve_abc */

  /* Template constructor bezier */
  template <typename Bezier, typename PointList, typename T_Point>
  Bezier* wrapBezierConstructorTemplate(const PointList& array, const real T_min =0., const real T_max =1.)
  {
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), T_min, T_max);
  }

  template <typename Bezier, typename PointList, typename T_Point, typename CurveConstraints>
  Bezier* wrapBezierConstructorConstraintsTemplate(const PointList& array, const CurveConstraints& constraints,
                                                   const real T_min =0., const real T_max =1.)
  {
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), constraints, T_min, T_max);
  }
  /* End Template constructor bezier */
  /* Helper converter constraintsX -> constraints 3 */
  curve_constraints3_t convertToConstraints3(curve_constraints_t constraintsX){
    curve_constraints3_t constraints3;
    constraints3.init_vel =point3_t(constraintsX.init_vel);
    constraints3.init_acc = point3_t(constraintsX.init_acc);
    constraints3.end_vel = point3_t(constraintsX.end_vel);
    constraints3.end_acc = point3_t(constraintsX.end_acc);
    return constraints3;
  }
  /* END helper converter constraintsX -> constraints 3 */

  /*3D constructors bezier */
  bezier3_t* wrapBezier3Constructor(const pointX_list_t& array)
  {
    return wrapBezierConstructorTemplate<bezier3_t, pointX_list_t, t_point3_t>(array) ;
  }
  bezier3_t* wrapBezier3ConstructorBounds(const pointX_list_t& array, const real T_min, const real T_max)
  {
    return wrapBezierConstructorTemplate<bezier3_t, pointX_list_t, t_point3_t>(array, T_min, T_max) ;
  }
  bezier3_t* wrapBezier3ConstructorConstraints(const pointX_list_t& array, const curve_constraints_t& constraints)
  {
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, pointX_list_t, t_point3_t, curve_constraints3_t>(array, convertToConstraints3(constraints)) ;
  }
  bezier3_t* wrapBezier3ConstructorBoundsConstraints(const pointX_list_t& array, const curve_constraints_t& constraints,
                                                   const real T_min, const real T_max)
  {
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, pointX_list_t, t_point3_t, curve_constraints3_t>(array, convertToConstraints3(constraints),T_min, T_max) ;
  }
  /*END 3D constructors bezier */

  /*constructors bezier */
  bezier_t* wrapBezierConstructor(const pointX_list_t& array)
  {
    return wrapBezierConstructorTemplate<bezier_t, pointX_list_t, t_pointX_t>(array) ;
  }
  bezier_t* wrapBezierConstructorBounds(const pointX_list_t& array, const real T_min, const real T_max)
  {
    return wrapBezierConstructorTemplate<bezier_t, pointX_list_t, t_pointX_t>(array, T_min, T_max) ;
  }
  bezier_t* wrapBezierConstructorConstraints(const pointX_list_t& array, const curve_constraints_t& constraints)
  {
    return wrapBezierConstructorConstraintsTemplate<bezier_t, pointX_list_t, t_pointX_t, curve_constraints_t>(array, constraints) ;
  }
  bezier_t* wrapBezierConstructorBoundsConstraints(const pointX_list_t& array, const curve_constraints_t& constraints,
                                                   const real T_min, const real T_max)
  {
    return wrapBezierConstructorConstraintsTemplate<bezier_t, pointX_list_t, t_pointX_t, curve_constraints_t>(array, constraints,
                                                                                                            T_min, T_max) ;
  }
  /*END constructors bezier */

  /* Wrap Cubic hermite spline */
  t_pair_pointX_tangent_t getPairsPointTangent(const pointX_list_t& points, const pointX_list_t& tangents)
  {
    t_pair_pointX_tangent_t res;
    if (points.size() != tangents.size())
    {
      throw std::length_error("size of points and tangents must be the same");
    }
    for(int i =0;i<points.cols();++i)
    {
      res.push_back(pair_pointX_tangent_t(points.col(i), tangents.col(i)));
    }
    return res;
  }

  cubic_hermite_spline_t* wrapCubicHermiteSplineConstructor(const pointX_list_t& points, const pointX_list_t& tangents,
                                                            const time_waypoints_t& time_pts)
  {
    t_pair_pointX_tangent_t ppt = getPairsPointTangent(points, tangents);
    std::vector<real> time_control_pts;
    for( int i =0; i<time_pts.size(); ++i)
    {
      time_control_pts.push_back(time_pts[i]);
    }
    return new cubic_hermite_spline_t(ppt.begin(), ppt.end(), time_control_pts);
  }
  /* End wrap Cubic hermite spline */

  /* Wrap polynomial */
  polynomial_t* wrapPolynomialConstructor1(const coeff_t& array, const real min, const real max)
  {
    return new polynomial_t(array, min, max);
  }
  polynomial_t* wrapPolynomialConstructor2(const coeff_t& array)
  {
    return new polynomial_t(array, 0., 1.);
  }
  polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree1(const pointX_t& init,const pointX_t& end,const real min, const real max)
  {
    return new polynomial_t(init,end,min,max);
  }
  polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree3(const pointX_t& init,const pointX_t& d_init,const pointX_t& end,const pointX_t& d_end,const real min, const real max)
  {
    return new polynomial_t(init,d_init,end,d_end,min,max);
  }
  polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree5(const pointX_t& init,const pointX_t& d_init,const pointX_t& dd_init,const pointX_t& end,const point_t& d_end,const point_t& dd_end,const real min, const real max)
  {
    return new polynomial_t(init,d_init,dd_init,end,d_end,dd_end,min,max);
  }
  /* End wrap polynomial */

  /* Wrap piecewise curve */
  piecewise_polynomial_curve_t* wrapPiecewisePolynomialCurveConstructor(const polynomial_t& pol)
  {
    return new piecewise_polynomial_curve_t(pol);
  }
  piecewise_polynomial_curve_t* wrapPiecewisePolynomialCurveEmptyConstructor()
  {
    return new piecewise_polynomial_curve_t();
  }
  piecewise_bezier_curve_t* wrapPiecewiseBezierCurveConstructor(const bezier_t& bc)
  {
    return new piecewise_bezier_curve_t(bc);
  }
  piecewise_bezier_curve_t* wrapPiecewiseBezierCurveEmptyConstructor()
  {
    return new piecewise_bezier_curve_t();
  }
  piecewise_cubic_hermite_curve_t* wrapPiecewiseCubicHermiteCurveConstructor(const cubic_hermite_spline_t& ch)
  {
    return new piecewise_cubic_hermite_curve_t(ch);
  }
  piecewise_cubic_hermite_curve_t* wrapPiecewiseCubicHermiteCurveEmptyConstructor()
  {
    return new piecewise_cubic_hermite_curve_t();
  }
  static piecewise_polynomial_curve_t discretPointToPolynomialC0(const pointX_list_t& points, const time_waypoints_t& time_points){
    t_pointX_t points_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points);
    t_time_t time_points_list = vectorFromEigenVector<time_waypoints_t,t_time_t>(time_points);
    return piecewise_polynomial_curve_t::convert_discrete_points_to_polynomial<polynomial_t>(points_list,time_points_list);
  }
  static piecewise_polynomial_curve_t discretPointToPolynomialC1(const pointX_list_t& points,const pointX_list_t& points_derivative, const time_waypoints_t& time_points){
    t_pointX_t points_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points);
    t_pointX_t points_derivative_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points_derivative);
    t_time_t time_points_list = vectorFromEigenVector<time_waypoints_t,t_time_t>(time_points);
    return piecewise_polynomial_curve_t::convert_discrete_points_to_polynomial<polynomial_t>(points_list,points_derivative_list,time_points_list);
  }
  static piecewise_polynomial_curve_t discretPointToPolynomialC2(const pointX_list_t& points,const pointX_list_t& points_derivative,const pointX_list_t& points_second_derivative, const time_waypoints_t& time_points){
    t_pointX_t points_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points);
    t_pointX_t points_derivative_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points_derivative);
    t_pointX_t points_second_derivative_list = vectorFromEigenArray<pointX_list_t,t_pointX_t>(points_second_derivative);

    t_time_t time_points_list = vectorFromEigenVector<time_waypoints_t,t_time_t>(time_points);
    return piecewise_polynomial_curve_t::convert_discrete_points_to_polynomial<polynomial_t>(points_list,points_derivative_list,points_second_derivative_list,time_points_list);
  }
  void addFinalPointC0(piecewise_polynomial_curve_t self,const pointX_t& end,const real time){
    if(self.is_continuous(1))
      std::cout<<"Warning: by adding this final point to the piecewise curve, you loose C1 continuity and only guarantee C0 continuity."<<std::endl;
    polynomial_t pol(self(self.max()),end,self.max(),time);
    self.add_curve(pol);
  }
  void addFinalPointC1(piecewise_polynomial_curve_t self,const pointX_t& end,const pointX_t& d_end,const real time){
    if(self.is_continuous(2))
      std::cout<<"Warning: by adding this final point to the piecewise curve, you loose C2 continuity and only guarantee C1 continuity."<<std::endl;
    if(!self.is_continuous(1))
      std::cout<<"Warning: the current piecewise curve is not C1 continuous."<<std::endl;
    polynomial_t pol(self(self.max()),self.derivate(self.max(),1),end,d_end,self.max(),time);
    self.add_curve(pol);
  }
  void addFinalPointC2(piecewise_polynomial_curve_t self,const pointX_t& end,const pointX_t& d_end,const pointX_t& dd_end,const real time){
    if(self.is_continuous(3))
      std::cout<<"Warning: by adding this final point to the piecewise curve, you loose C3 continuity and only guarantee C2 continuity."<<std::endl;
    if(!self.is_continuous(2))
      std::cout<<"Warning: the current piecewise curve is not C2 continuous."<<std::endl;
    polynomial_t pol(self(self.max()),self.derivate(self.max(),1),self.derivate(self.max(),2),end,d_end,dd_end,self.max(),time);
    self.add_curve(pol);
  }


  /* end wrap piecewise polynomial curve */

  /* Wrap exact cubic spline */
  t_waypoint_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp)
  {
    t_waypoint_t res;
    for(int i =0;i<array.cols();++i)
    {
      res.push_back(std::make_pair(time_wp(i), array.col(i)));
    }
    return res;
  }

  exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array, const time_waypoints_t& time_wp)
  {
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end());
  }

  exact_cubic_t* wrapExactCubicConstructorConstraint(const coeff_t& array, const time_waypoints_t& time_wp,
                                                     const curve_constraints_t& constraints)
  {
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end(), constraints);
  }

  /// For constraints XD
  point_t get_init_vel(const curve_constraints_t& c)
  {
    return c.init_vel;
  }

  point_t get_init_acc(const curve_constraints_t& c)
  {
    return c.init_acc;
  }

  point_t get_init_jerk(const curve_constraints_t& c)
  {
      return c.init_jerk;
  }

  point_t get_end_vel(const curve_constraints_t& c)
  {
    return c.end_vel;
  }

  point_t get_end_acc(const curve_constraints_t& c)
  {
    return c.end_acc;
  }

  point_t get_end_jerk(const curve_constraints_t& c)
  {
      return c.end_jerk;
  }

  void set_init_vel(curve_constraints_t& c, const point_t& val)
  {
    c.init_vel = val;
  }

  void set_init_acc(curve_constraints_t& c, const point_t& val)
  {
    c.init_acc = val;
  }

  void set_init_jerk(curve_constraints_t& c, const point_t& val)
  {
      c.init_jerk = val;
  }

  void set_end_vel(curve_constraints_t& c, const point_t& val)
  {
    c.end_vel = val;
  }

  void set_end_acc(curve_constraints_t& c, const point_t& val)
  {
    c.end_acc = val;
  }

  void set_end_jerk(curve_constraints_t& c, const point_t& val)
  {
      c.end_jerk = val;
  }

  matrix_vector* bezier_linear_variable_t_operator_call(const bezier_linear_variable_t* b, const double t)
  {
      bezier_linear_variable_t::point_t p = b->operator ()(t);
      matrix_vector* res = new matrix_vector(p.B(), p.c());
      return res;
  }

  bezier_t* bezier_linear_variable_t_evaluate(const bezier_linear_variable_t* b, const pointX_t& x)
  {
     return new bezier_t(evaluateLinear<bezier_t, bezier_linear_variable_t>(*b, x));
  }


  /* End wrap exact cubic spline */


  // TO DO : Replace all load and save function for serialization in class by using 
  //         SerializableVisitor in archive_python_binding.
  BOOST_PYTHON_MODULE(curves)
  {
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();
    eigenpy::enableEigenPySpecific<pointX_t,pointX_t>();
    eigenpy::enableEigenPySpecific<pointX_list_t,pointX_list_t>();
    eigenpy::enableEigenPySpecific<coeff_t,coeff_t>();
    eigenpy::enableEigenPySpecific<point_list_t,point_list_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/
    /** END eigenpy init**/
    class_<CurveWrapper,boost::noncopyable>("curve",no_init)
        .def("__call__", pure_virtual(&curve_abc_t::operator()),"Evaluate the curve at the given time.",args("self","t"))
        .def("derivate", pure_virtual(&curve_abc_t::derivate),"Evaluate the derivative of order N of curve at time t.",args("self","t","N"))
        .def("min", pure_virtual(&curve_abc_t::min), "Get the LOWER bound on interval definition of the curve.")
        .def("max", pure_virtual(&curve_abc_t::max),"Get the HIGHER bound on interval definition of the curve.")
        .def("dim", pure_virtual(&curve_abc_t::dim),"Get the dimension of the curve.")
        .def("saveAsText", pure_virtual(&curve_abc_t::saveAsText<curve_abc_t>),bp::args("filename"),"Saves *this inside a text file.")
        .def("loadFromText",pure_virtual(&curve_abc_t::loadFromText<curve_abc_t>),bp::args("filename"),"Loads *this from a text file.")
        .def("saveAsXML",pure_virtual(&curve_abc_t::saveAsXML<curve_abc_t>),bp::args("filename","tag_name"),"Saves *this inside a XML file.")
        .def("loadFromXML",pure_virtual(&curve_abc_t::loadFromXML<curve_abc_t>),bp::args("filename","tag_name"),"Loads *this from a XML file.")
        .def("saveAsBinary",pure_virtual(&curve_abc_t::saveAsBinary<curve_abc_t>),bp::args("filename"),"Saves *this inside a binary file.")
        .def("loadFromBinary",pure_virtual(&curve_abc_t::loadFromBinary<curve_abc_t>),bp::args("filename"),"Loads *this from a binary file.")
        ;

    /** BEGIN bezier3 curve**/
    class_<bezier3_t, bases<curve_abc_t> >("bezier3", init<>())
      .def("__init__", make_constructor(&wrapBezier3Constructor))
      .def("__init__", make_constructor(&wrapBezier3ConstructorBounds))
      .def("__init__", make_constructor(&wrapBezier3ConstructorConstraints))
      .def("__init__", make_constructor(&wrapBezier3ConstructorBoundsConstraints))
      .def("compute_derivate", &bezier3_t::compute_derivate)
      .def("compute_primitive", &bezier3_t::compute_primitive)
      .def("waypointAtIndex", &bezier3_t::waypointAtIndex)
      .def_readonly("degree", &bezier3_t::degree_)
      .def_readonly("nbWaypoints", &bezier3_t::size_)
      .def("saveAsText", &bezier3_t::saveAsText<bezier3_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&bezier3_t::loadFromText<bezier3_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&bezier3_t::saveAsXML<bezier3_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&bezier3_t::loadFromXML<bezier3_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&bezier3_t::saveAsBinary<bezier3_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&bezier3_t::loadFromBinary<bezier3_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<bezier_t>())
    ;
    /** END bezier3 curve**/
    /** BEGIN bezier curve**/
    class_<bezier_t, bases<curve_abc_t> >("bezier", init<>())
      .def("__init__", make_constructor(&wrapBezierConstructor))
      .def("__init__", make_constructor(&wrapBezierConstructorBounds))
      .def("__init__", make_constructor(&wrapBezierConstructorConstraints))
      .def("__init__", make_constructor(&wrapBezierConstructorBoundsConstraints))
      .def("compute_derivate", &bezier_t::compute_derivate)
      .def("compute_primitive", &bezier_t::compute_primitive)
      .def("waypointAtIndex", &bezier_t::waypointAtIndex)
      .def_readonly("degree", &bezier_t::degree_)
      .def_readonly("nbWaypoints", &bezier_t::size_)
      .def("saveAsText", &bezier_t::saveAsText<bezier_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&bezier_t::loadFromText<bezier_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&bezier_t::saveAsXML<bezier_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&bezier_t::loadFromXML<bezier_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&bezier_t::saveAsBinary<bezier_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&bezier_t::loadFromBinary<bezier_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<bezier_t>())
    ;
    /** END bezier curve**/
    /** BEGIN variable points bezier curve**/
    class_<matrix_pair>
        ("matrix_pair", no_init)
        .def_readonly("A", &matrix_pair::A)
        .def_readonly("b", &matrix_pair::b)
        ;

    class_<matrix_vector>
        ("matrix_vector", no_init)
        .def_readonly("A", &matrix_vector::A)
        .def_readonly("b", &matrix_vector::b)
        ;

    class_<LinearBezierVector>
    ("bezierVarVector", no_init)
      .def_readonly("size", &LinearBezierVector::size)
      .def("at", &LinearBezierVector::at, return_value_policy<manage_new_object>())
    ;
    class_<bezier_linear_variable_t>
    ("bezierVar", no_init)
        .def("__init__", make_constructor(&wrapBezierLinearConstructor))
        .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
        .def("min", &bezier_linear_variable_t::min)
        .def("max", &bezier_linear_variable_t::max)
        //.def("__call__", &bezier_linear_control_t::operator())
        .def("__call__", &bezier_linear_variable_t_operator_call, bp::return_value_policy<bp::manage_new_object>())
        .def("evaluate", &bezier_linear_variable_t_evaluate, bp::return_value_policy<bp::manage_new_object>())
        .def("derivate", &bezier_linear_variable_t::derivate)
        .def("compute_derivate", &bezier_linear_variable_t::compute_derivate)
        .def("compute_primitive", &bezier_linear_variable_t::compute_primitive)
        .def("split", &split_py, return_value_policy<manage_new_object>())
        .def("waypoints", &wayPointsToLists, return_value_policy<manage_new_object>())
        .def_readonly("degree", &bezier_linear_variable_t::degree_)
        .def_readonly("nbWaypoints", &bezier_linear_variable_t::size_)
        ;


    class_<quadratic_variable_t >
        ("cost", no_init)
        .add_property("A", &cost_t_quad)
        .add_property("b", &cost_t_linear)
        .add_property("c", &cost_t_constant)
        ;

    /** END variable points bezier curve**/
    /** BEGIN polynomial curve function**/
    class_<polynomial_t , bases<curve_abc_t> >("polynomial",  init<>())
      .def("__init__", make_constructor(&wrapPolynomialConstructor1,default_call_policies(),args("coeffs","min","max")),
           "Create polynomial spline from an Eigen matrix of coefficient defined for t \in [min,max]."
           " The matrix should contain one coefficient per column, from the zero order coefficient,up to the highest order."
           " Spline order is given by the number of the columns -1.")
      .def("__init__", make_constructor(&wrapPolynomialConstructor2,default_call_policies(),arg("coeffs")),
           "Create polynomial spline from an Eigen matrix of coefficient defined for t \in [0,1]."
           " The matrix should contain one coefficient per column, from the zero order coefficient,up to the highest order."
           " Spline order is given by the number of the columns -1.")
      .def("__init__", make_constructor(&wrapPolynomialConstructorFromBoundaryConditionsDegree1,
                                        default_call_policies(),args("init","end","min","max")),
           "Create a polynomial of degree 1 defined for t \in [min,max], "
           "such that c(min) == init and c(max) == end.")
      .def("__init__", make_constructor(&wrapPolynomialConstructorFromBoundaryConditionsDegree3,
                                        default_call_policies(),args("init","d_init","end","d_end","min","max")),
          "Create a polynomial of degree 3 defined for t \in [min,max], "
          "such that c(min) == init and c(max) == end"
          " dc(min) == d_init and dc(max) == d_end")
      .def("__init__", make_constructor(&wrapPolynomialConstructorFromBoundaryConditionsDegree5,
                                        default_call_policies(),
                                        args("init","d_init","dd_init","end","d_end","dd_end","min","max")),
           "Create a polynomial of degree 5 defined for t \in [min,max], "
           "such that c(min) == init and c(max) == end"
           " dc(min) == d_init and dc(max) == d_end"
           " ddc(min) == dd_init and ddc(max) == dd_end")
      .def("coeffAtDegree", &polynomial_t::coeffAtDegree)
      .def("coeff", &polynomial_t::coeff)
      .def("compute_derivate", &polynomial_t::compute_derivate,"Compute derivative of order N of curve at time t.")
      .def("saveAsText", &polynomial_t::saveAsText<polynomial_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&polynomial_t::loadFromText<polynomial_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&polynomial_t::saveAsXML<polynomial_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&polynomial_t::loadFromXML<polynomial_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&polynomial_t::saveAsBinary<polynomial_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&polynomial_t::loadFromBinary<polynomial_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    /** END polynomial function**/
    /** BEGIN piecewise curve function **/
    class_<piecewise_polynomial_curve_t , bases<curve_abc_t> >
    ("piecewise_polynomial_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewisePolynomialCurveConstructor,default_call_policies(),arg("curve")),
           "Create a peicewise-polynomial curve containing the given polynomial curve.")
      .def("FromPointsList",&discretPointToPolynomialC0,
           "Create a piecewise-polynomial connecting exactly all the given points at the given time. The created piecewise is C0 continuous.",args("points","time_points"))
       .def("FromPointsList",&discretPointToPolynomialC1,
             "Create a piecewise-polynomial connecting exactly all the given points at the given time and respect the given points derivative values. The created piecewise is C1 continuous.",args("points","points_derivative","time_points"))
       .def("FromPointsList",&discretPointToPolynomialC2,
             "Create a piecewise-polynomial connecting exactly all the given points at the given time and respect the given points derivative and second derivative values. The created piecewise is C2 continuous.",args("points","points_derivative","points_second_derivative","time_points"))
      .staticmethod("FromPointsList")
      .def("append",&addFinalPointC0,
           "Append a new polynomial curve of degree 1 at the end of the piecewise curve, defined between self.max() and time and connecting exactly self(self.max()) and end",args("self","end","time"))
       .def("append",&addFinalPointC1,
             "Append a new polynomial curve of degree 3 at the end of the piecewise curve, defined between self.max() and time and connecting exactly self(self.max()) and end. It guarantee C1 continuity and guarantee that self.derivate(time,1) == d_end",args("self","end","d_end","time"))
       .def("append",&addFinalPointC2,
              "Append a new polynomial curve of degree 5 at the end of the piecewise curve, defined between self.max() and time and connecting exactly self(self.max()) and end. It guarantee C2 continuity and guarantee that self.derivate(time,1) == d_end and self.derivate(time,2) == dd_end",args("self","end","d_end","d_end","time"))
      .def("compute_derivate",&piecewise_polynomial_curve_t::compute_derivate,"Return a piecewise_polynomial curve which is the derivate of this.",args("self","order"))
      .def("append", &piecewise_polynomial_curve_t::add_curve,
           "Add a new curve to piecewise curve, which should be defined in T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise_polynomial_curve_t::is_continuous,"Check if the curve is continuous at the given order.")
      .def("saveAsText", &piecewise_polynomial_curve_t::saveAsText<piecewise_polynomial_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_polynomial_curve_t::loadFromText<piecewise_polynomial_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_polynomial_curve_t::saveAsXML<piecewise_polynomial_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_polynomial_curve_t::loadFromXML<piecewise_polynomial_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_polynomial_curve_t::saveAsBinary<piecewise_polynomial_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_polynomial_curve_t::loadFromBinary<piecewise_polynomial_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    class_<piecewise_bezier_curve_t , bases<curve_abc_t> >
    ("piecewise_bezier_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewiseBezierCurveConstructor))
      .def("compute_derivate",&piecewise_polynomial_curve_t::compute_derivate,"Return a piecewise_polynomial curve which is the derivate of this.",args("self","order"))
      .def("add_curve", &piecewise_bezier_curve_t::add_curve)
      .def("is_continuous", &piecewise_bezier_curve_t::is_continuous)
      .def("saveAsText", &piecewise_bezier_curve_t::saveAsText<piecewise_bezier_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_bezier_curve_t::loadFromText<piecewise_bezier_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_bezier_curve_t::saveAsXML<piecewise_bezier_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_bezier_curve_t::loadFromXML<piecewise_bezier_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_bezier_curve_t::saveAsBinary<piecewise_bezier_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_bezier_curve_t::loadFromBinary<piecewise_bezier_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    class_<piecewise_cubic_hermite_curve_t, bases<curve_abc_t> >
    ("piecewise_cubic_hermite_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewiseCubicHermiteCurveConstructor))
      .def("add_curve", &piecewise_cubic_hermite_curve_t::add_curve)
      .def("is_continuous", &piecewise_cubic_hermite_curve_t::is_continuous)
      .def("saveAsText", &piecewise_cubic_hermite_curve_t::saveAsText<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_cubic_hermite_curve_t::loadFromText<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_cubic_hermite_curve_t::saveAsXML<piecewise_cubic_hermite_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_cubic_hermite_curve_t::loadFromXML<piecewise_cubic_hermite_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_cubic_hermite_curve_t::saveAsBinary<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_cubic_hermite_curve_t::loadFromBinary<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    /** END piecewise curve function **/
    /** BEGIN exact_cubic curve**/
    class_<exact_cubic_t, bases<curve_abc_t> >
    ("exact_cubic", init<>())
      .def("__init__", make_constructor(&wrapExactCubicConstructor))
      .def("__init__", make_constructor(&wrapExactCubicConstructorConstraint))
      .def("getNumberSplines", &exact_cubic_t::getNumberSplines)
      .def("getSplineAt", &exact_cubic_t::getSplineAt)
      .def("saveAsText", &exact_cubic_t::saveAsText<exact_cubic_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&exact_cubic_t::loadFromText<exact_cubic_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&exact_cubic_t::saveAsXML<exact_cubic_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&exact_cubic_t::loadFromXML<exact_cubic_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&exact_cubic_t::saveAsBinary<exact_cubic_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&exact_cubic_t::loadFromBinary<exact_cubic_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    /** END exact_cubic curve**/
    /** BEGIN cubic_hermite_spline **/
    class_<cubic_hermite_spline_t, bases<curve_abc_t> >
    ("cubic_hermite_spline", init<>())
      .def("__init__", make_constructor(&wrapCubicHermiteSplineConstructor))
      .def("saveAsText", &cubic_hermite_spline_t::saveAsText<cubic_hermite_spline_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&cubic_hermite_spline_t::loadFromText<cubic_hermite_spline_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&cubic_hermite_spline_t::saveAsXML<cubic_hermite_spline_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&cubic_hermite_spline_t::loadFromXML<cubic_hermite_spline_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&cubic_hermite_spline_t::saveAsBinary<cubic_hermite_spline_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&cubic_hermite_spline_t::loadFromBinary<cubic_hermite_spline_t>,bp::args("filename"),"Loads *this from a binary file.")
      ;

    /** END cubic_hermite_spline **/
    /** BEGIN curve constraints**/
    class_<curve_constraints_t>
        ("curve_constraints", init<>())
            .add_property("init_vel", &get_init_vel, &set_init_vel)
            .add_property("init_acc", &get_init_acc, &set_init_acc)
            .add_property("init_jerk", &get_init_jerk, &set_init_jerk)
            .add_property("end_vel", &get_end_vel, &set_end_vel)
            .add_property("end_acc", &get_end_acc, &set_end_acc)
            .add_property("end_jerk", &get_end_jerk, &set_end_jerk)
        ;
    /** END curve constraints**/
    /** BEGIN bernstein polynomial**/
    class_<bernstein_t>
    ("bernstein", init<const unsigned int, const unsigned int>())
      .def("__call__", &bernstein_t::operator())
    ;
    /** END bernstein polynomial**/

    /** BEGIN curves conversion**/
    def("polynomial_from_bezier", polynomial_from_curve<polynomial_t,bezier_t>);
    def("polynomial_from_hermite", polynomial_from_curve<polynomial_t,cubic_hermite_spline_t>);
    def("bezier_from_hermite", bezier_from_curve<bezier_t,cubic_hermite_spline_t>);
    def("bezier_from_polynomial", bezier_from_curve<bezier_t,polynomial_t>);
    def("hermite_from_bezier", hermite_from_curve<cubic_hermite_spline_t, bezier_t>);
    def("hermite_from_polynomial", hermite_from_curve<cubic_hermite_spline_t, polynomial_t>);
    /** END curves conversion**/

    optimization::python::exposeOptimization();

  } // End BOOST_PYTHON_MODULE
} // namespace curves
