/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

#include "Cartesian3d.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;

// ===================================== GMatNonLinearElastic ======================================

PYBIND11_MODULE(GMatNonLinearElastic, m) {

m.doc() = "Non-linear elastic material model";

// =============================== GMatNonLinearElastic.Cartesian3d ================================

{

// create sub-module
py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

// abbreviate name-space
namespace SM = GMatNonLinearElastic::Cartesian3d;

// -------------------------------------------------------------------------------------------------

sm.def("I"   , &SM::I   );
sm.def("II"  , &SM::II  );
sm.def("I4"  , &SM::I4  );
sm.def("I4rt", &SM::I4rt);
sm.def("I4s" , &SM::I4s );
sm.def("I4d" , &SM::I4d );

// -------------------------------------------------------------------------------------------------

py::class_<SM::NonLinearElastic>(sm, "NonLinearElastic")
  // constructor
  .def(
    py::init<double, double, double, double>(),
    "Non-linear elastic material point",
    py::arg("kappa"),
    py::arg("sig0"),
    py::arg("eps0"),
    py::arg("m")
  )
  // methods
  .def("kappa"  , &SM::NonLinearElastic::kappa)
  .def("sig0"   , &SM::NonLinearElastic::sig0)
  .def("eps0"   , &SM::NonLinearElastic::eps0)
  .def("m"      , &SM::NonLinearElastic::m)
  .def("Sig"    , &SM::NonLinearElastic::Sig)
  .def("Tangent", &SM::NonLinearElastic::Tangent)
  // print to screen
  .def("__repr__", [](const SM::NonLinearElastic &){
    return "<GMatNonLinearElastic.Cartesian3d.NonLinearElastic>"; });

// -------------------------------------------------------------------------------------------------

py::class_<SM::Matrix>(sm, "Matrix")
  // constructor
  .def(
    py::init<size_t, size_t>(),
    "Matrix of non-linear elastic material points",
    py::arg("nelem"),
    py::arg("nip")
  )
  // constructor
  .def(
    py::init<size_t, size_t, double, double, double, double>(),
    "Matrix of non-linear elastic material points",
    py::arg("nelem"),
    py::arg("nip"),
    py::arg("kappa"),
    py::arg("sig0"),
    py::arg("eps0"),
    py::arg("m")
  )
  // methods
  .def("nelem"  , &SM::Matrix::nelem)
  .def("nip"    , &SM::Matrix::nip)
  .def("kappa"  , &SM::Matrix::kappa)
  .def("sig0"   , &SM::Matrix::sig0)
  .def("eps0"   , &SM::Matrix::eps0)
  .def("m"      , &SM::Matrix::m)
  .def("check"  , &SM::Matrix::check)
  .def("set"    , &SM::Matrix::set)
  .def("I"      , &SM::Matrix::I)
  .def("II"     , &SM::Matrix::II)
  .def("I4"     , &SM::Matrix::I4)
  .def("I4rt"   , &SM::Matrix::I4rt)
  .def("I4s"    , &SM::Matrix::I4s)
  .def("I4d"    , &SM::Matrix::I4d)
  .def("Sig"    , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Sig    , py::const_), py::arg("Eps"))
  .def("Tangent", py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Tangent, py::const_), py::arg("Eps"))
  // print to screen
  .def("__repr__", [](const SM::Matrix &){
    return "<GMatNonLinearElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}

