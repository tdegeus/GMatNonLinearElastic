/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATNONLINEARELASTIC_ENABLE_ASSERT

// include library
#include "../include/GMatNonLinearElastic/Cartesian3d.h"

// abbreviate name-space
namespace py = pybind11;

// -------------------------------------------------------------------------------------------------

PYBIND11_MODULE(GMatNonLinearElastic, m) {

m.doc() = "Non-linear elastic material model";

// create submodule
py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

// abbreviate name-space
namespace SM = GMatNonLinearElastic::Cartesian3d;

// abbreviate types(s)
typedef SM::T2 T2;

// -------------------------------------------------------------------------------------------------

sm.def("I", &SM::I);
sm.def("II", &SM::II);
sm.def("I4", &SM::I4);
sm.def("I4rt", &SM::I4rt);
sm.def("I4s", &SM::I4s);
sm.def("I4d", &SM::I4d);

// -------------------------------------------------------------------------------------------------

sm.def("Hydrostatic",
  py::overload_cast<const T2&>(&SM::Hydrostatic),
  "Hydrostatic part of a 2nd-order tensor",
  py::arg("A"));

sm.def("Deviatoric",
  py::overload_cast<const T2&>(&SM::Deviatoric),
  "Deviatoric",
  py::arg("A"));

sm.def("Epseq",
  py::overload_cast<const T2&>(&SM::Epseq),
  "Equivalent strain deviator",
  py::arg("Eps"));

sm.def("Sigeq",
  py::overload_cast<const T2&>(&SM::Sigeq),
  "Equivalent stress deviator",
  py::arg("Sig"));

// -------------------------------------------------------------------------------------------------

sm.def("Hydrostatic",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Hydrostatic),
  "Hydrostatic part of a 2nd-order tensor",
  py::arg("A"));

sm.def("Deviatoric",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Deviatoric),
  "Deviatoric",
  py::arg("A"));

sm.def("Epseq",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Epseq),
  "Equivalent strain deviator",
  py::arg("Eps"));

sm.def("Sigeq",
  py::overload_cast<const xt::xtensor<double,4>&>(&SM::Sigeq),
  "Equivalent stress deviator",
  py::arg("Sig"));

// -------------------------------------------------------------------------------------------------

py::class_<SM::NonLinearElastic>(sm, "NonLinearElastic")

  .def(
    py::init<double, double, double, double>(),
    "Non-linear elastic material point",
    py::arg("kappa"),
    py::arg("sig0"),
    py::arg("eps0"),
    py::arg("m")
  )

  .def("kappa", &SM::NonLinearElastic::kappa)
  .def("sig0", &SM::NonLinearElastic::sig0)
  .def("eps0", &SM::NonLinearElastic::eps0)
  .def("m", &SM::NonLinearElastic::m)
  .def("Stress", &SM::NonLinearElastic::Stress, py::arg("Eps"))
  .def("Tangent", &SM::NonLinearElastic::Tangent, py::arg("Eps"))

  .def("__repr__", [](const SM::NonLinearElastic &){
    return "<GMatNonLinearElastic.Cartesian3d.NonLinearElastic>"; });

// -------------------------------------------------------------------------------------------------

py::module smm = sm.def_submodule("Type", "Type enumerator");

py::enum_<SM::Type::Value>(smm, "Type")
    .value("Unset", SM::Type::Unset)
    .value("NonLinearElastic", SM::Type::NonLinearElastic)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<SM::Matrix>(sm, "Matrix")

  .def(
    py::init<size_t, size_t>(),
    "Matrix of non-linear elastic material points",
    py::arg("nelem"),
    py::arg("nip")
  )

  .def(
    py::init<size_t, size_t, double, double, double, double>(),
    "Matrix of non-linear elastic material points",
    py::arg("nelem"),
    py::arg("nip"),
    py::arg("kappa"),
    py::arg("sig0"),
    py::arg("eps0"),
    py::arg("m"))

  .def("ndim", &SM::Matrix::ndim)
  .def("nelem", &SM::Matrix::nelem)
  .def("nip", &SM::Matrix::nip)

  .def("kappa", &SM::Matrix::kappa)
  .def("sig0", &SM::Matrix::sig0)
  .def("eps0", &SM::Matrix::eps0)
  .def("m", &SM::Matrix::m)

  .def("I", &SM::Matrix::I)
  .def("II", &SM::Matrix::II)
  .def("I4", &SM::Matrix::I4)
  .def("I4rt", &SM::Matrix::I4rt)
  .def("I4s", &SM::Matrix::I4s)
  .def("I4d", &SM::Matrix::I4d)

  .def("check", &SM::Matrix::check)

  .def("setNonLinearElastic",
    &SM::Matrix::setNonLinearElastic,
    py::arg("I"),
    py::arg("kappa"),
    py::arg("sig0"),
    py::arg("eps0"),
    py::arg("m"))

  .def("Stress",
    py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Stress, py::const_),
    py::arg("Eps"))

  .def("Tangent",
    py::overload_cast<const xt::xtensor<double,4> &>(&SM::Matrix::Tangent, py::const_),
    py::arg("Eps"))

  .def("__repr__", [](const SM::Matrix &){
    return "<GMatNonLinearElastic.Cartesian3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

