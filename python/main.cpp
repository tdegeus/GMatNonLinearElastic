/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATNONLINEARELASTIC_ENABLE_ASSERT

#include <GMatNonLinearElastic/Cartesian3d.h>

namespace py = pybind11;

PYBIND11_MODULE(GMatNonLinearElastic, m)
{

m.doc() = "Non-linear elastic material model";

// --------------------------------
// GMatNonLinearElastic.Cartesian3d
// --------------------------------

py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

namespace SM = GMatNonLinearElastic::Cartesian3d;

// Unit tensors

sm.def("I2", &SM::I2, "Second order unit tensor.");

sm.def("II", &SM::II, "Fourth order tensor with the result of the dyadic product II.");

sm.def("I4", &SM::I4, "Fourth order unit tensor.");

sm.def("I4rt", &SM::I4rt, "Fourth right-transposed order unit tensor.");

sm.def("I4s", &SM::I4s, "Fourth order symmetric projection tensor.");

sm.def("I4d", &SM::I4d, "Fourth order deviatoric projection tensor.");

// Tensor algebra

sm.def("Hydrostatic",
    py::overload_cast<const SM::Tensor2&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns scalar.",
    py::arg("A"));

sm.def("Hydrostatic",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns list of scalars.",
    py::arg("A"));

sm.def("Hydrostatic",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns matrix of scalars.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const SM::Tensor2&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns 2nd-order tensor.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns list 2nd-order tensors.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns matrix 2nd-order tensors.",
    py::arg("A"));

sm.def("Epseq",
    py::overload_cast<const SM::Tensor2&>(&SM::Epseq),
    "Equivalent strain deviator. Returns scalar.",
    py::arg("Eps"));

sm.def("Epseq",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Epseq),
    "Equivalent strain deviator. Returns list of scalars.",
    py::arg("Eps"));

sm.def("Epseq",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Epseq),
    "Equivalent strain deviator. Returns matrix of scalars.",
    py::arg("Eps"));

sm.def("Sigeq",
    py::overload_cast<const SM::Tensor2&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns scalar.",
    py::arg("Sig"));

sm.def("Sigeq",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns list of scalars.",
    py::arg("Sig"));

sm.def("Sigeq",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns matrix of scalars.",
    py::arg("Sig"));

// Material point: NonLinearElastic

py::class_<SM::NonLinearElastic>(sm, "NonLinearElastic")

    .def(py::init<double, double, double, double>(),
        "Non-linear elastic material point",
        py::arg("kappa"),
        py::arg("sig0"),
        py::arg("eps0"),
        py::arg("m"))

    .def("kappa", &SM::NonLinearElastic::kappa, "Returns the bulk modulus.")

    .def("sig0", &SM::NonLinearElastic::sig0, "Returns the reference stress.")

    .def("eps0", &SM::NonLinearElastic::eps0, "Returns the reference strain.")

    .def("m", &SM::NonLinearElastic::m, "Returns the exponent.")

    .def("Stress",
        &SM::NonLinearElastic::Stress,
        "Returns stress tensor, for a given strain tensor.",
        py::arg("Eps"))

    .def("Tangent",
        &SM::NonLinearElastic::Tangent,
        "Returns stress and tangent stiffness tensors, for a given strain tensor.",
        py::arg("Eps"))

    .def("__repr__", [](const SM::NonLinearElastic&) {
        return "<GMatNonLinearElastic.Cartesian3d.NonLinearElastic>";
    });

py::module smm = sm.def_submodule("Type", "Type enumerator");

py::enum_<SM::Type::Value>(smm, "Type")
    .value("Unset", SM::Type::Unset)
    .value("NonLinearElastic", SM::Type::NonLinearElastic)
    .export_values();

// Matrix

py::class_<SM::Matrix>(sm, "Matrix")

    .def(py::init<size_t, size_t>(),
        "Matrix of non-linear elastic material points",
        py::arg("nelem"),
        py::arg("nip"))

    .def(py::init<size_t, size_t, double, double, double, double>(),
        "Matrix of non-linear elastic material points",
        py::arg("nelem"),
        py::arg("nip"),
        py::arg("kappa"),
        py::arg("sig0"),
        py::arg("eps0"),
        py::arg("m"))

    .def("ndim", &SM::Matrix::ndim, "Return number of (tensor) dimensions.")

    .def("nelem", &SM::Matrix::nelem, "Return number of elements (matrix rows).")

    .def("nip", &SM::Matrix::nip, "Return number of integration points (matrix columns).")

    .def("kappa", &SM::Matrix::kappa, "Return matrix with bulk moduli.")

    .def("sig0" &SM::Matrix::sig0, "Return matrix with reference stresses.")

    .def("eps0" &SM::Matrix::eps0, "Return matrix with reference strains.")

    .def("m" &SM::Matrix::m, "Return matrix with exponents.")

    .def("I2", &SM::Matrix::I2, "Return matrix with second order unit tensors.")

    .def("II",
        &SM::Matrix::II,
        "Return matrix with fourth order tensors with the result of the dyadic product II.")

    .def("I4", &SM::Matrix::I4, "Return matrix with fourth order unit tensors.")

    .def("I4rt",
        &SM::Matrix::I4rt,
        "Return matrix with fourth right-transposed order unit tensors.")

    .def("I4s",
        &SM::Matrix::I4s,
        "Return matrix with fourth order symmetric projection tensors.")

    .def("I4d",
        &SM::Matrix::I4d,
        "Return matrix with fourth order deviatoric projection tensors.")

    .def("type", &SM::Matrix::type, "Return matrix with material types.")

    .def("check",
        &SM::Matrix::check,
        "Check that all matrix entries are set. Throws if any unset point is found.")

    .def("setNonLinearElastic",
        &SM::Matrix::setNonLinearElastic,
        py::arg("I"),
        py::arg("kappa"),
        py::arg("sig0"),
        py::arg("eps0"),
        py::arg("m"))

    .def("Stress",
        &SM::Matrix::Stress,
        "Returns matrix of stress tensors, for a given matrix of strain tensors.",
        py::arg("Eps"))

    .def("Tangent",
        &SM::Matrix::Tangent,
        "Returns matrices of stress tangent stiffness tensors, "
        "for a given matrix of strain tensors.",
        py::arg("Eps"))

    .def("__repr__", [](const SM::Matrix&) {
        return "<GMatNonLinearElastic.Cartesian3d.Matrix>";
    });
}
