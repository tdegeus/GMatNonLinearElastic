/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>

#include <GMatNonLinearElastic/version.h>
#include <GMatNonLinearElastic/Cartesian3d.h>

namespace py = pybind11;

template <class S, class T>
auto construct_Array(T& self)
{
    namespace SM = GMatNonLinearElastic::Cartesian3d;

    self.def(py::init<std::array<size_t, S::rank>>(), "Array of material points.", py::arg("shape"))

        .def("shape", &S::shape, "Shape of array.")
        .def("I2", &S::I2, "Array with 2nd-order unit tensors.")
        .def("II", &S::II, "Array with 4th-order tensors = dyadic(I2, I2).")
        .def("I4", &S::I4, "Array with 4th-order unit tensors.")
        .def("I4rt", &S::I4rt, "Array with 4th-order right-transposed unit tensors.")
        .def("I4s", &S::I4s, "Array with 4th-order symmetric projection tensors.")
        .def("I4d", &S::I4d, "Array with 4th-order deviatoric projection tensors.")
        .def("kappa", &S::kappa, "Array with kappa.")
        .def("sig0", &S::sig0, "Array with sig0.")
        .def("eps0", &S::eps0, "Array with eps0.")
        .def("m", &S::m, "Array with m.")
        .def("type", &S::type, "Array with material types.")

        .def(
            "setNonLinearElastic",
            &S::setNonLinearElastic,
            "Set specific entries 'NonLinearElastic'.",
            py::arg("I"),
            py::arg("kappa"),
            py::arg("sig0"),
            py::arg("eps0"),
            py::arg("m"))

        .def(
            "setStrain",
            &S::setStrain,
            "Set strain tensors (computes stress and optionally tangent).",
            py::arg("Eps"),
            py::arg("compute_tangent") = true)

        .def("Strain", &S::Strain, "Get strain tensors.")
        .def("Stress", &S::Stress, "Get stress tensors.")
        .def("Tangent", &S::Tangent, "Get stiffness tensors.")
        .def("getNonLinearElastic", &S::getNonLinearElastic, "Returns underlying NonLinearElastic model.")

        .def("__repr__", [](const S&) { return "<GMatNonLinearElastic.Cartesian3d.Array>"; });
}

template <class S, class T>
void add_deviatoric_overloads(T& module)
{
    module.def(
        "Deviatoric",
        static_cast<S (*)(const S&)>(&GMatNonLinearElastic::Cartesian3d::Deviatoric<S>),
        "Deviatoric part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_hydrostatic_overloads(T& module)
{
    module.def(
        "Hydrostatic",
        static_cast<R (*)(const S&)>(&GMatNonLinearElastic::Cartesian3d::Hydrostatic<S>),
        "Hydrostatic part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_epseq_overloads(T& module)
{
    module.def(
        "Epseq",
        static_cast<R (*)(const S&)>(
            &GMatNonLinearElastic::Cartesian3d::Epseq<S>),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_sigeq_overloads(T& module)
{
    module.def(
        "Sigeq",
        static_cast<R (*)(const S&)>(
            &GMatNonLinearElastic::Cartesian3d::Sigeq<S>),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"));
}


PYBIND11_MODULE(_GMatNonLinearElastic, m)
{
    xt::import_numpy();

    m.doc() = "Non-linear elastic material model";

    m.def("version",
          &GMatNonLinearElastic::version,
          "Return version string.");

    m.def("version_dependencies",
          &GMatNonLinearElastic::version_dependencies,
          "Return list of strings.");

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

    add_deviatoric_overloads<xt::xtensor<double, 4>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 3>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 2>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_epseq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_epseq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_epseq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);

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

        .def(
            "setStrain",
            &SM::NonLinearElastic::setStrain<xt::xtensor<double, 2>>,
            "Set current strain tensor (computes stress and optionally tangent).",
            py::arg("Eps"),
            py::arg("compute_tangent") = true)

        .def("Strain", &SM::NonLinearElastic::Strain, "Returns strain tensor.")
        .def("Stress", &SM::NonLinearElastic::Stress, "Returns stress tensor.")
        .def("Tangent", &SM::NonLinearElastic::Tangent, "Returns tangent stiffness.")

        .def("__repr__", [](const SM::NonLinearElastic&) {
            return "<GMatNonLinearElastic.Cartesian3d.NonLinearElastic>";
        });

    py::module smm = sm.def_submodule("Type", "Type enumerator");

    py::enum_<SM::Type::Value>(smm, "Type")
        .value("Unset", SM::Type::Unset)
        .value("NonLinearElastic", SM::Type::NonLinearElastic)
        .export_values();

    // Array

    py::class_<SM::Array<1>> array1d(sm, "Array1d");
    py::class_<SM::Array<2>> array2d(sm, "Array2d");
    py::class_<SM::Array<3>> array3d(sm, "Array3d");

    construct_Array<SM::Array<1>>(array1d);
    construct_Array<SM::Array<2>>(array2d);
    construct_Array<SM::Array<3>>(array3d);
}
