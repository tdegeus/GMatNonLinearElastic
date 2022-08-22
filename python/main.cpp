/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>
#include <xtensor-python/xtensor_python_config.hpp> // todo: remove for xtensor-python >0.26.1

#define GMATELASTIC_USE_XTENSOR_PYTHON
#define GMATNONLINEARELASTIC_USE_XTENSOR_PYTHON
#define GMATTENSOR_USE_XTENSOR_PYTHON
#include <GMatElastic/Cartesian3d.h>
#include <GMatElastic/version.h>
#include <GMatNonLinearElastic/Cartesian3d.h>
#include <GMatNonLinearElastic/version.h>
#include <GMatTensor/Cartesian3d.h>

namespace py = pybind11;

namespace my3d {

template <class S, class T>
auto NonLinearElastic(T& cls)
{
    cls.def(
        py::init<
            const xt::pytensor<double, S::rank>&,
            const xt::pytensor<double, S::rank>&,
            const xt::pytensor<double, S::rank>&,
            const xt::pytensor<double, S::rank>&>(),
        "Heterogeneous system.",
        py::arg("kappa"),
        py::arg("sig0"),
        py::arg("eps0"),
        py::arg("m"));

    cls.def_property_readonly("shape", &S::shape, "Shape of array.");
    cls.def_property_readonly("shape_tensor2", &S::shape_tensor2, "Array of rank 2 tensors.");
    cls.def_property_readonly("shape_tensor4", &S::shape_tensor4, "Array of rank 4 tensors.");
    cls.def_property_readonly("kappa", &S::kappa, "Bulk modulus.");
    cls.def_property_readonly("sig0", &S::sig0, "Reference stress.");
    cls.def_property_readonly("eps0", &S::eps0, "Reference strain.");
    cls.def_property_readonly("m", &S::m, "Exponent.");
    cls.def_property_readonly("Sig", &S::Sig, "Stress tensor.");
    cls.def_property_readonly("C", &S::C, "Tangent tensor.");

    cls.def_property(
        "Eps",
        static_cast<xt::pytensor<double, S::rank + 2>& (S::*)()>(&S::Eps),
        static_cast<void (S::*)(const xt::pytensor<double, S::rank + 2>&)>(&S::set_Eps),
        "Strain tensor");

    cls.def(
        "set_Eps",
        py::overload_cast<const xt::pytensor<double, S::rank + 2>&, bool>(
            &S::template set_Eps<xt::pytensor<double, S::rank + 2>>),
        "Overwrite strain tensor.",
        py::arg("arg"),
        py::arg("compute_tangent") = true);

    cls.def(
        "refresh", &S::refresh, "Recompute stress from strain.", py::arg("compute_tangent") = true);

    cls.def(
        "__repr__", [](const S&) { return "<GMatNonLinearElastic.Cartesian3d.NonLinearElastic>"; });
}

template <class R, class T, class M>
void Epseq(M& mod)
{
    mod.def(
        "Epseq",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Epseq),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void epseq(M& mod)
{
    mod.def(
        "epseq",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::epseq),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

template <class R, class T, class M>
void Sigeq(M& mod)
{
    mod.def(
        "Sigeq",
        static_cast<R (*)(const T&)>(&GMatElastic::Cartesian3d::Sigeq),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class T, class M>
void sigeq(M& mod)
{
    mod.def(
        "sigeq",
        static_cast<void (*)(const T&, R&)>(&GMatElastic::Cartesian3d::sigeq),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"),
        py::arg("ret"));
}

} // namespace my3d

/**
Overrides the `__name__` of a module.
Classes defined by pybind11 use the `__name__` of the module as of the time they are defined,
which affects the `__repr__` of the class type objects.
*/
class ScopedModuleNameOverride {
public:
    explicit ScopedModuleNameOverride(py::module m, std::string name) : module_(std::move(m))
    {
        original_name_ = module_.attr("__name__");
        module_.attr("__name__") = name;
    }
    ~ScopedModuleNameOverride()
    {
        module_.attr("__name__") = original_name_;
    }

private:
    py::module module_;
    py::object original_name_;
};

PYBIND11_MODULE(_GMatNonLinearElastic, m)
{
    ScopedModuleNameOverride name_override(m, "GMatNonLinearElastic");

    xt::import_numpy();

    m.doc() = "Elasto-plastic material model";
    m.def("version", &GMatElastic::version, "Return version string.");

    m.def(
        "version_dependencies",
        &GMatElastic::version_dependencies,
        "List of version strings, include dependencies.");

    // --------------------------------
    // GMatNonLinearElastic.Cartesian3d
    // --------------------------------

    py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

    namespace SM = GMatNonLinearElastic::Cartesian3d;

    // Tensor algebra

    my3d::Epseq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    my3d::Epseq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    my3d::Epseq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    my3d::epseq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    my3d::epseq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    my3d::epseq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    my3d::Sigeq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    my3d::Sigeq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    my3d::Sigeq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    my3d::sigeq<xt::pytensor<double, 2>, xt::pytensor<double, 4>>(sm);
    my3d::sigeq<xt::pytensor<double, 1>, xt::pytensor<double, 3>>(sm);
    my3d::sigeq<xt::pytensor<double, 0>, xt::pytensor<double, 2>>(sm);

    // NonLinearElastic

    py::class_<SM::NonLinearElastic<0>, GMatTensor::Cartesian3d::Array<0>> array0d(
        sm, "NonLinearElastic0d");

    py::class_<SM::NonLinearElastic<1>, GMatTensor::Cartesian3d::Array<1>> array1d(
        sm, "NonLinearElastic1d");

    py::class_<SM::NonLinearElastic<2>, GMatTensor::Cartesian3d::Array<2>> array2d(
        sm, "NonLinearElastic2d");

    py::class_<SM::NonLinearElastic<3>, GMatTensor::Cartesian3d::Array<3>> array3d(
        sm, "NonLinearElastic3d");

    my3d::NonLinearElastic<SM::NonLinearElastic<0>>(array0d);
    my3d::NonLinearElastic<SM::NonLinearElastic<1>>(array1d);
    my3d::NonLinearElastic<SM::NonLinearElastic<2>>(array2d);
    my3d::NonLinearElastic<SM::NonLinearElastic<3>>(array3d);
}
