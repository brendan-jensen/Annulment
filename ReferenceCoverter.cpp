#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <openbabel/obchargemodel.h>

static PyObject *
get_charge_mol(PyObject *self, PyObject *args)
{
    std::cout << "Hello World" << std::endl;
}