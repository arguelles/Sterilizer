#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/overloads.hpp>
#include "container_conversions.h"
#include "SterileSearch.h"

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

#include <nuSQuIDS/marray.h>

using namespace boost::python;
namespace bp = boost::python;

template<class T>
struct VecToList
{
  static PyObject* convert(const std::vector<T>& vec){
    boost::python::list* l = new boost::python::list();
    for(size_t i =0; i < vec.size(); i++)
      (*l).append(vec[i]);

    return l->ptr();
  }
};

// converting marray to numpy array and back
template<unsigned int DIM>
struct marray_to_numpyarray {
  static PyObject* convert( marray<double,DIM> const & iarray){
    // get the data from the marray
    double * data = iarray.size() ? const_cast<double*>(iarray.get_data()) : static_cast<double*>(NULL);
    // construct numpy object
    npy_intp size[DIM];
    for(unsigned int i = 0; i < DIM; i++)
      size[i] = iarray.extent(i);
    PyArrayObject * pyObj = (PyArrayObject*) PyArray_SimpleNew(DIM,size,PyArray_DOUBLE);
    memcpy(pyObj->data, data, sizeof(double) * iarray.size());

    return PyArray_Return(pyObj);
  }
};

template<typename T,unsigned int DIM>
static marray<T,DIM> numpyarray_to_marray(PyObject * iarray, NPY_TYPES type_num){
  // es un array de numpy
  if (! PyArray_Check(iarray) )
  {
    PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: Input is not a numpy array.");
    boost::python::throw_error_already_set();
  }
  // si es que fuera un array de numpy castearlo
  //PyArrayObject* numpy_array = (PyArrayObject*) iarray;
  // lets get the contiguos C-style array
  PyArrayObject* numpy_array = PyArray_GETCONTIGUOUS((PyArrayObject*)iarray);

  // revisemos que los tipos del array sean dobles o que
  if ( PyArray_DESCR(numpy_array)->type_num != type_num )
  {
    if ( PyArray_DESCR(numpy_array)->type_num == NPY_LONG &&
        PyArray_ITEMSIZE(numpy_array) == 4 && type_num == NPY_INT)
    {
      // numpy on 32 bits sets numpy.int32 to NPY_LONG. So its all ok.
    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: numpy type is not the same as the input array type.");
      boost::python::throw_error_already_set();
    }
  }

  // arrays vacios
  if (PyArray_SIZE(numpy_array) == 0){
      PyErr_SetString(PyExc_TypeError,"numpyarray_to_marray: empty numpy array.");
      boost::python::throw_error_already_set();
  }

  // create numpy iterator
  NpyIter* iter = NpyIter_New(numpy_array, NPY_ITER_READONLY|
                             NPY_ITER_EXTERNAL_LOOP|
                             NPY_ITER_REFS_OK,
                             NPY_KEEPORDER, NPY_NO_CASTING,
                             NULL);

  unsigned int array_dim = PyArray_NDIM(numpy_array);
  assert(DIM == array_dim && "No matching dimensions.");

  // get numpy array shape and create marray object
#ifdef NPY_1_7_API_VERSION
  npy_intp* array_shape = PyArray_SHAPE(numpy_array);
#else
  npy_intp* array_shape = PyArray_DIMS(numpy_array);
#endif
  std::vector<size_t> dimensions;
  for(unsigned int i = 0; i < array_dim; i++)
    dimensions.push_back(array_shape[i]);

  // construct output object
  marray<T,DIM> oarray;
  oarray.resize(dimensions);
  auto it = oarray.begin();

  NpyIter_IterNextFunc *iternext = NpyIter_GetIterNext(iter, NULL);
  char** dataptr = NpyIter_GetDataPtrArray(iter);
  npy_intp* strideptr = NpyIter_GetInnerStrideArray(iter);
  npy_intp* sizeptr = NpyIter_GetInnerLoopSizePtr(iter);
  npy_intp iop, nop = NpyIter_GetNOp(iter);

  // magic to make the int work
  bool magic = false;
  if ( type_num == NPY_INT or type_num == NPY_LONG )
    magic = true;

  do{
    char* data = *dataptr;
    npy_intp count = *sizeptr;
    npy_intp stride = *strideptr;

    while (count--)
    {
      for (iop = 0; iop < nop; ++iop, data+=stride){
        if (magic)
          *it++ = *(T*)(reinterpret_cast<int*>(data));
        else
          *it++ = *(T*)(data);
      }
    }
  } while(iternext(iter));

  NpyIter_Deallocate(iter);

  return oarray;
}

// SterileSearchPy module definitions
/*
static std::vector<double> wrap_llh(LVSearch* lv,std::vector<double> arg){
  if(arg.size() != 3)
    throw std::runtime_error("Number of arguments should be 3. You sent me " + std::to_string(arg.size()));
  std::array<double,3> argv;
  std::copy_n(arg.begin(),3,argv.begin());
  auto fr = lv->llh(argv);
  std::vector<double> result = fr.params;
  result.push_back(fr.likelihood);
  return result;
}

static double wrap_llhFull(LVSearch* lv,std::vector<double> arg){
  if(arg.size() != 9)
    throw std::runtime_error("Number of arguments should be 9. You sent me " + std::to_string(arg.size()));
  std::array<double,9> argv;
  std::copy_n(arg.begin(),9,argv.begin());
  double llh = lv->llhFull(argv);
  return llh;
}

marray<double,3> wrap_GetExpectedDistribution(LVSearch* lv,std::vector<double> arg){
  if(arg.size() != 9)
    throw std::runtime_error("Number of arguments should be 9. You sent me " + std::to_string(arg.size()));
  std::array<double,9> argv;
  std::copy_n(arg.begin(),9,argv.begin());
  auto data = lv->GetExpectationDistribution(argv);
  return data;
}
*/

BOOST_PYTHON_MODULE(SterileSearchPy)
{
  // import numpy array definitions
  import_array();
  import_ufunc();

  class_<DataPaths, boost::noncopyable,std::shared_ptr<DataPaths> >("DataPaths",init<>())
    .def_readwrite("compact_file_path",&DataPaths::compact_file_path)
    .def_readwrite("squids_files_path",&DataPaths::squids_files_path)
    .def_readwrite("prompt_squids_files_path",&DataPaths::prompt_squids_files_path)
    .def_readwrite("xs_spline_path",&DataPaths::xs_spline_path)
    .def_readwrite("data_path",&DataPaths::data_path)
    .def_readwrite("mc_path",&DataPaths::mc_path)
    .def_readwrite("oversize_path",&DataPaths::oversize_path)
    .def_readwrite("domeff_spline_path",&DataPaths::domeff_spline_path)
    .def_readwrite("flux_splines_path",&DataPaths::flux_splines_path)
  ;

  class_<SteeringParams, boost::noncopyable,std::shared_ptr<SteeringParams> >("SteeringParams",init<>())
    .def_readwrite("minFitEnergy",&SteeringParams::minFitEnergy)
    .def_readwrite("maxFitEnergy",&SteeringParams::maxFitEnergy)
    .def_readwrite("minCosth",&SteeringParams::minCosth)
    .def_readwrite("maxCosth ",&SteeringParams::maxCosth)
    .def_readwrite("logEbinEdge ",&SteeringParams::logEbinEdge)
    .def_readwrite("logEbinWidth ",&SteeringParams::logEbinWidth)
    .def_readwrite("cosThbinEdge ",&SteeringParams::cosThbinEdge)
    .def_readwrite("cosThbinWidth",&SteeringParams::cosThbinWidth)
    .def_readwrite("useFactorization",&SteeringParams::useFactorization)
    .def_readwrite("useBurnSample",&SteeringParams::useBurnSample)
    .def_readwrite("simToLoad",&SteeringParams::simToLoad)
    .def_readwrite("quiet",&SteeringParams::quiet)
    .def_readwrite("evalThreads",&SteeringParams::evalThreads)
    .def_readwrite("modelName",&SteeringParams::modelName)
    .def_readwrite("oversizeFunction",&SteeringParams::oversizeFunction)
    .def_readwrite("ReadCompact",&SteeringParams::ReadCompact)
    .def_readwrite("xs_model_name",&SteeringParams::xs_model_name)
    .def_readwrite("years",&SteeringParams::years)
    .def_readwrite("fullLivetime",&SteeringParams::fullLivetime)
    .def_readwrite("burnSampleLivetime",&SteeringParams::burnSampleLivetime)
  ;

  class_<SterileNuParams, boost::noncopyable,std::shared_ptr<SterileNuParams> >("SterileNuParams",init<>())
    .def_readwrite("modelId",&SterileNuParams::modelId)
    .def_readwrite("th14",&SterileNuParams::th14)
    .def_readwrite("th24",&SterileNuParams::th24)
    .def_readwrite("th34",&SterileNuParams::th34)
    .def_readwrite("del14",&SterileNuParams::del14)
    .def_readwrite("del24",&SterileNuParams::del24)
    .def_readwrite("dm41sq",&SterileNuParams::dm41sq)
  ;

  class_<Sterilizer, boost::noncopyable, std::shared_ptr<Sterilizer> >("Sterilizer", init<DataPaths,SteeringParams,SterileNuParams>())
    .def("CheckDataLoaded",&Sterilizer::CheckDataLoaded)
    .def("CheckSimulationLoaded",&Sterilizer::CheckSimulationLoaded)
    .def("CheckDOMEfficiencySplinesConstructed",&Sterilizer::CheckDOMEfficiencySplinesConstructed)
    .def("CheckCrossSectionWeighterConstructed",&Sterilizer::CheckCrossSectionWeighterConstructed)
    .def("CheckFluxWeighterConstructed",&Sterilizer::CheckFluxWeighterConstructed)
    .def("CheckLeptonWeighterConstructed",&Sterilizer::CheckLeptonWeighterConstructed)
    .def("CheckDataHistogramConstructed",&Sterilizer::CheckDataHistogramConstructed)
    .def("CheckLeptonWeighterConstructed",&Sterilizer::CheckLeptonWeighterConstructed)
    .def("CheckDataHistogramConstructed",&Sterilizer::CheckDataHistogramConstructed)
    .def("CheckSimulationHistogramConstructed",&Sterilizer::CheckSimulationHistogramConstructed)
    .def("CheckLikelihoodProblemConstruction",&Sterilizer::CheckLikelihoodProblemConstruction)
    .def("CheckDataPaths",&Sterilizer::CheckDataPaths)
    .def("ReportStatus",&Sterilizer::ReportStatus)
    .def("GetDataDistribution",&Sterilizer::GetDataDistribution)
    .def("GetExpectation",(marray<double,3>(Sterilizer::*)(Nuisance)const)&Sterilizer::GetExpectation)
    .def("GetRealization",(marray<double,3>(Sterilizer::*)(Nuisance,int)const)&Sterilizer::GetRealization)
    .def("EvalLLH",&Sterilizer::EvalLLH)
    .def("MinLLH",&Sterilizer::MinLLH)
    .def("SetSterileNuParams",&Sterilizer::SetSterileNuParams)
  ;

  // python container to vector<double> convertion
  using namespace scitbx::boost_python::container_conversions;
  from_python_sequence< std::vector<double>, variable_capacity_policy >();
  to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();
  to_python_converter< marray<double,1> , marray_to_numpyarray<1> >();
  to_python_converter< marray<double,2> , marray_to_numpyarray<2> >();
  to_python_converter< marray<double,3> , marray_to_numpyarray<3> >();
  to_python_converter< marray<double,4> , marray_to_numpyarray<4> >();
}
