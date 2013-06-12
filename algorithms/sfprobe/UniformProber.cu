
//We are just fowarding to uniform grid prober implementation
// but we need to set the device adapter first to be uniform prober

#  define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_CUDA
#  define BOOST_SP_DISABLE_THREADS
#include "UniformProber.cxx"
