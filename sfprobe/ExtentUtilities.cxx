#include "ExtentUtilities.h"

// C/C++ includes
#include <cassert>

namespace cosmologytools
{

ExtentUtilities::ExtentUtilities()
{
  // TODO Auto-generated constructor stub

}

//------------------------------------------------------------------------------
ExtentUtilities::~ExtentUtilities()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
bool ExtentUtilities::Is3DExtent(INTEGER ext[6])
{
  bool status = false;
  if( (IMAX(ext) > IMIN(ext) ) &&
       (JMAX(ext) > JMIN(ext) ) &&
       (KMAX(ext) > KMIN(ext)) )
   {
   status = true;
   }
 return( status );
}

//------------------------------------------------------------------------------
INTEGER ExtentUtilities::GetLinearIndex( INTEGER ijk[3], INTEGER ext[6])
{
  assert( "pre: given extent must be a 3-D extent" &&
          (ExtentUtilities::Is3DExtent(ext)) );

  // Get local indices staring from 0
  INTEGER li = ijk[0]-IMIN(ext);
  INTEGER lj = ijk[1]-JMIN(ext);
  INTEGER lk = ijk[2]-KMIN(ext);

  // Get Extent dimensions
  INTEGER N2  = JMAX(ext)-JMIN(ext)+1;
  INTEGER N1  = IMAX(ext)-IMIN(ext)+1;

  // Return linear index
  return( (lk*N2+lj)*N1+li );
}

//------------------------------------------------------------------------------
INTEGER ExtentUtilities::GetLinearIndex(
    const INTEGER i, const INTEGER j, const INTEGER k, INTEGER ext[6])
{
  INTEGER ijk[3];
  I(ijk) = i; J(ijk) = j; K(ijk) = k;
  return ExtentUtilities::GetLinearIndex(ijk,ext);
}

//------------------------------------------------------------------------------
void ExtentUtilities::GetStructuredCoordinates(
    INTEGER idx, INTEGER ext[6], INTEGER ijk[3])
{
  assert( "pre: given extent must be a 3-D extent" &&
          (ExtentUtilities::Is3DExtent(ext)) );
  // Get Extent dimensions
  INTEGER N2  = JMAX(ext)-JMIN(ext)+1;
  INTEGER N1  = IMAX(ext)-IMIN(ext)+1;
  INTEGER N12 = N1*N2;
  K(ijk)= idx/N12;
  J(ijk) = (idx-K(ijk)*N12)/N1;
  I(ijk) = idx-K(ijk)*N12-J(ijk)*N1;
}

//------------------------------------------------------------------------------
void ExtentUtilities::GetExtentDimensions(INTEGER ext[6], INTEGER dim[3])
{
  assert("pre: given extent must be a 3-D extent" &&
          ExtentUtilities::Is3DExtent(ext));

  for( int i=0; i < 3; ++i )
    {
    dim[i] = ext[i*2+1]-ext[i*2]+1;
    } // END for all dimensions
}

//------------------------------------------------------------------------------
INTEGER ExtentUtilities::ComputeNumberOfNodes( INTEGER ext[6] )
{
  assert( "pre: given extent must be a 3-D extent" &&
          (ExtentUtilities::Is3DExtent(ext)) );

  int N = 0;
  int idim = IMAX(ext)-IMIN(ext)+1;
  int jdim = JMAX(ext)-JMIN(ext)+1;
  int kdim = KMAX(ext)-KMIN(ext)+1;
  N        = idim * jdim * kdim;
  return( N );
}

//------------------------------------------------------------------------------
INTEGER ExtentUtilities::ComputeNumberOfCells( INTEGER ext[6] )
{
  assert( "pre: given extent must be a 3-D extent" &&
          (ExtentUtilities::Is3DExtent(ext)) );

  int N = 0;
  int idim = IMAX(ext)-IMIN(ext);
  int jdim = JMAX(ext)-JMIN(ext);
  int kdim = KMAX(ext)-KMIN(ext);
  N        = idim * jdim * kdim;
  return( N );
}

} /* namespace cosmologytools */
