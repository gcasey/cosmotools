/**
 * @class UniformProber
 * @brief A class that provides functionality for probing the StructureFormation
    with a uniform grid of points
 */

#ifndef UNIFORMPROBER_H_
#define UNIFORMPROBER_H_

#include "CosmoToolsMacros.h"

namespace cosmologytools {

//forward declare
class StructureFormationProbe;

class UniformProber
{
public:
  //probe the structure formation with the given uniform grid information
  UniformProber( REAL origin[3], REAL spacing[3], INTEGER ext[6] );

  virtual ~UniformProber();

  /**
   * @brief Executes the probing on the structureFormationProbe
      using the all the points in the uniform grid in serial
   */
  void RunSerialProber(cosmologytools::StructureFormationProbe * probe,
                       int timestep);

  /**
   * @brief Executes the probing on the structureFormationProbe
      using the all the points in the uniform grid in parallel
      using cuda or tbb depending on how the project was built.

      If Dax wasn't enabled during configuration this is a no op.
   */
  void RunDaxProber(cosmologytools::StructureFormationProbe * probe,
                    int timestep);

  /**
   * @brief Get macro for Number of Streams
   */
  virtual INTEGER* GetNumberOfStreams(){ return this->NumberOfStreams; }

  /**
   * @brief Get macro for Rho values
   */
  virtual REAL* GetRho(){ return this->Rho; }

  /**
   * @brief Get the number of points that the uniform prober has
   */
  GetMacro(NumberOfPoints,INTEGER);

  REAL* GetOrigin(){ return &(*Origin); }
  REAL* GetSpacing(){ return &(*Spacing); }
  INTEGER* GetExtents(){ return &(*Extents); }

protected:
  REAL Origin[3];
  REAL Spacing[3];
  INTEGER Extents[6];
  INTEGER NumberOfPoints;

  INTEGER* NumberOfStreams;
  REAL* Rho;

private:
  void ComputePoint(INTEGER index, REAL xyz[3]);


  DISABLE_COPY_AND_ASSIGNMENT(UniformProber);
};

}

#endif
