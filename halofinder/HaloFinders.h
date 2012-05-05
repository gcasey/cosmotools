/**
 * @file HaloFinders.h
 * @brief An enumerator for all available halo-finders.
 */
#ifndef HALOFINDERS_H_
#define HALOFINDERS_H_

namespace cosmologytools {

class HaloFinders {
public:
  enum
  {
  COSMO = 0,
  NUMBER_OF_HALO_FINDERS
  }

};

} /* end namespace hacctools */

#endif /* HALOFINDERS_H_ */
