#ifndef HALOTYPE_H_
#define HALOTYPE_H_

#include "CosmologyToolsMacros.h"

#include <sstream>

/**
 * @brief Corresponding string to halotypes
 */
static const char *HaloTypeName[] = {
    "<GHOST>",
    "<ZOMBIE>"
};

namespace cosmotk {


class HaloType
{
public:

  enum {
    GHOST,  /* halos that are owned by a remote (neighboring) process */
    ZOMBIE, /* halos that died in a previous timestep, but, are kept around */

    NUMBER_OF_HALO_TYPES
  } HaloTypes;

  /**
   * @brief Sets the type in the given type mask.
   * @param mask the user-supplied type mask that is modified.
   * @param type the type to set.
   * @post HaloType::IsType(mask,type)==true.
   */
  static void SetType(unsigned char &mask, const int type)
  {
    assert("pre: invalid type!" &&
           (type >= 0) && (type < HaloType::NUMBER_OF_HALO_TYPES) );
    assert("pre: insufficient number of bits to encode halo types!" &&
            (HaloType::NUMBER_OF_HALO_TYPES <= 8));
    mask |= (1 << type);
  }

  /**
   * @brief Unsets the type in the given type mask.
   * @param mask the user-supplied type mask to modify.
   * @param type the type to unset.
   * @post HaloType::IsType(mask,type)==false.
   */
  static void UnsetType(unsigned char &mask, const int type)
  {
    assert("pre: invalid type!" &&
           (type >= 0) && (type < 8) );
    assert("pre: insufficient number of bits to encode halo types!" &&
            (HaloType::NUMBER_OF_HALO_TYPES <= 8));
    mask &= ~(1 << type);
  }

  /**
   * @brief Check if the type in query is set in the given type mask.
   * @param mask the type mask to check.
   * @param type the type to check.
   * @return status true iff the type is set, else, false.
   */
  static bool IsType(unsigned char &mask, const int type)
  {
    assert("pre: invalid type!" &&
           (type >= 0) && (type < HaloType::NUMBER_OF_HALO_TYPES) );
    assert("pre: insufficient number of bits to encode halo types!" &&
            (HaloType::NUMBER_OF_HALO_TYPES <= 8));
    bool status = ( mask & (1 << type))? true : false;
    return( status );
  }

  /**
   * @brief Resets all bits in the type mask.
   * @param mask the mask to reset.
   * @post HaloType::IsType(mask,i)==false \f$ \forall i \in \[0,7\] \f$
   */
  static void Reset(unsigned char &mask)
  {
    for(int i=0; i < 8; ++i)
      {
      HaloType::UnsetType(mask,i);
      } // END for all bits in the mask
  }

  /**
   * @brief Returns a string representation for the given typemask.
   * @param mask the type mask in question.
   * @return s a string corresponding to the different types encoded in the
   * typemask.
   */
  static std::string GetTypeString(unsigned char mask)
  {
    std::ostringstream oss;
    for(int i=0; i < HaloType::NUMBER_OF_HALO_TYPES; ++i)
      {
      if(HaloType::IsType(mask,i))
        {
        oss << HaloTypeName[ i ];
        }
      } // END for all halos
    return( oss.str() );
  }

protected:
  HaloType();
  ~HaloType();

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloType);

};

}



#endif /* HALOTYPE_H_ */
