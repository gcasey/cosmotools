#ifndef COSMOLOGYTOOLSMACROS_H_
#define COSMOLOGYTOOLSMACROS_H_

// Select floating point precision
#ifdef USE_DOUBLE_PRECISION
 #define REAL double
#else
 #define REAL float
#endif

// Select integer precision
#ifdef USECPP0X
  #include <cstdint>
  #ifdef USEINT64
   #define INTEGER int32_t
  #else
   #define INTEGER int64_t
  #endif
#else
  #define INTEGER int
#endif

/**
 * @brief A macro to disable the copy constructor and assignment of a class.
 * @note Must be used within the private block of the class definition.
 */
#define DISABLE_COPY_AND_ASSIGNMENT(ClassName)    \
  ClassName(const ClassName&);                    \
  void operator=(const ClassName&);

/**
 * @brief A macro that expands a getter routine for the given ivar.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define GetMacro( ivar, type )       \
    virtual type Get##ivar() const { \
      return this->ivar;             \
    }

/**
 * @brief A macro that expands a setter routine for the given ivar.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define SetMacro( ivar, type )            \
    virtual void Set##ivar( type __arg) { \
        this->ivar = __arg;               \
    }

/**
 * @brief Same as SetMacro but, ivar is a vector of size 3.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define SetVector3Macro( ivar, type )      \
    virtual void Set##ivar( type* __arg ){ \
        this->ivar[0] = __arg[0];          \
        this->ivar[1] = __arg[1];          \
        this->ivar[2] = __arg[2];          \
    }

/**
 * @brief Same as SetMacro but, ivar is a vector of size 6.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define SetVector6Macro( ivar, type )      \
    virtual void Set##ivar( type* __arg ){ \
        this->ivar[0] = __arg[0];          \
        this->ivar[1] = __arg[1];          \
        this->ivar[2] = __arg[2];          \
        this->ivar[3] = __arg[3];          \
        this->ivar[4] = __arg[4];          \
        this->ivar[5] = __arg[5];          \
    }

/**
 * @brief A macro that expands a setter & getter routines for the given ivar.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define GetNSetMacro(ivar, type)          \
    virtual type Get##ivar() const {      \
      return this->ivar;                  \
    }                                     \
    virtual void Set##ivar( type __arg) { \
        this->ivar = __arg;               \
    }


#endif /* HACCTOOLSMACROS_H_ */