/**
 * @brief Base class for Generic I/O that encapsulates common functionality
 * for readers and writers.
 */
#ifndef GENERICIOBASE_H_
#define GENERICIOBASE_H_

#include "CosmologyToolsMacros.h"

// STL includes
#include <vector>
#include <limits>

// MPI include
#include <mpi.h>

namespace cosmotk
{


struct VariableInfo
{
  VariableInfo(const std::string &N, std::size_t S, bool IF, bool IS,
               bool PCX, bool PCY, bool PCZ, bool PG)
    : Name(N), Size(S), IsFloat(IF), IsSigned(IS),
      IsPhysCoordX(PCX), IsPhysCoordY(PCY), IsPhysCoordZ(PCZ),
      MaybePhysGhost(PG) {}

  std::string Name;
  std::size_t Size;
  bool IsFloat;
  bool IsSigned;
  bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
  bool MaybePhysGhost;
};


class GenericIOBase
{
public:

  enum VariableFlags {
    VarHasExtraSpace =  (1 << 0),
    VarIsPhysCoordX  =  (1 << 1),
    VarIsPhysCoordY  =  (1 << 2),
    VarIsPhysCoordZ  =  (1 << 3),
    VarMaybePhysGhost = (1 << 4)
  };

  enum FileIOStrategy {
    FileIOMPI,
    FileIOPOSIX,
    FileIOMPICollective,
    FileIOUndefined
  };

  GenericIOBase();
  virtual ~GenericIOBase();

  // Get and Set macros
  GetNSetMacro(FileName,std::string);
  GetNSetMacro(Comm,MPI_Comm);

  /**
   * @brief Adds a variable associated with given name and data array.
   * @param Name the name of the variable
   * @param Data the user-supplied array consisting of the data
   * @param Flags bit mask associated with this variable
   */
  template <typename T>
  void AddVariable(const std::string &Name, T *Data, unsigned Flags=0)
    {this->Vars.push_back(Variable(Name,Data,Flags)); }

  /**
   * @brief Adds a variable associated with the given name and STL vector.
   * @param Name the name of the variable
   * @param Data the STL vector with the data associated with the variable.
   * @param Flags bit mask associated with this variable
   */
  template <typename T, typename A>
  void AddVariable(const std::string &Name, std::vector<T,A> &Data,
                      unsigned Flags=0)
    {
    T *D = Data.empty() ? NULL : &Data[0];
    this->Vars.push_back( Variable(Name,D,Flags) );
    }

  /**
   * @brief Adds a variable associated with the given user-supplied
   * VariableInfo object and the data pointed by the void* pointer.
   * @param VI the user-supplied variable information object
   * @param Data pointer to the data corresponding to the data
   * @param Flags bit mask associated with this variable
   * @see VariableInfo
   */
  void AddVariable(const VariableInfo &VI, void *Data, unsigned Flags=0)
    { this->Vars.push_back( Variable(VI,Data,Flags) ); }

  /**
   * @brief Return the number of variables associated with this instance.
   * @return N the number of variable associated with this instance.
   */
  int GetNumberOfVariables()
    { return this->Vars.size(); }

  /**
   * @brief Clears all variables associated with this instance.
   * @post this->GetNumberOfVariables()==0.
   */
  void ClearVariables()
    { this->Vars.clear(); }

protected:

  struct Variable
  {
    template <typename T>
    Variable(const std::string &N, T* D, unsigned Flags = 0)
      : Name(N), Size(sizeof(T)),
        IsFloat(!std::numeric_limits<T>::is_integer),
        IsSigned(std::numeric_limits<T>::is_signed),
        Data((void *) D), HasExtraSpace(Flags & VarHasExtraSpace),
        IsPhysCoordX(Flags & VarIsPhysCoordX),
        IsPhysCoordY(Flags & VarIsPhysCoordY),
        IsPhysCoordZ(Flags & VarIsPhysCoordZ),
        MaybePhysGhost(Flags & VarMaybePhysGhost) {}

    Variable(const VariableInfo &VI, void *D, unsigned Flags = 0)
      : Name(VI.Name), Size(VI.Size), IsFloat(VI.IsFloat),
        IsSigned(VI.IsSigned), Data(D),
        HasExtraSpace(Flags & VarHasExtraSpace),
        IsPhysCoordX((Flags & VarIsPhysCoordX) || VI.IsPhysCoordX),
        IsPhysCoordY((Flags & VarIsPhysCoordY) || VI.IsPhysCoordY),
        IsPhysCoordZ((Flags & VarIsPhysCoordZ) || VI.IsPhysCoordZ),
        MaybePhysGhost((Flags & VarMaybePhysGhost) || VI.MaybePhysGhost) {}

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
    void *Data;
    bool HasExtraSpace;
    bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
    bool MaybePhysGhost;
  }; // END definition of Variable

  unsigned IOStrategy;
  std::string FileName;
  std::vector<Variable> Vars;
  MPI_Comm Comm;

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOBase);
}; // END GenericIOBase definition


} /* namespace cosmotk */
#endif /* GENERICIOBASE_H_ */
