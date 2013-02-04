#ifndef GENERICIO_H
#define GENERICIO_H

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <stdint.h>

#ifndef GENERICIO_NO_MPI
#include <mpi.h>
#else
#include <fstream>
#endif

namespace cosmotk {

// Forward declrations, defined in the implementation
struct GlobalHeader;
struct VariableHeader;
struct RankHeader;

class GenericIO {
protected:
  struct Variable {
    template <typename T>
    Variable(const std::string &N, T* D, bool ES)
      : Name(N), Size(sizeof(T)),
        IsFloat(!std::numeric_limits<T>::is_integer),
        IsSigned(std::numeric_limits<T>::is_signed),
        Data((void *) D), HasExtraSpace(ES) {}

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
    void *Data;
    bool HasExtraSpace;
  };

public:
#ifndef GENERICIO_NO_MPI
  GenericIO(const MPI_Comm &C, const std::string &FN)
    : NElems(0), Comm(C), FileName(FN), ShouldSwap(false) {}
#else
  GenericIO(const std::string &FN)
    : NElems(0), FileName(FN) {}
#endif

  ~GenericIO() {
    close();
  }

public:
  std::size_t requestedExtraSpace() const {
    return 8;
  }

  void setNumElems(std::size_t E) {
    NElems = E;
  }

  template <typename T>
  void addVariable(const std::string &Name, T *Data,
                   bool HasES = false) {
    Vars.push_back(Variable(Name, Data, HasES));
  }

  template <typename T, typename A>
  void addVariable(const std::string &Name, std::vector<T, A> &Data,
                   bool HasES = false) {
    T *D = Data.empty() ? 0 : &Data[0];
    addVariable(Name, D, HasES);
  }

#ifndef GENERICIO_NO_MPI
  // Writing
  void write();
#endif

  // Reading
  void openAndReadHeader(bool MustMatch = true);

  int readNRanks();
  void readDims(int Dims[3]);
  uint64_t readTotalNumElems();

  struct VariableInfo {
    VariableInfo(const std::string &N, std::size_t S, bool IF, bool IS)
      : Name(N), Size(S), IsFloat(IF), IsSigned(IS) {}

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
  };

  void clearVariables() { this->Vars.clear(); };

  int getNumberOfVariables() { return this->Vars.size(); };

  void getVariableInfo(std::vector<VariableInfo> &VI);

  std::size_t readNumElems(int EffRank = -1);
  void readCoords(int Coords[3], int EffRank = -1);

  void readData(int EffRank = -1, bool PrintStats = true, bool CollStats = true);

  void close() {
    FH.close();
  }

protected:
  std::vector<Variable> Vars;
  std::size_t NElems;
#ifndef GENERICIO_NO_MPI
  MPI_Comm Comm;
#endif
  std::string FileName;
  bool ShouldSwap;

  // Swaps the data in the global header
  void SwapGlobalHeader(GlobalHeader* GH);
  void SwapRankHeader(RankHeader* RH);

  // Swap endian of the given buffer pointing to a memory location of
  // Nb bytes.
  void SwapEndian(void *Addr, const int Nb);

  // This reference counting mechanism allows the the GenericIO class
  // to be used in a cursor mode. To do this, make a copy of the class
  // after reading the header but prior to adding the variables.
  struct FHManager {
#ifndef GENERICIO_NO_MPI
    typedef MPI_File FHType;
#else
    typedef std::fstream FHType;
#endif

    FHManager() : CountedFH(0) {
      allocate();
    }

    FHManager(const FHManager& F) {
      CountedFH = F.CountedFH;
      CountedFH->Cnt += 1;
    }

    ~FHManager() {
      close();
    }

    FHType &get() {
      if (!CountedFH)
        allocate();

      return CountedFH->FH;
    }

    std::vector<char> &getHeaderCache() {
      if (!CountedFH)
        allocate();

      return CountedFH->HeaderCache;
    }

    void allocate() {
      close();
      CountedFH = new FHWCnt;
    };

    void close() {
      if (CountedFH && CountedFH->Cnt == 1)
        delete CountedFH;
      else if (CountedFH)
        CountedFH->Cnt -= 1;

      CountedFH = 0;
    }

    struct FHWCnt {
#ifndef GENERICIO_NO_MPI
      FHWCnt() : FH(MPI_FILE_NULL), Cnt(1) {}
#else
      FHWCnt() : Cnt(1) {}
#endif

      ~FHWCnt() {
        close();
      }

protected:
      void close() {
#ifndef GENERICIO_NO_MPI
        (void) MPI_File_close(&FH);
#else
        FH.close();
#endif
      }

public:
      FHType FH;
      size_t Cnt;


      // Used for reading
      std::vector<char> HeaderCache;
    };

    FHWCnt *CountedFH;
  } FH;
};

}
#endif // GENERICIO_H

