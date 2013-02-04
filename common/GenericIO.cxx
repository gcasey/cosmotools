#define _XOPEN_SOURCE 600
#include "CRC64.h"
#include "GenericIO.h"
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <cstring>

#ifndef GENERICIO_NO_MPI
#include <ctime>
#endif

using namespace std;
namespace cosmotk {

static bool isBigEndian() {
  const uint32_t one = 1;
  return !(*((char *)(&one)));
}

static const size_t CRCSize = 8;

static const size_t MagicSize = 8;
static const char *MagicBE = "HACC01B";
static const char *MagicLE = "HACC01L";

struct GlobalHeader {
  char Magic[MagicSize];
  uint64_t HeaderSize;
  uint64_t NElems; // The global total
  uint64_t Dims[3];
  uint64_t NVars;
  uint64_t VarsSize;
  uint64_t VarsStart;
  uint64_t NRanks;
  uint64_t RanksSize;
  uint64_t RanksStart;
  uint64_t GlobalHeaderSize;
} __attribute__((packed));

enum {
  FloatValue =  (1 << 0),
  SignedValue = (1 << 1)
};

static const size_t NameSize = 256;
struct VariableHeader {
  char Name[NameSize];
  uint64_t Flags;
  uint64_t Size;
} __attribute__((packed));

struct RankHeader {
  uint64_t Coords[3];
  uint64_t NElems;
  uint64_t Start;
} __attribute__((packed));

#ifndef GENERICIO_NO_MPI
// Note: writing errors are not currently recoverable (one rank may fail
// while the others don't).
void GenericIO::write() {
  const char *Magic = isBigEndian() ? MagicBE : MagicLE;

  uint64_t FileSize = 0;

  int NRanks, Rank;
  MPI_Comm_rank(Comm, &Rank);
  MPI_Comm_size(Comm, &NRanks);

  RankHeader RHLocal;
  int Dims[3], Periods[3], Coords[3];
  MPI_Cart_get(Comm, 3, Dims, Periods, Coords);
  std::copy(Coords, Coords + 3, RHLocal.Coords);
  RHLocal.NElems = NElems;
  RHLocal.Start = 0;

  double StartTime = MPI_Wtime();

  if (Rank == 0) {
    uint64_t HeaderSize = sizeof(GlobalHeader) + Vars.size()*sizeof(VariableHeader) +
                          NRanks*sizeof(RankHeader) + CRCSize;
    vector<char> Header(HeaderSize, 0);
    GlobalHeader *GH = (GlobalHeader *) &Header[0];
    std::copy(Magic, Magic + MagicSize, GH->Magic);
    GH->HeaderSize = HeaderSize - CRCSize;
    GH->NElems = NElems; // This will be updated later
    std::copy(Dims, Dims + 3, GH->Dims);
    GH->NVars = Vars.size();
    GH->VarsSize = sizeof(VariableHeader);
    GH->VarsStart = sizeof(GlobalHeader);
    GH->NRanks = NRanks;
    GH->RanksSize = sizeof(RankHeader);
    GH->RanksStart = GH->VarsStart + Vars.size()*sizeof(VariableHeader);
    GH->GlobalHeaderSize = sizeof(GlobalHeader);

    uint64_t RecordSize = 0;
    VariableHeader *VH = (VariableHeader *) &Header[GH->VarsStart];
    for (size_t i = 0; i < Vars.size(); ++i, ++VH) {
      string VName(Vars[i].Name);
      VName.resize(NameSize);

      std::copy(VName.begin(), VName.end(), VH->Name);
      if (Vars[i].IsFloat)  VH->Flags |= FloatValue;
      if (Vars[i].IsSigned) VH->Flags |= SignedValue;
      RecordSize += VH->Size = Vars[i].Size;
    }

    MPI_Gather(&RHLocal, sizeof(RHLocal), MPI_BYTE,
               &Header[GH->RanksStart], sizeof(RHLocal),
               MPI_BYTE, 0, Comm);

    RankHeader *RH = (RankHeader *) &Header[GH->RanksStart];
    RH->Start = HeaderSize; ++RH;
    for (int i = 1; i < NRanks; ++i, ++RH) {
      uint64_t PrevNElems = RH[-1].NElems;
      uint64_t PrevData = PrevNElems*RecordSize + CRCSize*Vars.size();
      RH->Start = RH[-1].Start + PrevData;
      GH->NElems += RH->NElems;
    }

    // Compute the total file size.
    uint64_t LastNElems = RH[-1].NElems;
    uint64_t LastData = LastNElems*RecordSize + CRCSize*Vars.size();
    FileSize = RH[-1].Start + LastData;

    // Now that the starting offset has been computed, send it back to each rank.
    MPI_Scatter(&Header[GH->RanksStart], sizeof(RHLocal),
                MPI_BYTE, &RHLocal, sizeof(RHLocal),
                MPI_BYTE, 0, Comm);

    uint64_t HeaderCRC = crc64_omp(&Header[0], HeaderSize - CRCSize);
    crc64_invert(HeaderCRC, &Header[HeaderSize - CRCSize]);

    if (MPI_File_open(MPI_COMM_SELF, const_cast<char *>(FileName.c_str()),
                      MPI_MODE_WRONLY | MPI_MODE_CREATE,
                      MPI_INFO_NULL, &FH.get()) != MPI_SUCCESS)
      throw runtime_error("Unable to create the output file: " + FileName);

    MPI_File_set_size(FH.get(), FileSize);

    MPI_Status status;
    if (MPI_File_write_at(FH.get(), 0, &Header[0], HeaderSize, MPI_BYTE, &status) != MPI_SUCCESS) {
      close();
      throw runtime_error("Unable to write header to file: " + FileName);
    }

    close();
  } else {
    MPI_Gather(&RHLocal, sizeof(RHLocal), MPI_BYTE, 0, 0, MPI_BYTE, 0, Comm);
    MPI_Scatter(0, 0, MPI_BYTE, &RHLocal, sizeof(RHLocal), MPI_BYTE, 0, Comm);
  }

  MPI_Barrier(Comm);

  if (MPI_File_open(Comm, const_cast<char *>(FileName.c_str()), MPI_MODE_WRONLY,
      MPI_INFO_NULL, &FH.get()) != MPI_SUCCESS)
    throw runtime_error("Unable to create the output file: " + FileName);

  uint64_t Offset = RHLocal.Start;
  for (size_t i = 0; i < Vars.size(); ++i) {
    MPI_Status status;
    uint64_t WriteSize = NElems*Vars[i].Size;
    uint64_t CRC = crc64_omp(Vars[i].Data, WriteSize);
    char *CRCLoc = Vars[i].HasExtraSpace ?
      ((char *) Vars[i].Data) + WriteSize : (char *) &CRC;
    crc64_invert(CRC, CRCLoc);

    int e;
    if (Vars[i].HasExtraSpace) {
      e = MPI_File_write_at(FH.get(), Offset, Vars[i].Data, WriteSize + CRCSize,
                            MPI_BYTE, &status);
    } else {
      e =  MPI_File_write_at(FH.get(), Offset, Vars[i].Data, WriteSize,
                             MPI_BYTE, &status);
      if (e == MPI_SUCCESS)
        e = MPI_File_write_at(FH.get(), Offset + WriteSize, CRCLoc, CRCSize,
                              MPI_BYTE, &status);
    }

    if (e != MPI_SUCCESS) {
      close();
      throw runtime_error("Unable to write variable: " + Vars[i].Name +
                              " to file: " + FileName);
    }
    Offset += WriteSize + CRCSize;
  }

  close();
  MPI_Barrier(Comm);

  double EndTime = MPI_Wtime();
  double TotalTime = EndTime - StartTime;
  double MaxTotalTime;
  MPI_Reduce(&TotalTime, &MaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);

  if (Rank == 0) {
    double Rate = ((double) FileSize) / MaxTotalTime / (1024.*1024.);
    cout << "Wrote " << Vars.size() << " variables to " << FileName <<
            " (" << FileSize << " bytes) in " << MaxTotalTime << "s: " <<
            Rate << " MB/s" << endl;
  }
}
#endif // GENERICIO_NO_MPI

// Note: Errors from this function should be recoverable. This means that if
// one rank throws an exception, then all ranks should.
void GenericIO::openAndReadHeader(bool MustMatch) {
  const char *Magic = isBigEndian() ? MagicBE : MagicLE;
  const char *MagicInv = isBigEndian() ? MagicLE : MagicBE;

  int NRanks, Rank;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(Comm, &Rank);
  MPI_Comm_size(Comm, &NRanks);
#else
  Rank = 0;
  NRanks = 1;
#endif

  uint64_t HeaderSize;
  vector<char> Header;

  if (Rank == 0) {
#ifndef GENERICIO_NO_MPI
    char True = 1, False = 0;
    if (MPI_File_open(MPI_COMM_SELF, const_cast<char *>(FileName.c_str()),
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &FH.get()) != MPI_SUCCESS)
#else
    FH.get().open(FileName.c_str(), ios_base::in);
    if (!FH.get())
#endif
    {
#ifndef GENERICIO_NO_MPI
      MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
#endif
      throw runtime_error("Unable to open: " + FileName);
    }

    GlobalHeader GH;
#ifndef GENERICIO_NO_MPI
    MPI_Status status;
    if (MPI_File_read_at(FH.get(), 0, &GH, sizeof(GlobalHeader),
        MPI_BYTE, &status) != MPI_SUCCESS)
#else
    if (!FH.get().seekg(0) || !FH.get().read((char *) &GH, sizeof(GlobalHeader)))
#endif
    {
#ifndef GENERICIO_NO_MPI
      MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
#endif
      close();
      throw runtime_error("Unable to read global header: " + FileName);
    }

    if (string(GH.Magic, GH.Magic + MagicSize - 1) != Magic)
      {
      string Error;
      if (string(GH.Magic, GH.Magic + MagicSize - 1) == MagicInv)
        {
        this->ShouldSwap = true;
        for(int i=0; i < MagicSize; ++i )
          {
          this->SwapEndian(&GH.Magic[i],sizeof(char));
          }
        this->SwapEndian(&GH.HeaderSize,sizeof(uint64_t));
        this->SwapEndian(&GH.NElems,sizeof(uint64_t));
        this->SwapEndian(&GH.Dims[0],sizeof(uint64_t));
        this->SwapEndian(&GH.Dims[1],sizeof(uint64_t));
        this->SwapEndian(&GH.Dims[2],sizeof(uint64_t));
        this->SwapEndian(&GH.NVars,sizeof(uint64_t));
        this->SwapEndian(&GH.VarsSize,sizeof(uint64_t));
        this->SwapEndian(&GH.VarsStart,sizeof(uint64_t));
        this->SwapEndian(&GH.NRanks,sizeof(uint64_t));
        this->SwapEndian(&GH.RanksSize,sizeof(uint64_t));
        this->SwapEndian(&GH.RanksStart,sizeof(uint64_t));
        this->SwapEndian(&GH.GlobalHeaderSize,sizeof(uint64_t));
        }
      else
        {
        Error = "invalid file-type identifier";
        }

#ifndef GENERICIO_NO_MPI
      if( !this->ShouldSwap )
        {
        MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
        close();
        throw runtime_error("Won't read " + FileName + ": " + Error);
        }
#endif
      }

    if (MustMatch) {
      if (NRanks != (int) GH.NRanks) {
        stringstream ss;
        ss << "Won't read " << FileName << ": communicator-size mismatch: " <<
              "current: " << NRanks << ", file: " << GH.NRanks;
#ifndef GENERICIO_NO_MPI
        MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
#endif
        close();
        throw runtime_error(ss.str());
      }

#ifndef GENERICIO_NO_MPI
      int TopoStatus;
      MPI_Topo_test(Comm, &TopoStatus);
      if (TopoStatus == MPI_CART) {
        int Dims[3], Periods[3], Coords[3];
        MPI_Cart_get(Comm, 3, Dims, Periods, Coords);

        bool DimsMatch = true;
        for (int i = 0; i < 3; ++i) {
          if ((uint64_t) Dims[i] != GH.Dims[i]) {
            DimsMatch = false;
            break;
          }
        }

        if (!DimsMatch) {
          stringstream ss;
          ss << "Won't read " << FileName <<
                ": communicator-decomposition mismatch: " <<
                "current: " << Dims[0] << "x" << Dims[1] << "x" << Dims[2] <<
                ", file: " << GH.Dims[0] << "x" << GH.Dims[1] << "x" <<
                GH.Dims[2];
          MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);

          close();
          throw runtime_error(ss.str());
        }
      }
#endif
    }

    HeaderSize = GH.HeaderSize;
    Header.resize(HeaderSize + CRCSize, 0xFE /* poison */);
#ifndef GENERICIO_NO_MPI
    if (MPI_File_read_at(FH.get(), 0, &Header[0], HeaderSize + CRCSize,
        MPI_BYTE, &status) != MPI_SUCCESS)
#else
    if (!FH.get().seekg(0) || !FH.get().read((char *) &Header[0], HeaderSize + CRCSize))
#endif
    {
#ifndef GENERICIO_NO_MPI
      MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
#endif
      close();
      throw runtime_error("Unable to read header: " + FileName);
    }

    uint64_t CRC = crc64_omp(&Header[0], HeaderSize + CRCSize);
    if (CRC != (uint64_t) -1) {
#ifndef GENERICIO_NO_MPI
      MPI_Bcast(&False, 1, MPI_BYTE, 0, Comm);
#endif
      close();
      throw runtime_error("Header CRC check failed: " + FileName);
    }

#ifndef GENERICIO_NO_MPI
    close();
    MPI_Bcast(&True, 1, MPI_BYTE, 0, Comm);
#endif
  } else {
#ifndef GENERICIO_NO_MPI
    char Okay;
    MPI_Bcast(&Okay, 1, MPI_BYTE, 0, Comm);
    if (!Okay)
      throw runtime_error("Failure broadcast from rank 0");
#endif
  }

#ifndef GENERICIO_NO_MPI
  MPI_Bcast(&HeaderSize, 1, MPI_UINT64_T, 0, Comm);
#endif

  Header.resize(HeaderSize, 0xFD /* poison */);
#ifndef GENERICIO_NO_MPI
  MPI_Bcast(&Header[0], HeaderSize, MPI_BYTE, 0, Comm);
#endif

  FH.getHeaderCache().clear();
  FH.getHeaderCache().swap(Header);

#ifndef GENERICIO_NO_MPI
  MPI_Barrier(Comm);
  if (MPI_File_open(Comm, const_cast<char *>(FileName.c_str()),
                    MPI_MODE_RDONLY, MPI_INFO_NULL, &FH.get()) != MPI_SUCCESS)
    throw runtime_error("Unable to open: " + FileName);
#endif
}

void GenericIO::SwapEndian(void *Addr, const int Nb)
{
  assert("pre: pointer address is NULL!" && (Addr != NULL) );
  assert("NB > 0" && Nb > 0);

  // STEP 1: Allocate buffer where the data will be swapped
  char *Swapped = new char[Nb];
  assert("pre: Cannot allocate swapped buffer!" && (Swapped != NULL) );

  // STEP 2: Swap data
  for(int srcOffSet=Nb-1, idx=0; srcOffSet >= 0; --srcOffSet,++idx)
    {
    Swapped[idx] = *( (char*)Addr+srcOffSet);
    } // END for all bytes

  // STEP 3: Copy Swapped data to input buffer
  memcpy(Addr,(void*)Swapped,Nb);

  // STEP 4:  Clean dynamically allocated memory
  delete [] Swapped;
}

int GenericIO::readNRanks() {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  return (int) GH->NRanks;
}

void GenericIO::readDims(int Dims[3]) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  std::copy(GH->Dims, GH->Dims + 3, Dims);
}

uint64_t GenericIO::readTotalNumElems() {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  return GH->NElems;
}

size_t GenericIO::readNumElems(int EffRank) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  if (EffRank == -1) {
#ifndef GENERICIO_NO_MPI
    MPI_Comm_rank(Comm, &EffRank);
#else
    EffRank = 0;
#endif
  }

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  assert(EffRank <= (int) GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               EffRank*GH->RanksSize];
  return (size_t) RH->NElems;
}

void GenericIO::readCoords(int Coords[3], int EffRank) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  if (EffRank == -1) {
#ifndef GENERICIO_NO_MPI
    MPI_Comm_rank(Comm, &EffRank);
#else
    EffRank = 0;
#endif
  }

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  assert(EffRank <= (int) GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               EffRank*GH->RanksSize];

  std::copy(RH->Coords, RH->Coords + 3, Coords);
}

// Note: Errors from this function should be recoverable. This means that if
// one rank throws an exception, then all ranks should.
void GenericIO::readData(int EffRank, bool PrintStats, bool CollStats) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  int Rank;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(Comm, &Rank);
#else
  Rank = 0;
#endif

  if (EffRank == -1)
    EffRank = Rank;

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  assert(EffRank <= (int) GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               EffRank*GH->RanksSize];

  uint64_t TotalReadSize = 0;
#ifndef GENERICIO_NO_MPI
  double StartTime = MPI_Wtime();
#else
  double StartTime = double(clock())/CLOCKS_PER_SEC;
#endif

  int NErrs[2] = { 0, 0 };
  for (size_t i = 0; i < Vars.size(); ++i) {
    uint64_t Offset = RH->Start;
    bool VarFound = false;
    for (uint64_t j = 0; j < GH->NVars; ++j) {
      VariableHeader *VH = (VariableHeader *) &FH.getHeaderCache()[GH->VarsStart +
                                                           j*GH->VarsSize];

      string VName(VH->Name, VH->Name + NameSize);
      size_t VNameNull = VName.find('\0');
      if (VNameNull < NameSize)
        VName.resize(VNameNull);

      uint64_t ReadSize = RH->NElems*VH->Size + CRCSize;
      if (VName != Vars[i].Name) {
        Offset += ReadSize;
        continue;
      }

      VarFound = true;
      bool IsFloat = (bool) (VH->Flags & FloatValue),
           IsSigned = (bool) (VH->Flags & SignedValue);
      if (VH->Size != Vars[i].Size) {
        stringstream ss;
        ss << "Size mismatch for variable " << Vars[i].Name <<
              " in: " << FileName << ": current: " << Vars[i].Size <<
              ", file: " << VH->Size;
        throw runtime_error(ss.str());
      } else if (IsFloat != Vars[i].IsFloat) {
        string Float("float"), Int("integer");
        stringstream ss;
        ss << "Type mismatch for variable " << Vars[i].Name <<
              " in: " << FileName << ": current: " <<
              (Vars[i].IsFloat ? Float : Int) <<
              ", file: " << (IsFloat ? Float : Int);
        throw runtime_error(ss.str());
      } else if (IsSigned != Vars[i].IsSigned) {
        string Signed("signed"), Uns("unsigned");
        stringstream ss;
        ss << "Type mismatch for variable " << Vars[i].Name <<
              " in: " << FileName << ": current: " <<
              (Vars[i].IsSigned ? Signed : Uns) <<
              ", file: " << (IsSigned ? Signed : Uns);
        throw runtime_error(ss.str());
      }

      assert(Vars[i].HasExtraSpace && "Extra space required for reading");

#ifndef GENERICIO_NO_MPI
      MPI_Status status;
      if (MPI_File_read_at(FH.get(), Offset, Vars[i].Data, ReadSize,
                           MPI_BYTE, &status) != MPI_SUCCESS)
#else

      if (!FH.get().seekg(Offset) || !FH.get().read((char *) Vars[i].Data, ReadSize))
#endif
      {
        ++NErrs[0];
        break;
      }

      TotalReadSize += ReadSize;

      uint64_t CRC = crc64_omp(Vars[i].Data, ReadSize);
      if (CRC != (uint64_t) -1) {
        ++NErrs[1];
      }

      break;
    }

    if (!VarFound)
      throw runtime_error("Variable " + Vars[i].Name +
                          " not found in: " + FileName);

    if (NErrs[0] || NErrs[1])
      break;

    // Swap data elements
    if( this->ShouldSwap )
      {
      for(int idx=0; idx < this->NElems; ++idx )
        {
        size_t offSet = idx*Vars[i].Size;
        char* dataPtr = reinterpret_cast<char *>(Vars[i].Data)+offSet;
        this->SwapEndian(dataPtr,Vars[i].Size);
        }
      } // END if
  }

  int AllNErrs[2];
#ifndef GENERICIO_NO_MPI
  if(CollStats)
    {
    MPI_Allreduce(&NErrs, &AllNErrs, 2, MPI_INT, MPI_SUM, Comm);
    }
  else
    {
    AllNErrs[0] = NErrs[0]; AllNErrs[1] = NErrs[1];
    }
#else
  AllNErrs[0] = NErrs[0]; AllNErrs[1] = NErrs[1];
#endif

  if (AllNErrs[0] > 0 || AllNErrs[1] > 0) {
    stringstream ss;
    ss << "Experienced " << AllNErrs[0] << " I/O error(s) and " <<
          AllNErrs[1] << " CRC error(s) reading: " << FileName;
    throw runtime_error(ss.str());
  }

#ifndef GENERICIO_NO_MPI
  double EndTime = MPI_Wtime();
#else
  double EndTime = double(clock())/CLOCKS_PER_SEC;
#endif

  double TotalTime = EndTime - StartTime;
  double MaxTotalTime;
#ifndef GENERICIO_NO_MPI
  if(CollStats)
    {
    MPI_Reduce(&TotalTime, &MaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);
    }
#else
  MaxTotalTime = TotalTime;
#endif

  uint64_t AllTotalReadSize;
#ifndef GENERICIO_NO_MPI
  if(CollStats)
    {
    MPI_Reduce(&TotalReadSize, &AllTotalReadSize, 1, MPI_UINT64_T, MPI_SUM, 0, Comm);
    }
#else
  AllTotalReadSize = TotalReadSize;
#endif

  if (Rank == 0 && PrintStats) {
    double Rate = ((double) AllTotalReadSize) / MaxTotalTime / (1024.*1024.);
    cout << "Read " << Vars.size() << " variables from " << FileName <<
            " (" << AllTotalReadSize << " bytes) in " << MaxTotalTime << "s: " <<
            Rate << " MB/s [excluding header read]" << endl;
  }

}

void GenericIO::getVariableInfo(vector<VariableInfo> &VI) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  for (uint64_t j = 0; j < GH->NVars; ++j) {
    VariableHeader *VH = (VariableHeader *) &FH.getHeaderCache()[GH->VarsStart +
                                                         j*GH->VarsSize];

    string VName(VH->Name, VH->Name + NameSize);
    size_t VNameNull = VName.find('\0');
    if (VNameNull < NameSize)
      VName.resize(VNameNull);

    bool IsFloat = (bool) (VH->Flags & FloatValue),
         IsSigned = (bool) (VH->Flags & SignedValue);
    VI.push_back(VariableInfo(VName, (size_t) VH->Size, IsFloat, IsSigned));
  }
}

}
