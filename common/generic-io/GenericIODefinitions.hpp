/**
 * @brief Definitions of internal datal-structures and globals used by the
 * GenericIO framework.
 */
#ifndef GENERICIODEFINITIONS_HPP_
#define GENERICIODEFINITIONS_HPP_

namespace cosmotk {

static const size_t CRCSize   = 8;
static const size_t MagicSize = 8;
static const char *MagicBE    = "HACC01B";
static const char *MagicLE    = "HACC01L";
static const size_t NameSize  = 256;

struct VariableHeader {
  char Name[NameSize];
  uint64_t Flags;
  uint64_t Size;
} __attribute__((packed));

struct RankHeader {
  uint64_t Coords[3];
  uint64_t NElems;
  uint64_t Start;
  uint64_t GlobalRank;
} __attribute__((packed));

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
  double   PhysOrigin[3];
  double   PhysScale[3];
} __attribute__((packed));

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

}


#endif /* GENERICIODEFINITIONS_HPP_ */
