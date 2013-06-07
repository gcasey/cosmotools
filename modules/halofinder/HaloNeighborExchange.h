/**
 * @brief HaloNeighborExchange implements functionality for exchanging a set
 * of halos with neighboring processes over the DIY communication
 * infrastructure.
 * @note To use this class, DIY must be initialized and decomposed accordingly.
 */
#ifndef HALONEIGHBOREXCHANGE_H_
#define HALONEIGHBOREXCHANGE_H_

#include "CosmologyToolsMacros.h"
#include "Halo.h" // For Halo

// MPI
#include <mpi.h>  // For MPI_Comm

// C/C++ includes
#include <map>    // For STL map
#include <vector> // For STL vector

// Data-structure to store halos based on hash-code
typedef std::map<std::string,cosmotk::Halo> HaloHashMap;

namespace cosmotk
{

class HaloNeighborExchange
{
public:
  HaloNeighborExchange();
  virtual ~HaloNeighborExchange();

  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Enqueues a halo to send to neighboring processes
   * @param halo pointer to the halo to send
   * @pre halo != NULL
   * @note No communication is done at this point. Typically, once
   * the application enqueues the halos to send, ExchangeHalos is
   * called to carry out the communication task.
   */
  void EnqueueHalo(Halo *halo);

  /**
   * @brief Exchanges the halos with the neighboring processes.
   * @param localHalos array consisting of the halos in this process.
   * @param N the number of local halos in this process.
   * @param exchangedHalos output vector of halos
   * @param inclusive if true the output vector will consist of both local
   * and neighboring halos. If not specified, default is set to true.
   * @note This method is collective, all ranks must call it.
   */
  void ExchangeHalos(
      Halo *localHalos, const int N,
      std::vector<Halo>& exchangedHalos,
      bool inclusive=true);

  /**
   * @brief Exchanges pre-enqueued halos in this process with all neighbors.
   * @param exchangedHalos output vector of halos.
   * @param inclusive if true the output vector will consist of both local
   * and neighboring halos. If not specified, default is set to true.
   * @note This method is collective, all ranks must call it.
   */
  void ExchangeHalos(
      std::vector<Halo>& exchangedHalos,
      bool inclusive=true);

protected:
  MPI_Comm Communicator;

  std::vector< Halo  > EnqueuedHalos;

  /**
   * @brief Exchanges the HaloInformation object of each halo.
   * @param localHalos array consisting of the halos in this process.
   * @param N the number of local halos in this process.
   * @param neighborHalos data-structure wherein neighboring halos are stored.
   * @see Halo::GetHaloInformation
   * @see HaloNeighborExchange::ExchangeHalos
   */
  void ExchangeHaloInformation(
      Halo *localHalos, const int N, HaloHashMap& neighborHalos);

  /**
   * @brief Exchanges the particles IDs of the halos.
   * @param localHalos array consisting of the halos in this process.
   * @param N the number of local halos in this process.
   * @param neighborHalos data-structure wherein neighboring halos are stored.
   * @see Halo::GetHaloInformation
   * @see HaloNeighborExchange::ExchangeHalos
   */
  void ExchangeHaloParticles(
      Halo *localHalos, const int N, HaloHashMap& neighborHalos);

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloNeighborExchange);
};

} /* namespace cosmotk */
#endif /* HALONEIGHBOREXCHANGE_H_ */
