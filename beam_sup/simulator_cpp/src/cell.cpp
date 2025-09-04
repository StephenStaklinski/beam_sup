#include "cell.h"

// Constructor initializing all member variables with provided values
Cell::Cell(const IntVector& passengerMutations,
  int mutation,
  int anatomicalSite,
  int id,
  int parent,
  int generation)
  : _passengerMutations(passengerMutations)
  , _mutation(mutation)
  , _anatomicalSite(anatomicalSite)
  , _id(id)
  , _parent(parent)
  , _generation(generation)
{
}

// Default constructor initializing member variables with default values
Cell::Cell()
  : _passengerMutations()
  , _mutation(-1)
  , _anatomicalSite(-1)
  , _id(-1)
  , _parent(-1)
  , _generation(-1)
{
}

// Simulates one generation for the cell: possible replication (+mutation) or death
Cell::Outcome Cell::performGeneration(double logisticFactor,
          int nrDriverMutations) const
{

  // Uniform random number generator
  std::uniform_real_distribution<> unif(0, 1);

  // Selection coefficient based on logistic factor
  const double s_0 = 0.1 * logisticFactor;

  // Initial birth rate
  double birthRate = 0.5 * (1 + s_0);

  // Increase birth rate for each driver mutation
  for (int i = 0; i < nrDriverMutations; ++i)
  {
    birthRate *= (1 + s_0);
  }

  // Clamp birthRate to [0, 1]
  if (birthRate < 0)
    birthRate = 0;
  else if (birthRate > 1)
    birthRate = 1;

  // Generate random number in [0, 1]
  double r = unif(g_rng);

  // Cell replicates
  if (r < birthRate)
  {
    return REPLICATION;
  }
  // Cell dies
  else
  {
    return DEATH;
  }
}
