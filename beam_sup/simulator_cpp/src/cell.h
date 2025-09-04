#ifndef CELL_H
#define CELL_H

#include "utils.h"

/// This class models a cell, the central unit in the agent-based model
class Cell
{
public:
  /// Constructs a cell with specified passenger mutations, mutation, anatomical site, ID, parent ID, and generation
  Cell(const IntVector& passengerMutations,
       int mutation,
       int anatomicalSite,
       int id,
       int parent,
       int generation);

  /// Constructs a default cell with uninitialized members
  Cell();

  /// Possible outcomes of a cell cycle
  enum Outcome
  {
    /// Cell replicates
    REPLICATION,
    /// Cell dies
    DEATH
  };

  /// Simulates a cell cycle generation and returns the outcome
  /// @param logisticFactor Factor controlling replication probability
  /// @param nrDriverMutations Number of driver mutations in the cell
  Outcome performGeneration(double logisticFactor,
                            int nrDriverMutations) const;

  /// Migrates the cell to a new anatomical site
  /// @param newAnatomicalSite The anatomical site to migrate to
  void migrate(int newAnatomicalSite)
  {
    _anatomicalSite = newAnatomicalSite;
  }

  /// Returns the anatomical site of the cell
  int getAnatomicalSite() const
  {
    return _anatomicalSite;
  }

  /// Returns the set of passenger mutations in the cell
  const IntVector& getPassengerMutations() const
  {
    return _passengerMutations;
  }

  /// Returns the mutation introduced by this cell
  int getMutation() const
  {
    return _mutation;
  }

  /// Returns the unique ID of the cell
  int getID() const
  {
    return _id;
  }

  /// Returns the parent cell's ID
  int getParentID() const
  {
    return _parent;
  }

  /// Returns the generation number of the cell
  int getGeneration() const
  {
    return _generation;
  }

protected:
  /// Set of passenger mutations carried by the cell
  IntVector _passengerMutations;
  /// Mutation introduced by this cell
  int _mutation;
  /// Anatomical site where the cell resides
  int _anatomicalSite;
  /// Unique identifier for the cell
  int _id;
  /// Identifier of the parent cell
  int _parent;
  /// Generation number of the cell
  int _generation;
};

#endif // CELL_H
