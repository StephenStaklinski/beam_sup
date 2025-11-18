#ifndef SIMULATION_H
#define SIMULATION_H

#include "cell.h"
#include "tree.h"
#include <memory>
#include <unordered_set>
#include <map>
#include <lemon/bfs.h>

class Simulation
{
public:
  // Constructor
  Simulation(double K,
             std::vector<double> migrationRates,
             double mutationRate,
             double driverProb,
             double mutFreqThreshold,
             int maxNrAnatomicalSites,
             int maxGenerations,
             int downsampleCellNumber,
             int inputNumPossibleAnatomicalSites,
             std::vector<std::vector<double>> migrationTransitionProbs,
             int migrationStartGeneration,
             int migrationEndGeneration,
             bool resolvePolytomies);
  
  /// Destructor
  ~Simulation();
  
  /// Simulate a cell tree. Return true if simulation was successful and false otherwise
  bool simulate();

  /// Return simulated clone tree
  const Tree& getCellTree() const
  {
    assert(_pCellT);
    return *_pCellT;
  }

  /// Return vertex labeling
  const StringNodeMap& getCellVertexLabeling() const
  {
    assert(_cellVertexLabeling);
    return *_cellVertexLabeling;
  }

  typedef std::map<int, int> IntIntMap;
  
  /// Return cell to generation number map
  const IntIntMap getCellToGenerationNumber() const
  {
    return _cellToGenerationNumber;
  }

  /// Get newick string for a digraph tree with branch lengths as the number of total generations (same as cell divisions in this setup) for the mutations at the node branch point
  void writeNewick(std::ostream& os, const Tree& pCellT, const IntIntMap& mutationToGenerationNumber) const;

  /// Write migration history implied by the cell tree to output stream
  void writeMigrationHistory(std::ostream& os, const Tree& pCellT, const StringNodeMap& cellAnatomicalSites) const;
  
  /// Build cell Tree from final living cells and total allCells map
  void constructCellTree();

  /// Resolve polytomies in the cell tree
  void resolvePolytomiesInTree(Digraph& cellTree, Node root, StringNodeMap& cellLabels, StringNodeMap& cellAnatomicalSites);
  
private:
  /// Map from set of driver mutations to carrying capacity
  typedef std::map<IntSet, double> AnatomicalSiteMap;
  /// Anatomical site and drivers-specific carrying capacities
  typedef std::map<int, AnatomicalSiteMap> AnatomicalSiteFactorMap;
  /// Vector of cells
  typedef std::vector<Cell> CellVector;
  /// Map of set of driver mutations to cell vector (with those drivers)
  typedef std::map< int, std::map<IntSet, CellVector>> ClonalComposition;

  /// Return a new mutation
  int getNewMutation()
  {
    return _nrMutations++;
  }

  /// Initialize simulation
  void init();
  
  /// Update anatomical site factors
  void updateAnatomicalSiteFactors();
  
  /// Perform potential migrations
  void migrate();

  /// Return target anatomical site
  ///
  /// @param s Source anatomical site
  int getTargetAnatomicalSite(int s);

  /// Return the number of cells in the provided anatomical site.
  /// If weighted is true, weight by the number of driver mutations they possess.
  ///
  /// @param s Anatomical site
  /// @param weighted Whether to weight by number of driver mutations
  double getNrExtantCells(int s, bool weighted = false) const
  {
    double res = 0;
    for (const auto& kv : _extantCellsByDrivers.at(s))
    {
      if (weighted)
        res += kv.first.size() * kv.second.size();
      else
        res += kv.second.size();
    }
    return res;
  }
  
  /// Draw cells that will migrate
  ///
  /// @param s Anatomical site
  /// @param migrationCount Number of cells to draw
  IntVector draw(int s, int migrationCount) {
    // Get cell group sizes for anatomical site s
    const auto& groups = _extantCellsByDrivers.at(s);
    IntVector groupSizes;
    for (const auto& kv : groups) {
      groupSizes.push_back(kv.second.size());
    }

    // Compute total number of cells
    int totalCells = std::accumulate(groupSizes.begin(), groupSizes.end(), 0);

    // Compute cumulative proportions for each group
    DoubleVector cumProportions;
    double cumSum = 0;
    for (int sz : groupSizes)
    {
      cumSum += static_cast<double>(sz) / totalCells;
      cumProportions.push_back(cumSum);
    }

    // Initialize result vector to track drawn cells per group
    IntVector result(groupSizes.size(), 0);
    std::uniform_real_distribution<> unif(0, 1);

    // Draw cells according to group proportions, not exceeding group sizes
    for (int drawn = 0; drawn < migrationCount;)
    {
      double r = unif(g_rng);
      // Find group index for random value
      auto it = std::lower_bound(cumProportions.begin(), cumProportions.end(), r);
      int idx = std::distance(cumProportions.begin(), it);
      // Only draw if group not exhausted
      if (result[idx] < groupSizes[idx])
      {
        ++result[idx];
        ++drawn;
      }
    }
    return result;
  }
  
  /// Maximum number of cells
  const int _maxGenerations;
  /// Generation
  int _generation;
  /// Extant cells per anatomical site further split by drivers
  ClonalComposition _extantCellsByDrivers;
  /// All cells born in the simulation process
  CellVector _allCells;
  /// Number of extant cells
  int _nrExtantCells;
  /// Number of mutations
  int _nrMutations;
  /// Number of active anatomical site
  int _nrActiveAnatomicalSites;
  /// Is active anatomical site
  std::map<int, bool> _isActiveAnatomicalSite;
  /// Anatomical site factors
  AnatomicalSiteFactorMap _anatomicalSiteFactors;
  /// Anatomical site label
  std::map<int, std::string> _anatomicalSiteLabel;
  /// Carrying capacity for each site
  const double _K;
  /// Migration rate
  std::vector<double> _migrationRates;
  /// Mutation rate
  const double _mutationRate;
  /// Mutation frequency threshold
  const double _mutFreqThreshold;
  /// Driver mutation probability
  const double _driverProb;
  /// Maximum number of anatomical sites
  const int _maxNrAnatomicalSites;
  /// Set of driver mutations
  IntSet _driverMutations;
  /// Cell-tree mutation to observable mutation
  IntIntMap _cellTreeMutationToObservableMutation;
  /// Cell tree
  Tree* _pCellT;
  /// Sampled leaves
  NodeSet _sampledLeaves;
  /// Sampled mutations
  IntSet observableMutations;
  /// Cell vertex labeling
  StringNodeMap* _cellVertexLabeling;
  /// Map of cell ID to generation time
  IntIntMap _cellToGenerationNumber;
  /// Iterator for all cell IDs
  int _idIterator;
  /// Map for cell IDs to cell object
  std::map<int, Cell> _idCellMap;
  /// Desired number of cells to downsample until
  int _desiredNumberCells;
  /// Transition probability matrix for tissue migrations (must be a square matrix!)
  std::vector<std::vector<double>> _migrationTransitionProbs;
  /// Vector to track existing tissue sites
  IntSet _existingSites;
  /// Input number of possible anatomical sites (for default migration matrix when one is not provided as input)
  int _inputNumPossibleAnatomicalSites;
  /// Generation to start migrations
  int _migrationStartGeneration;
  /// Generation to stop migrations
  int _migrationEndGeneration;
  /// Whether to resolve polytomies in the cell tree
  bool _resolvePolytomies;
  };

#endif