#include "simulation.h"
#include "lemon/bfs.h"

Simulation::Simulation(double K,
                       std::vector<double> migrationRates,
                       double mutationRate,
                       double driverProb,
                       double mutFreqThreshold,
                       int maxNrAnatomicalSites,
                       int maxGenerations,
                       int downsampleCellNumber,
                       std::vector<std::vector<double>> migrationTransitionProbs)
  : _G()
  , _rootG(lemon::INVALID)
  , _anatomicalSiteMap(_G)
  , _indexToVertexG()
  , _generation(0)
  , _extantCellsByDrivers()
  , _nrExtantCells(0)
  , _nrMutations(0)
  , _nrActiveAnatomicalSites(0)
  , _isActiveAnatomicalSite()
  , _anatomicalSiteFactors()
  , _K(K)
  , _migrationRates(migrationRates)
  , _mutationRate(mutationRate)
  , _driverProb(driverProb)
  , _mutFreqThreshold(mutFreqThreshold)
  , _maxNrAnatomicalSites(maxNrAnatomicalSites)
  , _maxGenerations(maxGenerations)
  , _driverMutations()
  , _cellTreeMutationToObservableMutation()
  , _pCellT(NULL)
  , _cellVertexLabeling(NULL)
  ,  _idIterator()
  , _idCellMap()
  , _cellToGenerationNumber()
  , _desiredNumberCells(downsampleCellNumber)
  , _migrationTransitionProbs(migrationTransitionProbs)
  , _existingSites()
  , _anatomicalSiteLabel()
{}

Simulation::~Simulation()
{
  delete _pCellT;
  delete _cellVertexLabeling;
}

// Initialize the simulation state for a new run
void Simulation::init()
{
  // Set up the founder cell with no passenger mutations and mutation ID 0
  IntVector initialPassengerMutations; // Founder cell starts with no passenger mutations
  int initialMutation = 0; // Founder cell's mutation ID
  _driverMutations.clear(); // Clear any previous driver mutations
  _driverMutations.insert(initialMutation); // Insert founder mutation as the first driver mutation
  _nrMutations = 1; // Initialize mutation counter

  // Initialize anatomical site 0 as the starting site
  int initialAnatomicalSite = 0;
  _existingSites.clear(); // Remove any previous anatomical sites
  _existingSites.insert(initialAnatomicalSite); // Add site 0 as the initial primary site
  _isActiveAnatomicalSite[initialAnatomicalSite] = true; // Mark site 0 as active
  _nrActiveAnatomicalSites = 1; // Only one active site at start
  _generation = 0; // Start at generation 0
  _idIterator = 1; // Cell ID counter starts at 1

  int parent = 0; // Founder cell has no parent (ID 0)

  // Create the founder cell and register it in all relevant data structures
  Cell founder(initialPassengerMutations,
               initialMutation,
               initialAnatomicalSite,
               _idIterator,
               parent,
               _generation);

  _idCellMap[_idIterator] = founder; // Map cell ID to founder cell
  _cellToGenerationNumber[_idIterator] = _generation; // Map cell ID to generation number
  _extantCellsByDrivers[initialAnatomicalSite][_driverMutations].push_back(founder); // Add founder to extant cells by driver group

  // Initialize migration graph for anatomical sites
  _G.clear(); // Remove any previous graph structure
  _rootG = _G.addNode(); // Add root node for site 0
  _anatomicalSiteMap[_rootG] = initialAnatomicalSite; // Map root node to site 0
  _indexToVertexG.clear(); // Clear previous site-to-node mapping
  _indexToVertexG[initialAnatomicalSite] = _rootG; // Map site 0 to root node

  _nrExtantCells = 1; // Only the founder cell exists initially

  _generation++; // Advance to next generation (generation 1)

  // If migration transition probabilities are not provided, set up a default uniform matrix
  if (_migrationTransitionProbs.empty()) {
    int numPossibleAnatomicalSites = 10; // Default number of possible anatomical sites
    _migrationTransitionProbs.resize(numPossibleAnatomicalSites, std::vector<double>(numPossibleAnatomicalSites, 0.0));
    for (int i = 0; i < numPossibleAnatomicalSites; ++i) {
      double prob = 1.0 / (numPossibleAnatomicalSites - 1);
      for (int j = 0; j < numPossibleAnatomicalSites; ++j) {
        if (i != j) _migrationTransitionProbs[i][j] = prob;
      }
    }
  }

  int numPossibleAnatomicalSites = _migrationTransitionProbs.size(); // Number of sites from transition matrix

  if (_migrationRates.size() == 1) {
    // If only one rate is provided, replicate it for all sites
    _migrationRates.resize(numPossibleAnatomicalSites, _migrationRates[0]);
  } else if (_migrationRates.size() != numPossibleAnatomicalSites) {
    // If migration rates vector length does not match transition matrix, exit with error
    std::cout << "Error: Length of migrationRates input does not match the dimensions of the transition probability matrix." << std::endl;
    std::exit(1);
  }
}

// Updates anatomical site factors for all existing sites, determines active sites, and removes dead sites
void Simulation::updateAnatomicalSiteFactors()
{
  _anatomicalSiteFactors.clear();      // Reset anatomical site factors
  _nrActiveAnatomicalSites = 0;        // Reset count of active anatomical sites
  IntSet sitesToErase;                 // Set to track sites to be removed

  // Iterate over all existing anatomical sites
  for (auto s : _existingSites)
  {
    _isActiveAnatomicalSite[s] = false; // Mark site as inactive initially
    int nrCells_s = 0;                  // Total number of cells in site s

    // Remove empty driver mutation groups in site s
    auto& driverGroups = _extantCellsByDrivers.at(s);
    for (auto it = driverGroups.begin(); it != driverGroups.end(); ) {
      if (it->second.empty())
        it = driverGroups.erase(it);
      else
        ++it;
    }

    // Calculate anatomical site factors and count cells for each driver mutation group
    // For each driver mutation group X in anatomical site s
    for (const auto& kv : _extantCellsByDrivers.at(s)) {
      const IntSet& X = kv.first;           // Driver mutation set for this group
      const CellVector& cells = kv.second;  // Vector of extant cells in this group
      int nrCells_sX = cells.size();        // Number of cells in this driver group
      nrCells_s += nrCells_sX;              // Accumulate total cell count for site s

      // Calculate logistic growth factor for this driver group in site s
      // This factor determines the probability of cell replication vs. death
      // It is capped at 0 (no growth) if the group is at or above carrying capacity
      _anatomicalSiteFactors[s][X] = std::max(0.0, 1 - double(nrCells_sX) / (_K * X.size()));
    }

    // Mark site as active if it has more than 5000 cells (detectable threshold)
    if (nrCells_s > 5000) {
      ++_nrActiveAnatomicalSites;
      _isActiveAnatomicalSite[s] = true;
    } else if (nrCells_s == 0) {
      sitesToErase.insert(s); // Site has no cells, schedule it for removal
    }
  }

  // Remove dead anatomical sites and warn if migration events are lost
  for (auto s : sitesToErase) {
    _existingSites.erase(s);
    Node v = _indexToVertexG[s];
    if (lemon::countOutArcs(_G, v) > 0) {
      std::cout << "WARNING: Migration events with site " << s << " are lost since it died." << std::endl;
    }
    _G.erase(v);
  }
}

// Select a target anatomical site for migration from site s based on transition probabilities
int Simulation::getTargetAnatomicalSite(int s)
{
  // Prepare random number generator for discrete selection
  std::uniform_real_distribution<> unif(0, 1);

  // Get transition probabilities for site s
  const std::vector<double>& transitionProbs = _migrationTransitionProbs[s];

  // Use discrete_distribution to select a target site according to transition probabilities
  std::discrete_distribution<> dist(transitionProbs.begin(), transitionProbs.end());
  int t = dist(g_rng);

  // If the target anatomical site t does not exist yet, initialize it in the migration graph and related structures
  if (_existingSites.find(t) == _existingSites.end())
  {
    Node v_t = _G.addNode();           // Add a new node to the migration graph for site t
    _indexToVertexG[t] = v_t;          // Map anatomical site index t to the new graph node
    _anatomicalSiteMap[v_t] = t;       // Associate the new node with anatomical site t
    _isActiveAnatomicalSite[t] = false;// Mark site t as inactive initially
    _existingSites.insert(t);          // Add site t to the set of existing anatomical sites
  }

  return t;
}


// Handles migration of cells between anatomical sites based on migration rates and transition probabilities
void Simulation::migrate()
{
  std::uniform_real_distribution<> unif(0, 1);
  std::poisson_distribution<> poisson(1);

  // Iterate over all existing anatomical sites
  for (auto s : _existingSites)
  {
    // Skip sites with no possible migration targets based on input transition probabilities
    const auto& transitionProbs = _migrationTransitionProbs[s];
    if (std::all_of(transitionProbs.begin(), transitionProbs.end(), [](double prob) { return prob == 0; })) {
      continue;
    }

    // Probabilistically decide if migration occurs from this site using the weighted number of cells
    if (unif(g_rng) > getNrExtantCells(s, true) * _migrationRates[s]) {
      continue;
    }

    // Select migration target site
    int t = getTargetAnatomicalSite(s);
    if (t == s) continue; // Prevent self-migration
    Node v_s = _indexToVertexG[s];  // Get the migration graph node for anatomical site s
    Node v_t = _indexToVertexG[t];  // Get the migration graph node for target site t

    // Unweighted number of extant cells in site s
    int nrExtantCells_s = getNrExtantCells(s, false);
    // Determine number of cells to migrate (at least 1, up to available cells)
    int migrationCount = std::min(1 + poisson(g_rng), nrExtantCells_s);

    std::cout << "\tMigration of " << migrationCount << " cells from " << s << " to " << t << std::endl;

    // Draw cells to migrate, proportionally per driver mutation group
    IntVector migrationVector = draw(s, migrationCount); // migrationVector[idx] = number of cells to migrate from group idx
    int idx = 0;

    // Iterate over each driver mutation group in anatomical site s
    for (auto& kv : _extantCellsByDrivers.at(s))
    {
      const IntSet& X = kv.first; // Driver mutation set for this group
      CellVector& extantCellByDrivers_sX = kv.second; // Vector of extant cells in this group

      // Shuffle group if any cells are migrating from it, to randomize selection
      if (migrationVector[idx] > 0) {std::shuffle(extantCellByDrivers_sX.begin(), extantCellByDrivers_sX.end(), g_rng);}

      IntSet migratedClones; // Track mutation IDs of migrated clones to avoid duplicates
      int migrated = 0; // Count of successfully migrated cells from this group
      // Migrate up to migrationVector[idx] cells from this group
      for (size_t i = 0; i < extantCellByDrivers_sX.size() && migrated < migrationVector[idx]; )
      {
        Cell& cell = extantCellByDrivers_sX.back(); // Select cell from end of vector
        // Avoid migrating duplicate clones (by mutation ID)
        if (migratedClones.count(cell.getMutation()) == 0) {
          migratedClones.insert(cell.getMutation()); // Mark mutation as migrated
          cell.migrate(t); // Update cell's anatomical site to target t
          Arc a = _G.addArc(v_s, v_t); // Add migration event to migration graph
          _extantCellsByDrivers[t][X].push_back(cell); // Move cell to target site t, same driver group
          extantCellByDrivers_sX.pop_back(); // Remove cell from source site
          ++migrated;
        }
        else {
          // If duplicate, remove and continue (do not migrate)
          extantCellByDrivers_sX.pop_back();
        }
      }
      ++idx; // Move to next driver mutation group
    }
  }
}

// Simulate tumor growth and migration, returning true if a viable tumor exists
bool Simulation::simulate()
{
  std::uniform_real_distribution<> unif(0, 1);

  init();

  while (_nrExtantCells > 0)
  {

    // Logging
    std::cout << "Generation: " << _generation
          << ", Total cells: " << _nrExtantCells
          << ", Active anatomical sites: " << _nrActiveAnatomicalSites
          << std::endl;

    // Check generation and anatomical site limits to see if we should stop the simulation
    if (_maxGenerations != -1 && _generation >= _maxGenerations)
    {
      std::cout << "Maximum number of generations reached" << std::endl;
      break;
    }
    if (_maxNrAnatomicalSites != -1 && _nrActiveAnatomicalSites >= _maxNrAnatomicalSites)
    {
      std::cout << "Maximum number of active anatomical sites reached" << std::endl;
      break;
    }

    // Update anatomical site factors and determine active sites
    updateAnatomicalSiteFactors();

    // Prepare for next generation
    ClonalComposition newExtantCellsByDrivers;
    int newNrExtantCells = 0;

    // Iterate over all existing anatomical sites
    for (auto s : _existingSites)
    {
      // Create a new map for driver mutation groups in anatomical site s for the next generation
      std::map<IntSet, CellVector> newExtantCellsByDrivers_sD;
      newExtantCellsByDrivers[s] = newExtantCellsByDrivers_sD;
      auto& newExtantCellsByDrivers_s = newExtantCellsByDrivers.at(s);

      // Iterate over all driver mutation groups (clones) in anatomical site s
      for (auto& kv : _extantCellsByDrivers.at(s))
      {
        const IntSet& X = kv.first; // Driver mutation set for this group
        // Ensure anatomical site factor exists for this driver group
        assert(_anatomicalSiteFactors[s].count(X) == 1);
        double logisticFactor = _anatomicalSiteFactors[s][X]; // Growth factor for this clone

        // Prepare vector for new extant cells of this driver group in site s
        auto& newExtantCellsByDrivers_sX = newExtantCellsByDrivers_s[X];
        newExtantCellsByDrivers_sX.reserve(_K * X.size()); // Reserve space for efficiency

        // Process each cell in the current driver group
        for (const Cell& cell : kv.second)
        {
          // Simulate cell fate for this generation (replication or death)
          switch (cell.performGeneration(logisticFactor, X.size()))
          {
          case Cell::REPLICATION:
          {
            // Daughter 1 is identical to parent, add to new population
            newExtantCellsByDrivers_sX.push_back(cell);

            // Generate initial passenger mutations for daughter 2 (from mutation rate floor)
            int nrInitialPassengers = static_cast<int>(floor(_mutationRate));
            IntVector initialPassengerMutations;
            for (int i = 0; i < nrInitialPassengers; ++i) {
              initialPassengerMutations.push_back(getNewMutation());
            }

            // Draw a random number to determine mutation type for daughter 2
            double r = unif(g_rng);

            // Daughter 2 acquires a new driver mutation
            if (r < _driverProb * (X.size() + 1))
            {
              // Copy parent's passenger mutations and add new initial passengers
              IntVector passengerMutations2 = cell.getPassengerMutations();
              passengerMutations2.insert(passengerMutations2.end(),
                            initialPassengerMutations.begin(),
                            initialPassengerMutations.end());

              // Generate new driver mutation
              int new_mut = getNewMutation();

              // Create new driver mutation set for daughter 2
              IntSet driverMutations2 = X;
              driverMutations2.insert(new_mut);
              _driverMutations.insert(new_mut); // Track new driver mutation globally

              // Create new cell object for daughter 2
              _idIterator++;
              Cell daughterCell2(passengerMutations2, new_mut, cell.getAnatomicalSite(),
                        _idIterator, cell.getID(), _generation);

              // Add daughter cell to new population under new driver group
              newExtantCellsByDrivers_s[driverMutations2].push_back(daughterCell2);
              _idCellMap[_idIterator] = daughterCell2;
              _cellToGenerationNumber[_idIterator] = _generation;
            }
            // Daughter 2 acquires a new passenger mutation
            else if (r < _mutationRate - floor(_mutationRate))
            {
              // Generate new passenger mutation
              int new_mut = getNewMutation();

              // Copy parent's passenger mutations and add new initial passengers and new passenger
              IntVector passengerMutations2 = cell.getPassengerMutations();
              passengerMutations2.insert(passengerMutations2.end(),
                            initialPassengerMutations.begin(),
                            initialPassengerMutations.end());
              passengerMutations2.push_back(new_mut);

              // Create new cell object for daughter 2 (same driver group)
              _idIterator++;
              Cell daughterCell2(passengerMutations2, new_mut, cell.getAnatomicalSite(),
                        _idIterator, cell.getID(), _generation);

              // Add daughter cell to new population under same driver group
              newExtantCellsByDrivers_sX.push_back(daughterCell2);
              _idCellMap[_idIterator] = daughterCell2;
              _cellToGenerationNumber[_idIterator] = _generation;
            }
            // Daughter 2 acquires no new mutation (only initial passengers if any)
            else
            {
              // If there are initial passenger mutations, use the last one as mutation ID; otherwise, use parent's mutation
              int new_mut = initialPassengerMutations.empty() ?
                      cell.getMutation() : initialPassengerMutations.back();

              // Copy parent's passenger mutations and add new initial passengers
              IntVector passengerMutations2 = cell.getPassengerMutations();
              passengerMutations2.insert(passengerMutations2.end(),
                            initialPassengerMutations.begin(),
                            initialPassengerMutations.end());

              // Create new cell object for daughter 2
              _idIterator++;
              Cell daughterCell2(passengerMutations2, new_mut, cell.getAnatomicalSite(),
                        _idIterator, cell.getID(), _generation);

              // Add daughter cell to new population under same driver group
              newExtantCellsByDrivers_sX.push_back(daughterCell2);
              _idCellMap[_idIterator] = daughterCell2;
              _cellToGenerationNumber[_idIterator] = _generation;
            }

            // Two daughter cells are produced (replication event)
            newNrExtantCells += 2;
            break;
          }
          case Cell::DEATH:
            // Cell dies, do not add to new population (no action needed)
            break;
          }
        }
      }
    }

    _nrExtantCells = newNrExtantCells;
    std::swap(_extantCellsByDrivers, newExtantCellsByDrivers);

    // Handle cell migration between anatomical sites
    migrate();

    // Advance to next generation
    ++_generation;
  }

  updateAnatomicalSiteFactors();

  // If tumor exists, label anatomical sites and construct cell tree
  if (_nrExtantCells > 0)
  {
    char buf[32];
    for (auto s : _existingSites)
    {
      if (_isActiveAnatomicalSite.at(s))
      {
        if (s == 0)
          _anatomicalSiteLabel[s] = "P";
        else
        {
          snprintf(buf, sizeof(buf), "M%d", s);
          _anatomicalSiteLabel[s] = buf;
        }
      }
    }

    if (!_isActiveAnatomicalSite.at(0))
    {
      std::cerr << "No primary!" << std::endl;
      return false;
    }

    constructCellTree();
    return true;
  }
  // No tumor exists, return false
  return false;
}

// Returns a filtered migration graph containing only active anatomical sites and migration events between them
MigrationGraph Simulation::getMigrationGraph() const
{
  // Assign labels to nodes for active anatomical sites
  StringNodeMap label(_G);
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    int site = _anatomicalSiteMap[v];
    if (_isActiveAnatomicalSite.at(site)) {
      label[v] = _anatomicalSiteLabel.at(site);
    }
  }

  // Prepare node and arc filters: only keep nodes/arcs for active sites
  BoolNodeMap nodeFilter(_G, false);
  BoolArcMap arcFilter(_G, false);

  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    int site = _anatomicalSiteMap[v];
    nodeFilter[v] = _isActiveAnatomicalSite.at(site);
  }

  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    // Only keep arcs between active sites
    if (nodeFilter[u] && nodeFilter[v]) {
      arcFilter[a] = true;
    }
  }

  // Create a subgraph view with only active nodes/arcs
  SubDigraph subG(_G, nodeFilter, arcFilter);

  // Copy the filtered subgraph into a new graph object
  Digraph filteredGraph;
  Node filteredRoot;
  StringNodeMap filteredLabels(filteredGraph);

  lemon::digraphCopy(subG, filteredGraph)
    .node(_rootG, filteredRoot)
    .nodeMap(label, filteredLabels)
    .run();

  // Return the migration graph object
  return MigrationGraph(filteredGraph, filteredRoot, filteredLabels);
}


// Constructs the cell lineage tree for sampled extant cells, labels nodes, and updates global pointers
void Simulation::constructCellTree() {
  Digraph cellTree;
  std::map<int, Node> idNodeMap;

  // Recursively add parent nodes and arcs for a given cell ID
  std::function<void(int)> addChild = [&](int currentChildID) {
    const Cell& currentCell = _idCellMap.at(currentChildID);
    int parentID = currentCell.getParentID();
    // Recursively add parent if not already present
    if (idNodeMap.find(parentID) == idNodeMap.end()) {
      addChild(parentID);
    }
    Node parentNode = idNodeMap.at(parentID);
    Node childNode = cellTree.addNode();
    idNodeMap[currentChildID] = childNode;
    cellTree.addArc(parentNode, childNode);
  };

  // Gather all extant cells from active anatomical sites
  std::vector<Cell> allExtantCells;
  for (int s : _existingSites) {
    if (!_isActiveAnatomicalSite.at(s)) continue;
    for (const auto& pair : _extantCellsByDrivers.at(s)) {
      for (const Cell& cell : pair.second) {
        allExtantCells.push_back(cell);
      }
    }
  }

  // Downsample cells if requested
  std::vector<Cell> sampledCells;
  if (_desiredNumberCells == -1 || _desiredNumberCells >= _nrExtantCells) {
    sampledCells = allExtantCells;
  } else {
    std::shuffle(allExtantCells.begin(), allExtantCells.end(), g_rng);
    sampledCells.assign(allExtantCells.begin(), allExtantCells.begin() + _desiredNumberCells);
  }

  // Add root node (founder cell, assumed ID 1)
  Node root = cellTree.addNode();
  idNodeMap[1] = root;

  // Build tree by adding sampled cells and their ancestors
  for (const Cell& cell : sampledCells) {
    int id = cell.getID();
    if (idNodeMap.find(id) == idNodeMap.end()) {
      addChild(id);
    }
  }

  // For sampled internal nodes (with children), add a tip node for each to represent sampling
  std::map<int, int> tipToParentID;
  int tipID = -1;
  for (const Cell& cell : sampledCells) {
    int id = cell.getID();
    Node node = idNodeMap.at(id);
    if (lemon::countOutArcs(cellTree, node) != 0) {
      Node tip = cellTree.addNode();
      cellTree.addArc(node, tip);
      idNodeMap[tipID] = tip;
      tipToParentID[tipID] = id;
      --tipID;
    }
  }

  // Collapse nodes with a single child to simplify the tree
  std::vector<Node> nodesToRemove;
  std::vector<int> idsToRemove;
  for (const auto& pair : idNodeMap) {
    int id = pair.first;
    Node v = pair.second;
    if (v != root && lemon::countOutArcs(cellTree, v) == 1) {
      Node child = cellTree.target(OutArcIt(cellTree, v));
      Node parent = cellTree.source(InArcIt(cellTree, v));
      cellTree.addArc(parent, child);
      nodesToRemove.push_back(v);
      idsToRemove.push_back(id);
    }
  }
  for (Node v : nodesToRemove) cellTree.erase(v);
  for (int id : idsToRemove) idNodeMap.erase(id);

  // Assign string labels and anatomical site labels to nodes
  StringNodeMap cellLabels(cellTree, "");
  StringNodeMap cellAnatomicalSites(cellTree, "");
  for (const auto& pair : idNodeMap) {
    Node node = pair.second;
    int id = pair.first;
    if (id > 0) {
      cellLabels[node] = std::to_string(id);
      cellAnatomicalSites[node] = _anatomicalSiteLabel.at(_idCellMap.at(id).getAnatomicalSite());
    } else {
      int parentID = tipToParentID.at(id);
      cellLabels[node] = std::to_string(parentID) + "_sampled";
      cellAnatomicalSites[node] = _anatomicalSiteLabel.at(_idCellMap.at(parentID).getAnatomicalSite());
    }
  }

  // Update global cell tree pointer
  delete _pCellT;
  _pCellT = new Tree(cellTree, root, cellLabels, cellAnatomicalSites);

  // Update global cell vertex labeling pointer
  delete _cellVertexLabeling;
  _cellVertexLabeling = new StringNodeMap(_pCellT->tree(), "");
  for (NodeIt v(_pCellT->tree()); v != lemon::INVALID; ++v) {
    std::string label = _pCellT->label(v);
    int id;
    try {
      id = std::stoi(label);
    } catch (...) {
      continue; // skip sampled tips
    }
    _cellVertexLabeling->set(v, cellAnatomicalSites[idNodeMap[id]]);
  }

  // Output observable cell mutations for sampled clones (used for clone collapsing)
  std::vector<std::string> usedClones;
  for (int s : _existingSites) {
    if (!_isActiveAnatomicalSite.at(s)) continue;
    for (const auto& kv : _extantCellsByDrivers.at(s)) {
      for (const Cell& cell : kv.second) {
        int id = cell.getID();
        if (idNodeMap.find(id) != idNodeMap.end()) {
          IntSet allMutations = kv.first;
          allMutations.insert(cell.getPassengerMutations().begin(), cell.getPassengerMutations().end());
          IntSet Z;
          std::set_intersection(observableMutations.begin(), observableMutations.end(),
                                allMutations.begin(), allMutations.end(),
                                std::inserter(Z, Z.begin()));
          std::string mutationString;
          for (auto it = Z.begin(); it != Z.end(); ++it) {
            mutationString += std::to_string(_cellTreeMutationToObservableMutation[*it]);
            if (std::next(it) != Z.end()) mutationString += "/";
          }
          if (!mutationString.empty() &&
              std::find(usedClones.begin(), usedClones.end(), mutationString) == usedClones.end()) {
            usedClones.push_back(mutationString);
          }
        }
      }
    }
  }
}

// Writes the Tree in Newick format to the provided output stream
// Branch lengths are calculated using generation times / cell divisions
void Simulation::writeNewick(std::ostream& os, const Tree* pCloneT, const IntIntMap& mutationToGenerationNumber) const {
  const Digraph& tree = pCloneT->tree();

  // Helper to parse mutation IDs from a label string
  auto parseMutations = [](const std::string& label) -> std::vector<int> {
    std::vector<int> mutations;
    std::istringstream iss(label);
    int val;
    while (iss >> val) {
      mutations.push_back(val);
      if (iss.peek() == '/') iss.ignore();
    }
    return mutations;
  };

  // Recursive function to traverse the tree and write Newick format
  std::function<void(Digraph::Node, int)> traverse = [&](Digraph::Node node, int parent_time) {
    const std::string label = pCloneT->label(node);

    // Determine generation time for this node
    int node_time = parent_time;
    std::vector<int> mutations = parseMutations(label);
    for (int mut : mutations) {
      auto it = mutationToGenerationNumber.find(mut);
      if (it != mutationToGenerationNumber.end() && it->second > node_time)
        node_time = it->second;
    }

    if (lemon::countOutArcs(tree, node) == 0) {
      // Leaf node: write label and branch length
      double branchLength = static_cast<double>(_generation - parent_time);
      os << label << ":" << branchLength;
    } else {
      // Internal node: write children recursively
      os << "(";
      bool first = true;
      for (Digraph::OutArcIt arc(tree, node); arc != lemon::INVALID; ++arc) {
        if (!first) os << ",";
        first = false;
        traverse(tree.target(arc), node_time);
      }
      // Branch length for internal node (except for root "GL")
      double branchLength = (label != "0") ? static_cast<double>(node_time - parent_time) : 0.0;
      if (label != "GL") {
        os << ")" << label << ":" << branchLength;
      }
    }
  };

  // Start traversal from the root node
  Digraph::Node rootNode = pCloneT->root();
  traverse(rootNode, 0);
  os << ";" << std::endl;
}
