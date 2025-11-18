#include "simulation.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

// Global variables for migration matrix filename and migration rate string
std::string filenameMigrationMatrix;

// Output results of the simulation to files in the specified directory
void output(const Simulation& simulation,
      const std::string& outputDirectory,
      int seed)
{
  char buf[1024];

  // Check if output directory exists; create it if not
  struct stat info;
  if (stat(outputDirectory.c_str(), &info) != 0) {
    // Directory does not exist, attempt to create it
    if (mkdir(outputDirectory.c_str(), 0755) != 0) {
      std::cerr << "ERROR: Could not create output directory: " << outputDirectory << std::endl;
      exit(1);
    }
  } else if (!(info.st_mode & S_IFDIR)) {
    // Path exists but is not a directory
    std::cerr << "ERROR: Output path exists but is not a directory: " << outputDirectory << std::endl;
    exit(1);
  }

  // Get cell tree from the simulation to trace cell lineage
  const Tree& CT = simulation.getCellTree();

  // Output cell tree in Newick format for phylogenetic analysis
  snprintf(buf, 1024, "%s/cell_tree_seed%d.nwk", outputDirectory.c_str(), seed);
  std::ofstream outCT_newick(buf);
  const std::map<int, int> cellToGeneration = simulation.getCellToGenerationNumber();
  simulation.writeNewick(outCT_newick, &CT, cellToGeneration);
  outCT_newick.close();

  // Output leaf labeling for the cell tree
  snprintf(buf, 1024, "%s/cell_tree_seed%d.labeling", outputDirectory.c_str(), seed);
  std::ofstream outCT_labeling(buf);
  CT.writeLeafLabeling(outCT_labeling);
  outCT_labeling.close();

  // Output vertex labeling (all internal nodes and leaves) for the cell tree
  snprintf(buf, 1024, "%s/cell_tree_seed%d.vertex.labeling", outputDirectory.c_str(), seed);
  std::ofstream outCT_vertex_labeling(buf);
  const StringNodeMap& cPlus = simulation.getCellVertexLabeling();
  CT.writeVertexLabeling(outCT_vertex_labeling, cPlus);
  outCT_vertex_labeling.close();
}

int main(int argc, char** argv)
{
  // Simulation parameters with default values
  int seed = 0;
  double K = 5e4; // Carrying capacity
  std::string migrationRateString = "1e-6"; // Migration rate
  double mutFreqThreshold = 0.05; // Mutation frequency threshold
  double mutationRate = 0.1; // Mutation rate
  double driverProb = 2e-7; // Driver mutation probability
  int maxNrAnatomicalSites = -1; // Maximum number of anatomical sites
  int maxGenerations = 250; // Maximum number of cell division generations
  std::string outputDirectory = "."; // Output directory
  int N = -1; // Number of successful simulations
  int downsampleCellNumber = 100; // Number of cells to downsample to
  int inputNumPossibleAnatomicalSites = 10; // Number of possible anatomical sites for default migration matrix
  int migrationStartGeneration = 0; // Generation to start migrations
  int migrationEndGeneration = -1; // Last generation to perform migrations
  bool keepPolytomies = false; // Whether to keep polytomies in the cell tree
  
  // Parse command line arguments using Lemon ArgParser
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
  .refOption("ns", "Number of possible anatomical sites for default uniform migration matrix when not providing one as input (default: 10)", inputNumPossibleAnatomicalSites)
  .refOption("m", "Filepath to a csv file of a square migration probability matrix with rows sum to 1 and diagonal elements 0 (default: Uniform probabilities)", filenameMigrationMatrix, false)
  .refOption("d", "Driver probability (default: 1e-7)", driverProb)
  .refOption("f", "Mutation frequency threshold (default: 0.05)", mutFreqThreshold)
  .refOption("k", "Carrying capacity (default: 5e4)", K)
  .refOption("mig", "Comma-sep migration rates (default: 1e-6 uniform rate for all source tissues)", migrationRateString)
  .refOption("mut", "Mutation rates (default: 0.1)", mutationRate)
  .refOption("n", "Number of successful simulations (default: -1)", N)
  .refOption("c", "Maximum number of cell division generations to end the simulation after (default: 250)", maxGenerations)
  .refOption("m", "Maximum number of detectable anatomical sites (default: -1)", maxNrAnatomicalSites)
  .refOption("d", "Number of cells to randomly downsample cell tree until reached (default: -1)", downsampleCellNumber)
  .refOption("gs", "Generation to start migrations (default: 0, meaning migrations start from the beginning)", migrationStartGeneration)
  .refOption("ge", "Last generation for migrations (default: -1, meaning migrations continue until the end)", migrationEndGeneration)
  .refOption("p", "Keep polytomies in the resulting cell tree (default: false)", keepPolytomies)
  .refOption("o", "Output directory (default: '.')", outputDirectory);
  ap.parse();

  // Round about way to set resolvePolytomies due to ArgParser behavior with bool variables always being false by default
  bool resolvePolytomies = true;
  if (keepPolytomies) {
    resolvePolytomies = false;
  }

  // Read migration transition probabilities from file if provided
  std::vector<std::vector<double>> migrationTransitionProbs;
  if (!filenameMigrationMatrix.empty())
  {
    std::ifstream file(filenameMigrationMatrix);
    if (!file.is_open()) {
      std::cerr << "Could not open the file " << filenameMigrationMatrix << std::endl;
      return 1;
    } else {
      std::cout << "Found user input migration matrix file, using provided transition probabilities." << std::endl;
      // Parse CSV file into migrationTransitionProbs matrix
      std::string line;
      while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        std::string token;
        while (std::getline(ss, token, ',')) {
          row.push_back(std::stod(token));
        }
        migrationTransitionProbs.push_back(row);
      }

      // Verify that diagonal elements are zero (no self-migration)
      for (int i = 0; i < migrationTransitionProbs.size(); ++i) {
      if (migrationTransitionProbs[i][i] != 0) {
        // Print matrix and error message
        for (const auto& row : migrationTransitionProbs) {
          for (const auto& value : row) {
            std::cout << value << " ";
          }
          std::cout << std::endl;
        }
        std::cout << "ERROR: Self migrations are not allowed! Fix all transition probabilities on the diagonal to be 0.0." << std::endl;
        return 1;
      }
      }

      // Verify that each row sums to 1 (valid probability distribution)
      for (const auto& row : migrationTransitionProbs) {
        double sum = 0.0;
        for (const auto& value : row) {
          sum += value;
        }
        if (sum != 0.0) {
          if (std::abs(sum - 1.0) > 1e-6) {
            // Print matrix and error message
            for (const auto& row : migrationTransitionProbs) {
              for (const auto& value : row) {
              std::cout << value << " ";
              }
              std::cout << std::endl;
            }
            std::cout << "sum: " << sum << std::endl;
            std::cout << "ERROR: Transition probabilities in a row do not sum to 1.0." << std::endl;
            return 1;
          }
        }
      }
      // Print transition matrix for confirmation
      for (const auto& row : migrationTransitionProbs) {
        for (const auto& value : row) {
          std::cout << value << " ";
        }
        std::cout << std::endl;
      }
    }
  } else {
    // No migration matrix file provided; use uniform transition probabilities
    std::cout << "No user input migration matrix file, using uniform transition probabilities." << std::endl;
  }

  // Read tissue-specific migration rates
  std::vector<double> migrationRates;
  std::stringstream ss(migrationRateString);
  std::string token;
  while (std::getline(ss, token, ',')) {
    migrationRates.push_back(std::stod(token));
  }

  // Repeat simulation if it fails (e.g., all cells die)
  bool run = true;
  while (run == true) {
    // Set random seed
    g_rng = std::mt19937(seed);

    std::cout << "Starting simulation with seed " << seed << "." << std::endl;

    // Create simulation object with parameters
    Simulation simulation(K,
              migrationRates,
              mutationRate,
              driverProb,
              mutFreqThreshold,
              maxNrAnatomicalSites,
              maxGenerations,
              downsampleCellNumber,
              inputNumPossibleAnatomicalSites,
              migrationTransitionProbs,
              migrationStartGeneration,
              migrationEndGeneration,
              resolvePolytomies);

    // Run simulation
    if (simulation.simulate())
    {
      // Output results if simulation succeeds
      std::cout << "Simulation successful, outputting results." << std::endl;
      output(simulation, outputDirectory, seed);
      run = false;
    } else {
      // If simulation fails, draw a new positive seed and retry
      std::random_device rd;
      seed = rd();
      while (seed < 0) {
      seed = rd();
      }
      std::cout << "Simulation failed, drawing new seed and retrying." << std::endl;
    }
  }
  
  return 0;
}
