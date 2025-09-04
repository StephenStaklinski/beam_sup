#ifndef MIGRATIONGRAPH_H
#define MIGRATIONGRAPH_H

#include "tree.h"

/// @brief This class models a migration graph, representing nodes and directed edges.
/// It provides functionality for construction, copying, and serialization.
class MigrationGraph
{
public:
  /// @brief Default constructor. Initializes an empty migration graph.
  MigrationGraph();

  /// @brief Constructs a migration graph from a directed graph, root node, and node identifiers.
  /// @param G Directed graph representing the migration structure.
  /// @param root Root node of the graph.
  /// @param id Mapping from nodes to their string identifiers.
  MigrationGraph(const Digraph& G,
                 Node root,
                 const StringNodeMap& id);

  /// @brief Copy constructor. Creates a deep copy of another MigrationGraph.
  /// @param other MigrationGraph to copy from.
  MigrationGraph(const MigrationGraph& other);

  /// @brief Serializes the migration graph to an output stream.
  /// @param out Output stream to write the graph data.
  /// @brief Serializes the migration graph to an output stream.
  /// @param out Output stream to write the graph data.
  void write(std::ostream& out) const;

private:
  // Directed graph representing migration paths.
  Digraph _G;
  // Root node of the migration graph.
  Node _root;
  // Maps nodes to their string identifiers.
  StringNodeMap _nodeToId;
  // Maps string identifiers to nodes for fast lookup.
  StringToNodeMap _idToNode;
};

#endif
