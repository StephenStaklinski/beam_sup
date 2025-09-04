#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <random>
#include <cassert>

// Type alias for directed graph using Lemon library
typedef lemon::ListDigraph Digraph;
// Defines types for nodes, arcs, etc., for Digraph
DIGRAPH_TYPEDEFS(Digraph);
// Map from Digraph nodes to strings
typedef Digraph::NodeMap<std::string> StringNodeMap;
// Map from strings to Digraph nodes
typedef std::map<std::string, Node> StringToNodeMap;
// Set of Digraph nodes
typedef std::set<Node> NodeSet;
// Map from Digraph nodes to sets of nodes
typedef Digraph::NodeMap<NodeSet> NodeNodeSetMap;
// Vector of integers
typedef std::vector<int> IntVector;
// Vector of doubles
typedef std::vector<double> DoubleVector;
// Set of integers
typedef std::set<int> IntSet;
// Subgraph view of a Digraph (const version)
typedef lemon::SubDigraph<const Digraph> SubDigraph;
/// Random number generator
extern std::mt19937 g_rng;

#endif
