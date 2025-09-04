#include "migrationgraph.h"

// Default constructor for MigrationGraph
MigrationGraph::MigrationGraph()
  : _G()
  , _root(lemon::INVALID)
  , _nodeToId(_G)
  , _idToNode()
{}

// Constructor that copies a given graph, root node, and node ID map
MigrationGraph::MigrationGraph(const Digraph& G,
                               Node root,
                               const StringNodeMap& id)
: _G()
, _root(lemon::INVALID)
, _nodeToId(_G)
, _idToNode()
{
  // Copy the input graph, root node, and node ID map into this object
  lemon::digraphCopy(G, _G)
    .node(root, _root)
    .nodeMap(id, _nodeToId)
    .run();
  
  // Build reverse mapping from node ID string to node object
  for (NodeIt u(_G); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    assert(_idToNode.count(str) == 0);
    _idToNode[str] = u;
  }
}

// Write the migration graph as a CSV to the output stream
void MigrationGraph::write(std::ostream& out) const
{
  out << "source" << "," << "recipient" << std::endl;
  // Iterate over all arcs and output their source and target node IDs
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    out << _nodeToId[u] << "," << _nodeToId[v] << std::endl;
  }
}

