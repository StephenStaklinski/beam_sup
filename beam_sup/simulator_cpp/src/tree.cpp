#include "tree.h"
#include <lemon/bfs.h>

// Default constructor for Tree
Tree::Tree()
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
  , _leafLabels(_tree)
{
}

// Constructor for Tree from an existing Digraph, root node, and label maps
Tree::Tree(const Digraph& T,
           Node root,
           const StringNodeMap& label,
           const StringNodeMap& l)
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
  , _leafLabels(_tree)
{
  // Copy the input digraph T into _tree, mapping nodes and labels
  lemon::digraphCopy(T, _tree)
    .node(root, _root)
    .nodeMap(label, _nodeToId)
    .run();

  // Build reverse mapping from label to node
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    assert(_idToNode.count(str) == 0);
    _idToNode[str] = u;
  }

  // Initialize leaf sets and other properties
  init();

  // Assign leaf labels if provided (l should be defined on T)
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    // Check if node v is a leaf in T
    if (OutArcIt(T, v) == lemon::INVALID) // check if leaf
    {
      const std::string& label_v = label[v];
      Node u = getNodeByLabel(label_v);
      if (u != lemon::INVALID)
      {
        _leafLabels[u] = l[v];
      }
    }
  }
}

// Write the tree structure as edge list to output stream
void Tree::write(std::ostream& out) const
{
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    out << _nodeToId[u] << " " << _nodeToId[v] << std::endl;
  }
}

// Write leaf node labels to output stream
void Tree::writeLeafLabeling(std::ostream& out) const
{
  for (Node u : _leafSet)
  {
    out << _nodeToId[u] << " " << _leafLabels[u] << std::endl;
  }
}

// Write vertex labels to output stream using provided label map
void Tree::writeVertexLabeling(std::ostream& out,
                               const StringNodeMap& lPlus) const
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    out << _nodeToId[u] << " " << lPlus[u] << std::endl;
  }
}

// Get node in tree by its label string
Node Tree::getNodeByLabel(const std::string& lbl) const
{
  auto it = _idToNode.find(lbl);
  if (it == _idToNode.end())
  {
    return lemon::INVALID;
  }
  else
  {
    return it->second;
  }
}

// Initialize leaf sets, leaf subset, and node levels
void Tree::init()
{
  _leafSet.clear();
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _leafSubset[u].clear();
    OutArcIt a(_tree, u);
    // If node u is a leaf, add to leaf set
    if (a == lemon::INVALID)
    {
      _leafSet.insert(u);
      _isLeaf[u] = true;
    }
  }

  // Compute levels (distance from root) using BFS
  lemon::bfs(_tree).distMap(_level).run(_root);
  // Initialize leaf subsets for all nodes
  initLeafSubset(_root);
}

// Recursively initialize leaf subset for node u
void Tree::initLeafSubset(Node u)
{
  NodeSet merged;
  for (OutArcIt a(_tree, u); a != lemon::INVALID; ++a)
  {
    Node v = _tree.target(a);
    initLeafSubset(v);
    merged.insert(_leafSubset[v].begin(), _leafSubset[v].end());
  }

  // If no children, u is a leaf
  if (merged.empty())
  {
    merged.insert(u);
  }

  _leafSubset[u] = merged;
}
