#ifndef TREE_H
#define TREE_H

#include "utils.h"

/// Models a rooted directed tree with unique node identifiers.
/// Optionally supports leaf labels (e.g., tissue or clone types).
class Tree
{
public:
  /// Constructs an empty tree.
  Tree();

  /// Constructs a tree from a directed graph, root node, node identifiers,
  /// and optional leaf labels.
  ///
  /// @param T Directed graph representing the tree topology.
  /// @param root Root node of the tree.
  /// @param id Map from nodes to unique identifier strings.
  /// @param l (Optional) Map from nodes to leaf labels.
  Tree(const Digraph& T,
      Node root,
      const StringNodeMap& id,
      const StringNodeMap& l);

  /// Copy constructor: creates a deep copy of another Tree.
  Tree(const Tree& other);

  /// Serializes the tree topology and node labels to an output stream.
  void write(std::ostream& out) const;

  /// Writes leaf labeling (if provided) to an output stream.
  void writeLeafLabeling(std::ostream& out) const;

  /// Writes arbitrary vertex labeling to an output stream.
  /// @param lPlus Map from nodes to labels to be written.
  void writeVertexLabeling(std::ostream& out,
                           const StringNodeMap& lPlus) const;

  /// Returns the underlying LEMON directed graph representing the tree.
  const Digraph& tree() const { return _tree; }

  /// Returns the root node of the tree.
  Node root() const { return _root; }

  /// Returns the identifier string of a given node.
  /// @param u Node whose label is requested.
  const std::string& label(Node u) const { return _nodeToId[u]; }

  /// Returns the node corresponding to a given identifier string,
  /// or lemon::INVALID if it does not exist.
  /// @param lbl Identifier string.
  Node getNodeByLabel(const std::string& lbl) const;

  /// Initializes auxiliary data structures (leaf status, leaf sets, etc.).
  void init();

protected:
  /// Recursively initializes the leaf subset rooted at the given node.
  /// @param node Node for which to initialize the leaf subset.
  void initLeafSubset(Node node);
  /// Directed graph representing the tree topology.
  Digraph _tree;
  /// Root node of the tree.
  Node _root;
  /// Map from nodes to unique identifier strings.
  StringNodeMap _nodeToId;
  /// Reverse map: identifier string to node.
  StringToNodeMap _idToNode;
  /// Boolean map indicating whether a node is a leaf.
  BoolNodeMap _isLeaf;
  /// Set containing all leaf nodes.
  NodeSet _leafSet;
  /// Map from nodes to the set of leaf nodes in their subtree.
  NodeNodeSetMap _leafSubset;
  /// Map from nodes to their BFS level (distance from root).
  IntNodeMap _level;
  /// Optional: map from leaf nodes to their labels (e.g., tissue type).
  StringNodeMap _leafLabels;
};

#endif
