//
// Created by Andrey Bzikadze on 11/10/21.
//

#pragma once

#include "error_correction/multiplicity_estimation.hpp"
#include "mdbg_topology.hpp"
#include "paths.hpp"

namespace repeat_resolution {

class MultiplexDBG
    : public graph_lite::Graph<
          /*typename NodeType=*/RRVertexType,
          /*typename NodePropType=*/RRVertexProperty,
          /*typename EdgePropType=*/RREdgeProperty,
          /*EdgeDirection direction=*/graph_lite::EdgeDirection::DIRECTED,
          /*MultiEdge multi_edge=*/graph_lite::MultiEdge::ALLOWED,
          /*SelfLoop self_loop=*/graph_lite::SelfLoop::ALLOWED,
          /*Map adj_list_spec=*/graph_lite::Map::UNORDERED_MAP,
          /*Container neighbors_container_spec=*/
          graph_lite::Container::MULTISET> {
  friend class MultiplexDBGIncreaser;
  uint64_t next_edge_index{0};
  uint64_t next_vert_index{0};
  uint64_t n_iter{0};
  std::unordered_map<RRVertexType, RREdgeProperty> isolate_properties;

  void freeze_isolated_loops();

  void assert_validity() const;

public:
  // TODO remove public
  RRPaths *rr_paths;

  // This constructor is for testing purposes
  MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges, uint64_t start_k,
               RRPaths *rr_paths);

  MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *rr_paths, uint64_t start_k,
               UniqueClassificator &classificator, bool debug,
               const std::experimental::filesystem::path &dir,
               logging::Logger &logger);

  MultiplexDBG(const MultiplexDBG &) = delete;
  MultiplexDBG(MultiplexDBG &&) = default;
  MultiplexDBG &operator=(const MultiplexDBG &) = delete;
  MultiplexDBG &operator=(MultiplexDBG &&) = default;

  void serialize_to_dot(const std::experimental::filesystem::path &path) const;

  [[nodiscard]] bool is_frozen() const;

  RRVertexType get_new_vertex(uint64_t len);

  void freeze_vertex(const RRVertexType &vertex) {
    RRVertexProperty &prop = node_prop(vertex);
    prop.freeze();
  }

  void move_edge(const RRVertexType &s1, NeighborsIterator e1_it,
                 const RRVertexType &s2, const RRVertexType &e2);

  void merge_edges(const RRVertexType &s1, NeighborsIterator e1_it,
                   const RRVertexType &s2, NeighborsIterator e2_it,
                   uint64_t overlap_len);

  EdgeIndexType add_connecting_edge(NeighborsIterator e1_it,
                                    const RRVertexType &s2,
                                    NeighborsIterator e2_it);
};

} // End namespace repeat_resolution