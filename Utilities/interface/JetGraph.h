#ifndef JETGRAPH
#define JETGRAPH

#include <vector>
#include <map>
#include <memory>

namespace jetgraph {
  class DiGraph {
  public:
    class Node {
    public:
      Node(int p0_, int p1_, float pt_, float eta_, float phi_, float e_) : 
        p0(p0_), p1(p1_),
        pt(pt_), eta(eta_), phi(phi_), e(e_) { }
      ~Node() { }

      int p0, p1;
      float pt, eta, phi, e;
    };

    DiGraph() { }
    ~DiGraph() { }

    void addNode(int p0, int p1, float pt, float eta, float phi, float e) {
      nodes.emplace_back(p0, p1, pt, eta, phi, e);
    }

    const Node& getnode(int idx) const { return nodes.at(idx); }

    void clear() { nodes.clear(); }

  private:
    std::vector<Node> nodes;
  }; 
};

#endif
