#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Circulator.hh"

curvenetwork::Vertex::Vertex(int id_)
    : id        (id_)
    , is_deleted()
    , flag      ()
{}
bool curvenetwork::Vertex::is_openend() const {
    return halfedge->prev == 0;
}
bool curvenetwork::Vertex::is_corner() const {
    for (ConstVIHIter h(this); h; ++h)
        if (!h->is_deleted && h->is_corner) return true;
    return false;
}
bool curvenetwork::Vertex::is_boundary() const {
    for (curvenetwork::ConstVOCIter c(this); c; ++c) {
        if (c->is_boundary())
            return true;
    }
    return false;
}
int  curvenetwork::Vertex::valence() const {
    int result = 0;
    for (ConstVOHIter h(this); h; ++h)
        if (!h->is_deleted) ++result;
    return result;
}
