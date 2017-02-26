#include "Patch.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include <numeric>
#include <kt84/util.h>
using namespace std;
using namespace Eigen;
using namespace kt84;
using curvenetwork::Patch;
using curvenetwork::Halfchain;

double              Patch::quadSize;
bool                Patch::use_even_num_subdiv = true;
bool                Patch::prefer_rect_proc3   = true;

bool Patch::is_valid_for_patch(Halfchain*& halfchain) {
    Halfchain* c_corner = nullptr;
    for (auto c = halfchain; ; ) {
        if (!c)
            return false;
        if (c->is_corner())
            c_corner = c;
        c = c->next();
        if (c == halfchain)
            break;
    }
    if (!c_corner)          // Patch must have at least one corner
        return false;
    halfchain = c_corner->next();
    return true;
}
int Patch::num_corners(Halfchain* halfchain) {
    int result = 0;
    for (auto c = halfchain; ; ) {
        if (c->is_corner())
            ++result;
        c = c->next();
        if (c == halfchain)
            break;
    }
    return result;
}
int Patch::num_subdiv(Halfchain* halfchain, int side_index) {
    int current_index = 0;
    int result        = 0;
    for (auto c = halfchain; ; ) {
        result += c->num_subdiv();
        if (c->is_corner()) {
            if (side_index == current_index)
                return result;
            ++current_index;
            result = 0;
        }
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    return 0;
}
VectorXi Patch::num_subdiv(Halfchain* halfchain) {
    int n = num_corners(halfchain);
    VectorXi l(n);
    for (int i = 0; i < n; ++i)
        l[i] = num_subdiv(halfchain, i);
    return l;
}
int     Patch::num_free_halfchains(Halfchain* halfchain, int side_index) {
    int current_index = 0;
    int result        = 0;
    for (auto c = halfchain; ; ) {
        if (c->num_subdiv() == 0)
            ++result;
        if (c->is_corner()) {
            if (current_index == side_index)
                return result;
            ++current_index;
            result = 0;
        }
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    return -1;
}
double  Patch::side_length        (Halfchain* halfchain, int side_index) {
    int current_index = 0;
    double result     = 0;
    for (auto c = halfchain; ; ) {
        result += c->length();
        if (c->is_corner()) {
            if (side_index == current_index)
                return result;
            
            ++current_index;
            result = 0;
        }
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    return 0;
}
void    Patch::change_num_subdiv  (Halfchain* halfchain, int side_index, int target_num_subdiv) {
    int delta_target = target_num_subdiv - num_subdiv(halfchain, side_index);
    if (delta_target <= 0)
        return;
    vector<Halfchain*> c_free;
    int current_index = 0;
    for (auto c = halfchain; ; ) {
        if (current_index == side_index && c->num_subdiv() == 0)
            c_free.push_back(c);
        if (c->is_corner()) {
            if (current_index == side_index)
                break;
            ++current_index;
        }
        c = c->next();
        if (c == halfchain)
            break;      // this should not happen!
    }
    assert(!c_free.empty());
    double length_sum = accumulate(
        c_free.begin(), c_free.end(), 0.0,
        [] (double result, Halfchain* c) { return result + c->length(); });
    for_each(c_free.begin(), c_free.end(),
        [delta_target, length_sum] (Halfchain* c) {
            c->edgechain->num_subdiv = static_cast<int>(floor(0.5 + delta_target * c->length() / length_sum));
            if (use_even_num_subdiv && c->edgechain->num_subdiv % 2)
                // ensure even num_subdiv option
                ++c->edgechain->num_subdiv;
    });
    // not sure if this is really useful...
    //delta_target = target_num_subdiv - num_subdiv(side_index);
    //
    //if (delta_target != 0)
    //    c_free.front()->edgechain->num_subdiv += delta_target;
}
void Patch::set_default_num_subdiv(Halfchain* halfchain) {
    // set default num_subdiv for unspecified edgechains
    int num_sides = num_corners(halfchain);
    for (int i = 0; i < num_sides; ++i) {
        int j = (i + 2) % num_sides;
        bool is_fixed_i = num_free_halfchains(halfchain, i) == 0;
        bool is_fixed_j = num_free_halfchains(halfchain, j) == 0;
        if (!is_fixed_i && !is_fixed_j) {
            int target_num_subdiv = util::max<int>(
                static_cast<int>(side_length(halfchain, i) / quadSize),
                static_cast<int>(side_length(halfchain, j) / quadSize),
                num_subdiv(halfchain, i) + num_free_halfchains(halfchain, i),
                num_subdiv(halfchain, j) + num_free_halfchains(halfchain, j));
            change_num_subdiv(halfchain, i, target_num_subdiv);
            change_num_subdiv(halfchain, j, target_num_subdiv);
        } else if (!is_fixed_i) {
            change_num_subdiv(halfchain, i, num_subdiv(halfchain, j));
        } else if (!is_fixed_j) {
            change_num_subdiv(halfchain, j, num_subdiv(halfchain, i));
        }
    }
}
