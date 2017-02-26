#pragma once
#include "PatchVertexTag.h"

namespace patchgen {
    struct PatchVertexTraits {
        struct Data {
            int corner_index;
			int irregular_index;
			int i_polychord_he_h;

            PatchVertexTag tag;            // later used when finding vertices in the core pattern
            Data()
                : corner_index(-1)
                , tag(PatchVertexTag::None)
            {}
        } patchgen;
    };
}