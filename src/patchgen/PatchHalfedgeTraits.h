#pragma once
#include <kt84/util.h>

namespace patchgen {
	struct PatchHalfedgeTraits {
		struct Data {
			int polychord_id;
			bool polyline_lock;

			Data(){
				polyline_lock = false;
			}

		} patchgen;
	};
}
