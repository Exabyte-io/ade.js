import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import {
    type RuntimeItemsInMemoryEntity,
    runtimeItemsMixin,
} from "@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import { type FlavorMixin, flavorMixin, flavorStaticMixin } from "./flavorMixin";

type Base = typeof NamedDefaultableInMemoryEntity &
    Constructor<FlavorMixin> &
    Constructor<RuntimeItemsInMemoryEntity>;

export default class Flavor extends (NamedDefaultableInMemoryEntity as Base) {}

// Apply mixins
flavorMixin(Flavor.prototype);
runtimeItemsMixin(Flavor.prototype);
flavorStaticMixin(Flavor);
