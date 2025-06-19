import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import { runtimeItemsMixin } from "@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import { type FlavorMixin, flavorMixin } from "./flavorMixin";

type Base = typeof NamedDefaultableInMemoryEntity & Constructor<FlavorMixin>;

export default class Flavor extends (NamedDefaultableInMemoryEntity as Base) {}

// Apply mixins
flavorMixin(Flavor.prototype);
runtimeItemsMixin(Flavor.prototype);
