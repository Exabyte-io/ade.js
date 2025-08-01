import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import { runtimeItemsMixin } from "@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import { type ExecutableMixin, executableMixin } from "./executableMixin";

type Base = typeof NamedDefaultableInMemoryEntity & Constructor<ExecutableMixin>;

export default class Executable extends (NamedDefaultableInMemoryEntity as Base) {}

// Apply mixins
runtimeItemsMixin(Executable.prototype);
executableMixin(Executable.prototype);
