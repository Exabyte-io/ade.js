import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import {
    type ApplicationMixin,
    applicationMixin,
    applicationStaticMixin,
} from "./applicationMixin";

type Base = typeof NamedDefaultableInMemoryEntity & Constructor<ApplicationMixin>;

export default class Application extends (NamedDefaultableInMemoryEntity as Base) {}

applicationMixin(Application.prototype);
applicationStaticMixin(Application);
