import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import {
    type ApplicationMixin,
    applicationMixin,
    applicationStaticMixin,
} from "./applicationMixin";
import { type CreateApplicationConfig, getApplicationConfig } from "./tree";

type Base = typeof NamedDefaultableInMemoryEntity & Constructor<ApplicationMixin>;

export default class Application extends (NamedDefaultableInMemoryEntity as Base) {
    constructor(config: CreateApplicationConfig) {
        const staticConfig = getApplicationConfig(config);
        super({ ...staticConfig, ...config });
    }
}

applicationMixin(Application.prototype);
applicationStaticMixin(Application);
