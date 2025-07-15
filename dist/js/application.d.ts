import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type ApplicationMixin } from "./applicationMixin";
import { type CreateApplicationConfig } from "./tree";
type Base = typeof NamedDefaultableInMemoryEntity & Constructor<ApplicationMixin>;
declare const Application_base: Base;
export default class Application extends Application_base {
    constructor(config: CreateApplicationConfig);
}
export {};
