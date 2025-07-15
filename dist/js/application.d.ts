import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type ApplicationMixin, type ApplicationStaticMixin } from "./applicationMixin";
type Base = typeof NamedDefaultableInMemoryEntity & Constructor<ApplicationMixin> & ApplicationStaticMixin;
declare const Application_base: Base;
export default class Application extends Application_base {
}
export {};
