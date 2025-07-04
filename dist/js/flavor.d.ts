import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import { type RuntimeItemsInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type FlavorMixin } from "./flavorMixin";
type Base = typeof NamedDefaultableInMemoryEntity & Constructor<FlavorMixin> & Constructor<RuntimeItemsInMemoryEntity>;
declare const Flavor_base: Base;
export default class Flavor extends Flavor_base {
}
export {};
