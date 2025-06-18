import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type FlavorMixin } from "./flavorMixin";
type Base = typeof NamedDefaultableInMemoryEntity & Constructor<FlavorMixin>;
declare const Flavor_base: Base;
export default class Flavor extends Flavor_base {
}
export {};
