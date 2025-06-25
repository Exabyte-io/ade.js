import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { FlavorMixin } from "./flavorMixin";
type BaseFlavor = FlavorMixin & NamedInMemoryEntity & InMemoryEntity;
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;
export declare function executableMixin(item: Base): Base;
export type CreateExecutableConfig = {
    name: string;
    flavors?: Record<string, any>;
};
export type BaseConstructor = Constructor<Base> & {
    constructCustomFlavor?: (config: object) => BaseFlavor;
};
export type ExecutableMixin = {
    toJSON: (exclude?: string[]) => object;
};
export {};
