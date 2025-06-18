import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
export type MaterialsContextMixinType = {
    materials: (MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet)[];
    initMaterialsContextMixin: () => void;
};
export declare function materialsContextMixin(item: ContextProvider & {
    _materials: (MaterialMixin & InMemoryEntity)[];
}): void;
export declare function MaterialsContextMixin<T extends Constructor<ContextProvider & {
    _materials: (MaterialMixin & InMemoryEntity)[];
}>>(superclass: T): T & Constructor<MaterialsContextMixinType>;
