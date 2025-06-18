import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
export type MaterialsSetContextMixinType = {
    materialsSet: any;
    _materialsSet: any;
    sortMaterialsByIndexInSet: (materials?: OrderedInMemoryEntityInSet[]) => OrderedInMemoryEntityInSet[];
    initMaterialsSetContextMixin: () => void;
};
export declare function materialsSetContextMixin(item: ContextProvider & {
    _materialsSet: any;
}): void;
export declare function MaterialsSetContextMixin<T extends Constructor<ContextProvider & {
    _materialsSet: any;
}>>(superclass: T): T & Constructor<MaterialsSetContextMixinType>;
