import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import { compareEntitiesInOrderedSetForSorting } from "@mat3ra/code/dist/js/entity/set/ordered/utils";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

export type MaterialsSetContextMixinType = {
    materialsSet: any;
    _materialsSet: any;
    sortMaterialsByIndexInSet: (
        materials?: OrderedInMemoryEntityInSet[],
    ) => OrderedInMemoryEntityInSet[];
    initMaterialsSetContextMixin: () => void;
};

export function materialsSetContextMixin(item: ContextProvider & { _materialsSet: any }) {
    const properties = {
        _materialsSet: undefined,

        get materialsSet() {
            return this._materialsSet;
        },

        initMaterialsSetContextMixin() {
            // @ts-ignore
            this._materialsSet = this.config.context?.materialsSet;
        },

        sortMaterialsByIndexInSet(materials: OrderedInMemoryEntityInSet[] = []) {
            // DO NOT SORT IN PLACE AS IT CHANGES THE ORDER IN `this.materials` AND HAS SIDE EFFECTS (MaterialViewer).
            return materials.concat().sort((a, b) => {
                return compareEntitiesInOrderedSetForSorting(a, b, this.materialsSet._id, false);
            });
        },
    } as MaterialsSetContextMixinType & typeof item;

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}

export function MaterialsSetContextMixin<
    T extends Constructor<ContextProvider & { _materialsSet: any }>,
>(superclass: T) {
    materialsSetContextMixin(superclass.prototype);
    return superclass as T & Constructor<MaterialsSetContextMixinType>;
}
