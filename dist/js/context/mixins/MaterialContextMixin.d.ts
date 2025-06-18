import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import type { ApplicationMixin } from "src/js/applicationMixin";
export type MaterialContextMixinType = {
    assertMaterial: () => void;
    isEditedIsSetToFalseOnMaterialUpdate?: boolean;
    updateMaterialHash: () => void;
    isMaterialCreatedDefault: boolean;
    isMaterialUpdated: boolean;
    material: MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet;
    extraData?: {
        materialHash: string;
    };
    initMaterialContextMixin: () => void;
    _application: ApplicationMixin;
};
export declare function materialContextMixin(item: ContextProvider & {
    _material?: MaterialMixin & InMemoryEntity;
}): void;
export declare function MaterialContextMixin<T extends Constructor<ContextProvider & {
    _material: MaterialMixin & InMemoryEntity;
}>>(superclass: T): T & Constructor<MaterialContextMixinType>;
