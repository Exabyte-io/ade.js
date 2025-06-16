import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import type { ApplicationMixin } from "src/js/applicationMixin";

import Application from "../../application";

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

export function materialContextMixin(
    item: ContextProvider & { _material?: MaterialMixin & InMemoryEntity },
) {
    const properties = {
        assertMaterial() {
            if (!this.material) {
                throw new Error("Material is not set");
            }
        },

        updateMaterialHash() {
            this.assertMaterial();
            if (this.isEditedIsSetToFalseOnMaterialUpdate) this.isEdited = false;
            this.extraData = { materialHash: this.material.hash };
        },

        // Workaround: Material.createDefault() used to initiate workflow reducer and hence here too
        //  does not have an id. Here we catch when such material is used and avoid resetting isEdited
        get isMaterialCreatedDefault() {
            this.assertMaterial();
            return !this.material.id;
        },

        get isMaterialUpdated() {
            this.assertMaterial();
            return Boolean(this.extraData && this.extraData.materialHash !== this.material.hash);
        },

        get material() {
            this.assertMaterial();
            return this._material;
        },

        initMaterialContextMixin() {
            // @ts-ignore
            const ConstructorApplication = this.constructor.Application as typeof Application;

            if (!ConstructorApplication) {
                throw Error("ApplicationContextMixin: Application is undefined");
            }

            // @ts-ignore
            const application = this.config.context?.application as ApplicationMixin;

            this._application = application || ConstructorApplication.createDefault();
        },
    } as MaterialContextMixinType & typeof item;

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}

export function MaterialContextMixin<
    T extends Constructor<ContextProvider & { _material: MaterialMixin & InMemoryEntity }>,
>(superclass: T) {
    materialContextMixin(superclass.prototype);
    return superclass as T & Constructor<MaterialContextMixinType>;
}
