"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.materialContextMixin = materialContextMixin;
exports.MaterialContextMixin = MaterialContextMixin;
function materialContextMixin(item) {
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
            var _a;
            // @ts-ignore
            const ConstructorApplication = this.constructor.Application;
            if (!ConstructorApplication) {
                throw Error("ApplicationContextMixin: Application is undefined");
            }
            // @ts-ignore
            const application =
                (_a = this.config.context) === null || _a === void 0 ? void 0 : _a.application;
            this._application = application || ConstructorApplication.createDefault();
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
function MaterialContextMixin(superclass) {
    materialContextMixin(superclass.prototype);
    return superclass;
}
