"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.flavorMixin = flavorMixin;
function flavorMixin(item) {
    // @ts-ignore
    const properties = {
        get input() {
            return this.prop("input", []);
        },
        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
