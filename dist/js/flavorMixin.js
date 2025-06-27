"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.flavorMixin = flavorMixin;
// TODO: should we add fields from esse schema (executableId, executableName, applicationName)?
function flavorMixin(item) {
    // @ts-expect-error
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
