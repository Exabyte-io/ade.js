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
        get executableId() {
            return this.prop("executableId", "");
        },
        get executableName() {
            return this.prop("executableName", "");
        },
        get applicationName() {
            return this.prop("applicationName", "");
        },
        get supportedApplicationVersions() {
            return this.prop("supportedApplicationVersions");
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
