"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.flavorMixin = flavorMixin;
exports.flavorStaticMixin = flavorStaticMixin;
const JSONSchemasInterface_1 = __importDefault(require("@mat3ra/esse/dist/js/esse/JSONSchemasInterface"));
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
function flavorStaticMixin(Flavor) {
    const properties = {
        get jsonSchema() {
            return JSONSchemasInterface_1.default.getSchemaById("software/flavor");
        },
    };
    Object.defineProperties(Flavor, Object.getOwnPropertyDescriptors(properties));
}
