"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
exports.executableStaticMixin = executableStaticMixin;
const JSONSchemasInterface_1 = __importDefault(require("@mat3ra/esse/dist/js/esse/JSONSchemasInterface"));
function executableMixin(item) {
    // @ts-expect-error
    const properties = {
        get applicationId() {
            return this.prop("applicationId", []);
        },
        set applicationId(value) {
            this.setProp("applicationId", value);
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
function executableStaticMixin(Executable) {
    const properties = {
        get jsonSchema() {
            return JSONSchemasInterface_1.default.getSchemaById("software/executable");
        },
    };
    Object.defineProperties(Executable, Object.getOwnPropertyDescriptors(properties));
}
