"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.flavorMixin = flavorMixin;
const template_1 = __importDefault(require("./template"));
function flavorMixin(item) {
    // @ts-ignore
    const properties = {
        get input() {
            return this.prop("input", []);
        },
        // TODO : prevent this from running in client
        get inputAsTemplates() {
            return this.input.map((input) => {
                const template = template_1.default.fromFlavor(this.prop("applicationName", ""), this.prop("executableName", ""), input.templateName || input.name);
                template.name = input.name;
                return template;
            });
        },
        getInputAsRenderedTemplates(context) {
            return this.inputAsTemplates.map((t) => t.getRenderedJSON(context));
        },
        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
