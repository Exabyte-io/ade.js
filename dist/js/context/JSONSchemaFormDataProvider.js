"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const ContextProvider_1 = __importDefault(require("./ContextProvider"));
/**
 * @summary Provides jsonSchema and uiSchema for generating react-jsonschema-form
 *          See https://github.com/mozilla-services/react-jsonschema-form for Form UI.
 *          Form generation example:
 * ```
 * <Form schema={provider.jsonSchema}
 *      uiSchema={provider.uiSchema}
 *      formData={provider.getData(unit.important)} />
 * ```
 */
class JSONSchemaFormDataProvider extends ContextProvider_1.default {
    constructor(config) {
        super(config);
        this.isUsingJinjaVariables = Boolean(config === null || config === void 0 ? void 0 : config.isUsingJinjaVariables);
    }
    get jsonSchema() {
        throw new Error("Not implemented.");
    }
    get uiSchema() {
        throw new Error("Not implemented.");
    }
    get fields() {
        return {};
    }
    get defaultFieldStyles() {
        return {};
    }
    get uiSchemaStyled() {
        const schema = this.uiSchema;
        return Object.fromEntries(Object.entries(schema).map(([key, value]) => [
            key,
            {
                ...value,
                ...this.defaultFieldStyles,
                classNames: `${value.classNames || ""}`,
            },
        ]));
    }
}
exports.default = JSONSchemaFormDataProvider;
