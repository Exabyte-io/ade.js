/* eslint-disable class-methods-use-this */
import type { UiSchema } from "react-jsonschema-form";

import ContextProvider, { type ContextProviderConfig } from "./ContextProvider";

interface JSONSchemaFormDataProviderConfig extends ContextProviderConfig {
    isUsingJinjaVariables?: boolean;
}

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
export default class JSONSchemaFormDataProvider extends ContextProvider {
    isUsingJinjaVariables: boolean;

    constructor(config: JSONSchemaFormDataProviderConfig) {
        super(config);
        this.isUsingJinjaVariables = Boolean(config?.isUsingJinjaVariables);
    }

    get jsonSchema() {
        throw new Error("Not implemented.");
    }

    get uiSchema(): UiSchema {
        throw new Error("Not implemented.");
    }

    get fields() {
        return {};
    }

    get defaultFieldStyles() {
        return {};
    }

    get uiSchemaStyled(): UiSchema {
        const schema = this.uiSchema;
        return Object.fromEntries(
            Object.entries(schema).map(([key, value]) => [
                key,
                {
                    ...value,
                    ...this.defaultFieldStyles,
                    classNames: `${value.classNames || ""}`,
                },
            ]),
        );
    }
}
