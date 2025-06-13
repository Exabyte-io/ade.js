import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";

import { Template } from "./template";

type Base = InMemoryEntity & NamedInMemoryEntity;

export type FlavorMixin = {
    input: { name: string; templateName?: string }[];
    inputAsTemplates: Template[];
    disableRenderMaterials: boolean;
    getInputAsRenderedTemplates: (context: object) => object[];
};

export function flavorMixin(item: Base) {
    // @ts-ignore
    const properties: FlavorMixin & Base = {
        get input() {
            return this.prop<{ name: string; templateName?: string }[]>("input", []);
        },

        // TODO : prevent this from running in client
        get inputAsTemplates() {
            return this.input.map((input) => {
                const template = Template.fromFlavor(
                    this.prop("applicationName"),
                    this.prop("executableName"),
                    input.templateName || input.name,
                );
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
