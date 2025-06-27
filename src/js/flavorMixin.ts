import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { FlavorSchema } from "@mat3ra/esse/dist/js/types";

type Base = InMemoryEntity & NamedInMemoryEntity;

type Input = Required<FlavorSchema>["input"];

export type FlavorMixin = {
    input: Input;
    disableRenderMaterials: boolean;
    getInputAsRenderedTemplates: (context: Record<string, unknown>) => Record<string, unknown>[];
};

// TODO: should we add fields from esse schema (executableId, executableName, applicationName)?
export function flavorMixin(item: Base) {
    // @ts-expect-error
    const properties: FlavorMixin & Base = {
        get input() {
            return this.prop<Input>("input", []);
        },

        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return properties;
}
