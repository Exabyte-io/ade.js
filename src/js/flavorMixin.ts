import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";

type Base = InMemoryEntity & NamedInMemoryEntity;

export type FlavorMixin = {
    input: { name: string; templateName?: string }[];
    disableRenderMaterials: boolean;
    getInputAsRenderedTemplates: (context: Record<string, unknown>) => Record<string, unknown>[];
};

export function flavorMixin(item: Base) {
    // @ts-ignore
    const properties: FlavorMixin & Base = {
        get input() {
            return this.prop<{ name: string; templateName?: string }[]>("input", []);
        },

        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return properties;
}
