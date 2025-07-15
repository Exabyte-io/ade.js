import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { TemplateBase, TemplateMixin } from "./templateMixin";
type Base = InMemoryEntity & NamedInMemoryEntity;
export type FlavorMixin = {
    input: {
        name: string;
        templateName?: string;
    }[];
    inputAsTemplates: (TemplateMixin & TemplateBase)[];
    disableRenderMaterials: boolean;
    getInputAsRenderedTemplates: (context: Record<string, unknown>) => Record<string, unknown>[];
};
export declare function flavorMixin(item: Base): FlavorMixin & InMemoryEntity & NamedInMemoryEntity;
export {};
