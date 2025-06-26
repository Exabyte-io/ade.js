import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
type Base = InMemoryEntity & NamedInMemoryEntity;
export type FlavorMixin = {
    input: {
        name: string;
        templateName?: string;
    }[];
    disableRenderMaterials: boolean;
    getInputAsRenderedTemplates: (context: Record<string, unknown>) => Record<string, unknown>[];
};
export declare function flavorMixin(item: Base): FlavorMixin & InMemoryEntity & {
    setName(name: string): void;
    name: string;
};
export {};
