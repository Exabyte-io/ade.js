import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { FlavorSchema } from "@mat3ra/esse/dist/js/types";
type Base = InMemoryEntity & NamedInMemoryEntity;
type Input = Required<FlavorSchema>["input"];
export type FlavorMixin = {
    input: Input;
    disableRenderMaterials: boolean;
    executableId: string;
    executableName: string;
    applicationName: string;
    supportedApplicationVersions?: string[];
    getInputAsRenderedTemplates: (context: Record<string, unknown>) => Record<string, unknown>[];
};
export declare function flavorMixin(item: Base): FlavorMixin & InMemoryEntity & NamedInMemoryEntity;
export {};
