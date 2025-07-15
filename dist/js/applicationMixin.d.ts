import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import Executable from "./executable";
import { CreateApplicationConfig } from "./tree";
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;
export declare function applicationMixin(item: Base): Base;
export type BaseConstructor = Constructor<Base> & {
    constructCustomExecutable?: (config: object) => Executable;
};
export declare function applicationStaticMixin<T extends BaseConstructor>(Application: T): ApplicationStaticProperties;
export type ApplicationStaticProperties = {
    defaultConfig: {
        name: string;
        shortName: string;
        version: string;
        summary: string;
        build: string;
    };
    create: (config: CreateApplicationConfig) => Base;
    createFromNameVersionBuild: (config: CreateApplicationConfig) => Base;
    getUniqueAvailableNames: () => string[];
    constructExecutable: (config: object) => Executable;
};
export type ApplicationConstructor = Constructor<ApplicationMixin> & ApplicationStaticProperties;
export type ApplicationMixin = {
    defaultExecutable: Executable;
    summary: string | undefined;
    version: string;
    build: string | undefined;
    shortName: string;
    executables: Executable[];
    hasAdvancedComputeOptions: boolean;
    isLicensed: boolean;
    isUsingMaterial: boolean;
    getExecutableByName: (name?: string) => Executable;
    getExecutableByConfig: (config?: {
        name: string;
    }) => Executable;
    getExecutables: () => Executable[];
    getBuilds: () => string[];
    getVersions: () => string[];
};
export {};
