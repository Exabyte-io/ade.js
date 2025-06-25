import { type ApplicationName } from "@exabyte-io/application-flavors.js";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import Executable from "./executable";
import { CreateApplicationConfig } from "./tree";
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;
export type BaseConstructor = Constructor<Base> & {
    constructCustomExecutable?: (config: object) => Executable;
};
export type ApplicationConstructor = Constructor<ApplicationMixin> & ApplicationStaticMixin;
export type ApplicationMixin = {
    summary: string | undefined;
    version: string;
    build: string | undefined;
    shortName: string;
    name: ApplicationName;
    hasAdvancedComputeOptions: boolean;
    isLicensed: boolean;
    isUsingMaterial: boolean;
};
export type ApplicationStaticMixin = {
    defaultConfig: {
        name: string;
        shortName: string;
        version: string;
        summary: string;
        build: string;
    };
    create: (config: CreateApplicationConfig) => Base;
    getUniqueAvailableNames: () => string[];
};
export declare function applicationMixin(item: Base): Base;
export declare function applicationStaticMixin<T extends BaseConstructor>(Application: T): ApplicationStaticMixin;
export {};
