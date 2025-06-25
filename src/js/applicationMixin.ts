import { type ApplicationName, allApplications } from "@exabyte-io/application-flavors.js";
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
    // defaultExecutable: Executable;
    summary: string | undefined;
    version: string;
    build: string | undefined;
    shortName: string;
    name: ApplicationName;
    // executables: Executable[];
    hasAdvancedComputeOptions: boolean;
    isLicensed: boolean;
    isUsingMaterial: boolean;
    // getExecutableByName: (name?: string) => Executable;
    // getExecutableByConfig: (config?: { name: string }) => Executable;
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

export function applicationMixin(item: Base) {
    // @ts-expect-error
    const properties: ApplicationMixin & Base = {
        get summary() {
            return this.prop("summary");
        },

        get version() {
            return this.prop("version", "");
        },

        get build() {
            return this.prop("build");
        },

        get shortName() {
            return this.prop("shortName", this.name);
        },

        get hasAdvancedComputeOptions() {
            return this.prop("hasAdvancedComputeOptions", false);
        },

        get isLicensed() {
            return this.prop("isLicensed", false);
        },

        get isUsingMaterial() {
            const materialUsingApplications = ["vasp", "nwchem", "espresso", "exabyteml"];
            return materialUsingApplications.includes(this.name);
        },

        // get executables() {
        //     const tree = getAppTree(this.name as ApplicationName);
        //     return Object.keys(tree)
        //         .filter((key) => {
        //             const { supportedApplicationVersions } = tree[key];
        //             return (
        //                 !supportedApplicationVersions ||
        //                 supportedApplicationVersions.includes(this.version)
        //             );
        //         })
        //         .map((key) => {
        //             return (
        //                 this.constructor as unknown as ApplicationStaticMixin
        //             ).constructExecutable({
        //                 ...tree[key],
        //                 name: key,
        //             });
        //         });
        // },

        // get defaultExecutable() {
        //     return this.getExecutableByName();
        // },

        // getExecutableByName(name?: string) {
        //     return (this.constructor as unknown as ApplicationStaticMixin).constructExecutable(
        //         getExecutableConfig({
        //             appName: this.name as ApplicationName,
        //             execName: name,
        //         }),
        //     );
        // },

        // getExecutableByConfig(config?: { name: string }) {
        //     return config ? this.getExecutableByName(config.name) : this.defaultExecutable;
        // },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return item;
}

export function applicationStaticMixin<T extends BaseConstructor>(Application: T) {
    const properties: ApplicationStaticMixin = {
        get defaultConfig() {
            return {
                name: "espresso",
                shortName: "qe",
                version: "6.3",
                summary: "Quantum Espresso",
                build: "Default",
            };
        },

        create({ name, version = null, build = "Default" }: CreateApplicationConfig) {
            return new Application({ name, version, build });
        },

        getUniqueAvailableNames() {
            return allApplications as string[];
        },
    };

    Object.defineProperties(Application, Object.getOwnPropertyDescriptors(properties));

    return properties;
}
