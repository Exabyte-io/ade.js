import {
    type ApplicationName,
    allApplications,
    getAppData,
    getAppTree,
} from "@exabyte-io/application-flavors.js";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import lodash from "lodash";

import { Executable } from "./executable";
import { CreateApplicationConfig, getExecutableConfig } from "./tree";

type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;

export function applicationMixin(item: Base) {
    // @ts-ignore
    const properties: ApplicationMixin & Base = {
        get defaultExecutable() {
            return this.getExecutableByName();
        },

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

        get executables() {
            const tree = getAppTree(this.name);
            return Object.keys(tree)
                .filter((key) => {
                    const { supportedApplicationVersions } = tree[key];
                    return (
                        !supportedApplicationVersions ||
                        supportedApplicationVersions.includes(this.version)
                    );
                })
                .map((key) => {
                    return (
                        this.constructor as unknown as ApplicationStaticProperties
                    ).constructExecutable({
                        ...tree[key],
                        name: key,
                    });
                });
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

        getExecutableByName(name?: string) {
            return (this.constructor as unknown as ApplicationStaticProperties).constructExecutable(
                getExecutableConfig({
                    appName: this.name,
                    execName: name,
                }),
            );
        },

        getExecutableByConfig(config?: { name: string }) {
            return config ? this.getExecutableByName(config.name) : this.defaultExecutable;
        },

        getExecutables() {
            return this.executables;
        },

        getBuilds() {
            const data = getAppData(this.name as ApplicationName);
            const { versions } = data;
            const builds = ["Default"];
            versions.map((v) => v.build && builds.push(v.build));
            return lodash.uniq(builds);
        },

        getVersions() {
            const data = getAppData(this.name as ApplicationName);
            const { versions } = data;
            const these = versions.map((v) => v.version);
            return lodash.uniq(these);
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return item;
}

export type BaseConstructor = Constructor<Base> & {
    constructCustomExecutable?: (config: object) => Executable;
};

export function applicationStaticMixin<T extends BaseConstructor>(Application: T) {
    const properties: ApplicationStaticProperties = {
        get defaultConfig() {
            return {
                name: "espresso",
                shortName: "qe",
                version: "6.3",
                summary: "Quantum Espresso",
                build: "Default",
            };
        },

        create(config: CreateApplicationConfig) {
            return this.createFromNameVersionBuild(config);
        },

        createFromNameVersionBuild({
            name,
            version = null,
            build = "Default",
        }: CreateApplicationConfig) {
            return new Application({ name, version, build });
        },

        getUniqueAvailableNames() {
            return allApplications as string[];
        },

        constructExecutable(this: BaseConstructor & typeof properties, config: object) {
            if (this.constructCustomExecutable) {
                return this.constructCustomExecutable(config);
            }
            return new Executable(config);
        },
    };

    Object.defineProperties(Application, Object.getOwnPropertyDescriptors(properties));

    return properties;
}

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
    getExecutableByConfig: (config?: { name: string }) => Executable;
    getExecutables: () => Executable[];
    getBuilds: () => string[];
    getVersions: () => string[];
};
