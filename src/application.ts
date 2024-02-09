// @ts-expect-error application-flavors.js is not typed
import { allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { NamedDefaultableInMemoryEntity } from "@exabyte-io/code.js/dist/entity";

import lodash from "lodash";

import { Executable as AdeExecutable } from "./executable";
import { getApplicationConfig, getExecutableConfig } from "./tree";
import { ApplicationConfig, ApplicationData } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";

const Base = NamedDefaultableInMemoryEntity;
abstract class ApplicationBaseEntity extends Base {};

export function ApplicationMixin<
    E extends Constructor<AdeExecutable> = Constructor<AdeExecutable>,
    T extends Constructor<ApplicationBaseEntity> = Constructor<ApplicationBaseEntity>,
>(superclass: T, Executable: E) {
    return class Application extends superclass {
        static Executable = Executable;

        constructor(...args: any[]) {
            const config = args[0] as ApplicationConfig;
            if (!config || typeof config.name !== "string") throw new Error("Invalid application configuration object.");
            const staticConfig = getApplicationConfig(config);
            super({ ...staticConfig, ...config });
        }

        // TODO: extract this from application-flavors "global" default config for espresso 5.4.0
        static get defaultConfig() {
            return {
                name: "espresso",
                shortName: "qe",
                version: "6.3",
                summary: "Quantum Espresso",
                build: "Default",
            };
        }

        static create(config: {
            name: string,
            version?: string,
            build?: string
        }) {
            return this.createFromNameVersionBuild(config);
        }

        static createFromNameVersionBuild({
            name,
            version = undefined,
            build = "Default"
        }: {
            name: string,
            version?: string,
            build?: string
        }) {
            return new Application({ name, version, build });
        }

        getExecutables() {
            return this.executables;
        }

        getBuilds() {
            const data = getAppData(this.prop("name")) as ApplicationData;
            const { versions } = data;
            const builds = ["Default"];
            versions.map((v) => v.build && builds.push(v.build));
            return lodash.uniq(builds);
        }

        getVersions() {
            const data = getAppData(this.prop("name")) as ApplicationData;
            const { versions } = data;
            const these: string[] = versions.map((v) => v.version);
            return lodash.uniq(these);
        }

        static getUniqueAvailableNames() {
            return allApplications;
        }

        getExecutableByName(name?: string) {
            return new Application.Executable(
                getExecutableConfig({
                    appName: this.prop("name"),
                    execName: name,
                }),
            );
        }

        getExecutableByConfig(config: {name: string} | null | undefined = null) {
            return config ? this.getExecutableByName(config.name) : this.defaultExecutable;
        }

        get defaultExecutable() {
            return this.getExecutableByName();
        }

        // override upon inheritance
        // eslint-disable-next-line class-methods-use-this
        get allowedModelTypes() {
            return [];
        }

        get summary() {
            return this.prop<string>("summary");
        }

        get version() {
            return this.prop<string>("version");
        }

        get build() {
            return this.prop<string>("build");
        }

        get shortName() {
            return this.prop<string>("shortName", this.prop<string>("name"));
        }

        get executables() {
            const tree = getAppTree(this.prop("name"));
            return Object.keys(tree)
                .filter((key) => {
                    const { supportedApplicationVersions } = tree[key];
                    return (
                        !supportedApplicationVersions ||
                        supportedApplicationVersions.includes(this.prop("version"))
                    );
                })
                .map((key) => {
                    return new Application.Executable({ ...tree[key], name: key });
                });
        }

        get hasAdvancedComputeOptions() {
            return this.prop<boolean>("hasAdvancedComputeOptions", false);
        }

        get isLicensed() {
            return this.prop<boolean>("isLicensed", false);
        }

        get isUsingMaterial() {
            const materialUsingApplications = ["vasp", "nwchem", "espresso", "exabyteml"];
            return materialUsingApplications.includes(this.name);
        }
    }
}

export const Application = ApplicationMixin(Base, AdeExecutable);

export type Application = typeof Application;
