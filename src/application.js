import lodash from "lodash";
import { getAppTree, getAppData, allApplications } from "@exabyte-io/application-flavors.js";
import { NamedDefaultableInMemoryEntity } from "@exabyte-io/code.js/dist/entity";

import { Executable } from "./executable";
import { getApplicationConfig, getExecutableConfig } from "./tree";


export class Application extends NamedDefaultableInMemoryEntity {
    static Executable = Executable;

    constructor(config) {
        const staticConfig = getApplicationConfig(config);
        super({...staticConfig, ...config});
    }

    // TODO: extract this from application-flavors "global" default config for espresso 5.4.0
    static get defaultConfig() {
        return {
            name: 'espresso',
            shortName: 'qe',
            version: '5.4.0',
            summary: 'Quantum Espresso',
            build: 'Default',
        }
    }

    static create(config) {
        return this.createFromNameVersionBuild(config);
    }

    static createFromNameVersionBuild({name, version = null, build = "Default"}) {
        return new Application({name, version, build});
    }

    getExecutables() {
        return this.executables;
    }

    getBuilds() {
        const data = getAppData(this.prop("name"));
        const {versions} = data;
        const builds = ["Default"];
        versions.map((v) => v.build && builds.push(v.build));
        return lodash.uniq(builds);
    }

    getVersions() {
        const data = getAppData(this.prop("name"));
        const {versions} = data;
        const these = versions.map((v) => v.version);
        return lodash.uniq(these);
    }

    static getUniqueAvailableNames() {
        return allApplications;
    }

    getExecutableByName(name = null) {
        return new this.constructor.Executable(
            getExecutableConfig({
                appName: this.prop("name"),
                execName: name,
            })
        );
    }

    getExecutableByConfig(config) {
        return config ? this.getExecutableByName(config.name) : this.defaultExecutable;
    }

    get defaultExecutable() {
        return this.getExecutableByName();
    }

    // override upon inheritance
    get allowedModelTypes() {
        return []
    }

    get summary() {
        return this.prop("summary");
    }

    get version() {
        return this.prop("version");
    }

    get build() {
        return this.prop("build");
    }

    get shortName() {
        return this.prop("shortName", this.prop("name"));
    }

    get executables() {
        const tree = getAppTree(this.prop("name"));
        return Object.keys(tree).map(key => {
            return new this.constructor.Executable(
                Object.assign({}, tree[key], {name: key})
            )
        })
    }

    get hasAdvancedComputeOptions() {
        return this.prop("hasAdvancedComputeOptions");
    }

    get isLicensed() {
        return this.prop("isLicensed");
    }

}