"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.applicationMixin = applicationMixin;
exports.applicationStaticMixin = applicationStaticMixin;
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
const lodash_1 = __importDefault(require("lodash"));
const executable_1 = __importDefault(require("./executable"));
const tree_1 = require("./tree");
function applicationMixin(item) {
    // @ts-ignore
    const properties = {
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
            const tree = (0, application_flavors_js_1.getAppTree)(this.name);
            return Object.keys(tree)
                .filter((key) => {
                const { supportedApplicationVersions } = tree[key];
                return (!supportedApplicationVersions ||
                    supportedApplicationVersions.includes(this.version));
            })
                .map((key) => {
                return this.constructor.constructExecutable({
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
        getExecutableByName(name) {
            return this.constructor.constructExecutable((0, tree_1.getExecutableConfig)({
                appName: this.name,
                execName: name,
            }));
        },
        getExecutableByConfig(config) {
            return config ? this.getExecutableByName(config.name) : this.defaultExecutable;
        },
        getExecutables() {
            return this.executables;
        },
        getBuilds() {
            const data = (0, application_flavors_js_1.getAppData)(this.name);
            const { versions } = data;
            const builds = ["Default"];
            versions.map((v) => v.build && builds.push(v.build));
            return lodash_1.default.uniq(builds);
        },
        getVersions() {
            const data = (0, application_flavors_js_1.getAppData)(this.name);
            const { versions } = data;
            const these = versions.map((v) => v.version);
            return lodash_1.default.uniq(these);
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return item;
}
function applicationStaticMixin(Application) {
    const properties = {
        get defaultConfig() {
            return {
                name: "espresso",
                shortName: "qe",
                version: "6.3",
                summary: "Quantum Espresso",
                build: "Default",
            };
        },
        create(config) {
            return this.createFromNameVersionBuild(config);
        },
        createFromNameVersionBuild({ name, version = null, build = "Default", }) {
            return new Application({ name, version, build });
        },
        getUniqueAvailableNames() {
            return application_flavors_js_1.allApplications;
        },
        constructExecutable(config) {
            if (this.constructCustomExecutable) {
                return this.constructCustomExecutable(config);
            }
            return new executable_1.default(config);
        },
    };
    Object.defineProperties(Application, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
