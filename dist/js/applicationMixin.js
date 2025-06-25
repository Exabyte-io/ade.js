"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.applicationMixin = applicationMixin;
exports.applicationStaticMixin = applicationStaticMixin;
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
function applicationMixin(item) {
    // @ts-expect-error
    const properties = {
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
        // create(config: CreateApplicationConfig) {
        //     return this.createFromNameVersionBuild(config);
        // },
        create({ name, version = null, build = "Default" }) {
            return new Application({ name, version, build });
        },
        getUniqueAvailableNames() {
            return application_flavors_js_1.allApplications;
        },
        // constructExecutable(this: BaseConstructor & typeof properties, config: object) {
        //     if (this.constructCustomExecutable) {
        //         return this.constructCustomExecutable(config);
        //     }
        //     return new Executable(config);
        // },
    };
    Object.defineProperties(Application, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
