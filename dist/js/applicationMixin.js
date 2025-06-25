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
        create({ name, version = null, build = "Default" }) {
            return new Application({ name, version, build });
        },
        getUniqueAvailableNames() {
            return application_flavors_js_1.allApplications;
        },
    };
    Object.defineProperties(Application, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
