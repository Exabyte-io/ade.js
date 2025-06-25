"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
const object_1 = require("@mat3ra/code/dist/js/utils/object");
const application_1 = __importDefault(require("./application"));
const executable_1 = __importDefault(require("./executable"));
const flavor_1 = __importDefault(require("./flavor"));
class AdeFactory {
    // applications
    static createApplication({ name, version = null, build = "Default" }) {
        return new application_1.default({ name, version, build });
    }
    static getApplicationExecutables(application) {
        const tree = (0, application_flavors_js_1.getAppTree)(application.name);
        return Object.keys(tree)
            .filter((key) => {
            const { supportedApplicationVersions } = tree[key];
            return (!supportedApplicationVersions ||
                supportedApplicationVersions.includes(application.version));
        })
            .map((key) => new executable_1.default({ ...tree[key], name: key }));
    }
    static getApplicationExecutableByName(application, name) {
        const appTree = (0, application_flavors_js_1.getAppTree)(application.name);
        Object.entries(appTree).forEach(([name, exec]) => {
            exec.name = name;
        });
        const config = name
            ? appTree[name]
            : (0, object_1.getOneMatchFromObject)(appTree, "isDefault", true);
        return new executable_1.default(config);
    }
    // TODO: remove this method and use getApplicationExecutableByName directly
    static getApplicationExecutableByConfig(application, config) {
        return this.getApplicationExecutableByName(application, config === null || config === void 0 ? void 0 : config.name);
    }
    // executables
    static getFlavorsByApplicationVersion(executable, version) {
        const filteredFlavors = this.getExecutableFlavors(executable).filter((flavor) => {
            const supportedApplicationVersions = flavor.prop("supportedApplicationVersions");
            return !supportedApplicationVersions || supportedApplicationVersions.includes(version);
        });
        return filteredFlavors;
    }
    static getExecutableFlavors(executable) {
        const flavorsTree = executable.prop("flavors", {});
        return Object.keys(flavorsTree).map((key) => {
            return new flavor_1.default({
                ...flavorsTree[key],
                name: key,
                executable,
            });
        });
    }
    static getFlavorByName(executable, name) {
        return this.getExecutableFlavors(executable).find((flavor) => name ? flavor.name === name : flavor.isDefault);
    }
    static getFlavorByConfig(executable, config) {
        return this.getFlavorByName(executable, config === null || config === void 0 ? void 0 : config.name);
    }
}
exports.default = AdeFactory;
