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
const template_1 = __importDefault(require("./template"));
class ApplicationRegistry {
    static createApplication({ name, version = null, build = "Default" }) {
        const staticConfig = ApplicationRegistry.getApplicationConfig({ name, version, build });
        return new application_1.default({ ...staticConfig, name, version, build });
    }
    static getUniqueAvailableApplicationNames() {
        return application_flavors_js_1.allApplications;
    }
    /**
     * @summary Return all applications as both a nested object of Applications and an array of config objects
     * @returns containing applications and applicationConfigs
     */
    static getAllApplications() {
        if (this.applicationsTree && this.applicationsArray) {
            return {
                applicationsTree: this.applicationsTree,
                applicationsArray: this.applicationsArray,
            };
        }
        const applicationsTree = {};
        const applicationsArray = [];
        application_flavors_js_1.allApplications.forEach((appName) => {
            const { versions, defaultVersion, ...appData } = (0, application_flavors_js_1.getAppData)(appName);
            const appTreeItem = { defaultVersion };
            versions.forEach((versionInfo) => {
                const { version, build = "Default" } = versionInfo;
                const appVersion = version in appTreeItem && typeof appTreeItem[version] === "object"
                    ? appTreeItem[version]
                    : {};
                appTreeItem[version] = appVersion;
                const applicationConfig = {
                    ...appData,
                    build,
                    ...versionInfo,
                };
                appVersion[build] = applicationConfig;
                applicationsArray.push(applicationConfig);
            });
            applicationsTree[appName] = appTreeItem;
        });
        this.applicationsTree = applicationsTree;
        this.applicationsArray = applicationsArray;
        return {
            applicationsTree,
            applicationsArray: this.applicationsArray,
        };
    }
    /**
     * @summary Get an application from the constructed applications
     * @param name name of the application
     * @param version version of the application (optional, defaults to defaultVersion)
     * @param build  the build to use (optional, defaults to Default)
     * @return an application
     */
    static getApplicationConfig({ name, version = null, build = "Default", }) {
        var _a;
        const { applicationsTree } = this.getAllApplications();
        const app = applicationsTree[name];
        if (!app) {
            throw new Error(`Application ${name} not found`);
        }
        const version_ = version || app.defaultVersion;
        const appVersion = app[version_];
        if (!appVersion || typeof appVersion === "string") {
            console.warn(`Version ${version_} not available for ${name} !`);
            return null;
        }
        return (_a = appVersion[build]) !== null && _a !== void 0 ? _a : null;
    }
    static getExecutables({ name, version }) {
        const tree = (0, application_flavors_js_1.getAppTree)(name);
        return Object.keys(tree)
            .filter((key) => {
            const executable = tree[key];
            const { supportedApplicationVersions } = executable;
            return (!supportedApplicationVersions ||
                (version && supportedApplicationVersions.includes(version)));
        })
            .map((key) => new executable_1.default({ ...tree[key], name: key }));
    }
    static getExecutableByName(appName, execName) {
        const appTree = (0, application_flavors_js_1.getAppTree)(appName);
        Object.entries(appTree).forEach(([name, exec]) => {
            exec.name = name;
        });
        const config = execName
            ? appTree[execName]
            : (0, object_1.getOneMatchFromObject)(appTree, "isDefault", true);
        return new executable_1.default(config);
    }
    // TODO: remove this method and use getApplicationExecutableByName directly
    static getExecutableByConfig(appName, config) {
        return this.getExecutableByName(appName, config === null || config === void 0 ? void 0 : config.name);
    }
    static getExecutableFlavors(executable) {
        const flavorsTree = executable.prop("flavors", {});
        return Object.keys(flavorsTree).map((key) => {
            return new flavor_1.default({
                ...flavorsTree[key],
                name: key,
            });
        });
    }
    static getFlavorByName(executable, name) {
        return this.getExecutableFlavors(executable).find((flavor) => name ? flavor.name === name : flavor.isDefault);
    }
    static getFlavorByConfig(executable, config) {
        return this.getFlavorByName(executable, config === null || config === void 0 ? void 0 : config.name);
    }
    // flavors
    static getInputAsTemplates(flavor) {
        const appName = flavor.prop("applicationName", "");
        const execName = flavor.prop("executableName", "");
        return flavor.input.map((input) => {
            const inputName = input.templateName || input.name;
            const filtered = application_flavors_js_1.allTemplates.filter((temp) => temp.applicationName === appName &&
                temp.executableName === execName &&
                temp.name === inputName);
            if (filtered.length !== 1) {
                console.log(`found ${filtered.length} templates for app=${appName} exec=${execName} name=${inputName} expected 1`);
            }
            return new template_1.default({ ...filtered[0], name: input.name });
        });
    }
    static getInputAsRenderedTemplates(flavor, context) {
        return this.getInputAsTemplates(flavor).map((template) => {
            return template.getRenderedJSON(context);
        });
    }
    static getAllFlavorsForApplication(appName, version) {
        const allExecutables = this.getExecutables({ name: appName, version });
        return allExecutables.flatMap((executable) => this.getExecutableFlavors(executable));
    }
}
exports.default = ApplicationRegistry;
