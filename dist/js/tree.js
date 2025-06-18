"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.getAllApplications = getAllApplications;
exports.getApplication = getApplication;
exports.getApplicationConfig = getApplicationConfig;
exports.getExecutableConfig = getExecutableConfig;
/* eslint-disable new-cap */
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
const utils_1 = require("@mat3ra/code/dist/js/utils");
/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
function getAllApplications(cls = null) {
    const applicationsTree = {};
    const applicationsArray = [];
    application_flavors_js_1.allApplications.forEach((appName) => {
        const { versions, defaultVersion, build = "Default", ...appData } = (0, application_flavors_js_1.getAppData)(appName);
        const appTreeItem = { defaultVersion };
        versions.forEach((options) => {
            const { version } = options;
            const appVersion = version in appTreeItem && typeof appTreeItem[version] === "object"
                ? appTreeItem[version]
                : {};
            appTreeItem[version] = appVersion;
            const config = { ...appData, build, ...options };
            if (cls) {
                appVersion[build] = new cls(config);
                applicationsArray.push(new cls(config));
            }
            else {
                appVersion[build] = config;
                applicationsArray.push(config);
            }
        });
        applicationsTree[appName] = appTreeItem;
    });
    return { applicationsTree, applicationsArray };
}
/**
 * @summary Get an application from the constructed applications
 * @param applicationsTree See getAllApplications applicationsTree object structure
 * @param name name of the application
 * @param version version of the application (optional, defaults to defaultVersion)
 * @param build  the build to use (optional, defaults to Default)
 * @return an application
 */
function getApplication({ applicationsTree, name, version = null, build = "Default", }) {
    var _a;
    const app = applicationsTree[name];
    if (!app) {
        throw new Error(`Application ${name} not found`);
    }
    const version_ = version || app.defaultVersion;
    const appVersion = app[version_];
    if (!appVersion || typeof appVersion === "string") {
        console.log(`Version ${version_} not available for ${name} !`);
    }
    if (typeof appVersion === "string") {
        return null;
    }
    return (_a = appVersion[build]) !== null && _a !== void 0 ? _a : null;
}
const { applicationsTree } = getAllApplications(null);
/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 */
function getApplicationConfig({ name, version = null, build = "Default", }) {
    return getApplication({
        applicationsTree,
        name,
        version,
        build,
    });
}
/**
 * @summary Get executable config
 * @param appName name of application to get executable for
 * @param execName  if not provided, find the executable with isDefault === true
 */
function getExecutableConfig({ appName, execName, }) {
    const appTree = (0, application_flavors_js_1.getAppTree)(appName);
    Object.entries(appTree).forEach(([name, exec]) => {
        exec.name = name;
    });
    if (!execName) {
        return (0, utils_1.getOneMatchFromObject)(appTree, "isDefault", true);
    }
    return appTree[execName];
}
