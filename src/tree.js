/* eslint-disable new-cap */
import { allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@mat3ra/code/dist/js/utils";

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls {*} optional class to use to create applications
 * @returns {Object} containing applications and applicationConfigs
 */
export function getAllApplications(cls = null) {
    const applicationsTree = {};
    const applicationsArray = [];
    allApplications.forEach((appName) => {
        applicationsTree[appName] = {};
        const { versions, defaultVersion, build = "Default", ...appData } = getAppData(appName);
        applicationsTree[appName].defaultVersion = defaultVersion;
        versions.forEach((options) => {
            const { version } = options;
            if (!(version in applicationsTree[appName])) applicationsTree[appName][version] = {};
            const config = { ...appData, build, ...options };
            if (cls) {
                applicationsTree[appName][version][build] = new cls(config);
                applicationsArray.push(new cls(config));
            } else {
                applicationsTree[appName][version][build] = config;
                applicationsArray.push(config);
            }
        });
    });
    return { applicationsTree, applicationsArray };
}

/**
 * @summary Get an application from the constructed applications
 * @param applicationsTree {Object} See getAllApplications applicationsTree object structure
 * @param name {String} name of the application
 * @param version {String|null} version of the application (optional, defaults to defaultVersion)
 * @param build {String} the build to use (optional, defaults to Default)
 * @return {*} an application
 */
export function getApplication({ applicationsTree, name, version = null, build = "Default" }) {
    const app = applicationsTree[name];
    const version_ = version || app.defaultVersion;
    if (!app[version_]) console.log(`Version ${version_} not available for ${name} !`);
    return app[version_]?.[build];
}

const { applicationsTree } = getAllApplications(null);

/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 * @param name
 * @param version {String|null}
 * @param build
 * @returns {*}
 */
export function getApplicationConfig({ name, version = null, build = "Default" }) {
    return getApplication({
        applicationsTree,
        name,
        version,
        build,
    });
}

/**
 * @summary Get executable config
 * @param appName {String} name of application to get executable for
 * @param execName {String|null} if not provided, find the executable with isDefault === true
 * @returns {*}
 */
export function getExecutableConfig({ appName, execName }) {
    const appTree = getAppTree(appName);
    Object.entries(appTree).forEach(([name, exec]) => {
        exec.name = name;
    });
    if (!execName) return getOneMatchFromObject(appTree, "isDefault", true);
    return appTree[execName];
}

/**
 * @summary Get flavor config
 * @param appName
 * @param execName
 * @param flavorName
 */
// eslint-disable-next-line no-unused-vars
export function getFlavorConfig({ appName, execName, flavorName }) {
    // TODO : reduce redundancy of object construction in getting flavors from executable
}
