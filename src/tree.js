import { allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@exabyte-io/code.js/utils";

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param Cls {*} optional class to use to create applications
 * @returns {Object} containing applications and applicationConfigs
 */
export function getAllApplications(Cls = null) {
    const applicationsTree = {};
    const applicationsArray = [];
    allApplications.map((appName) => {
        applicationsTree[appName] = {};
        const { versions, defaultVersion, build = "Default", ...appData } = getAppData(appName);
        applicationsTree[appName].defaultVersion = defaultVersion;
        versions.map((options) => {
            const { version } = options;
            if (!(version in applicationsTree[appName])) applicationsTree[appName][version] = {};
            const config = { ...appData, build, ...options };
            if (Cls) {
                applicationsTree[appName][version][build] = new Cls(config);
                applicationsArray.push(new Cls(config));
            } else {
                applicationsTree[appName][version][build] = config;
                applicationsArray.push(config);
            }
            return null;
        });
        return null;
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
    let current = version;
    if (!version) current = app.defaultVersion;
    return app[current][build];
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
    return getApplication({ applicationsTree, name, version, build });
}

/**
 * @summary Get executable config
 * @param appName {String} name of application to get executable for
 * @param execName {String|null} if not provided, find the executable with isDefault === true
 * @returns {*}
 */
export function getExecutableConfig({ appName, execName }) {
    const appTree = getAppTree(appName);
    Object.entries(appTree).map(([name, exec]) => {
        exec.name = name;
        return null;
    });
    if (!execName) return getOneMatchFromObject(appTree, "isDefault", true);
    return appTree[execName];
}
