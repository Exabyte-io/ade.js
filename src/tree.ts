// @ts-expect-error application-flavors.js is not typed
import { allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@exabyte-io/code.js/dist/utils";
import { AppTree, ExecutableData } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
export function getAllApplications(cls?: Constructor<any>) {
    const applicationsTree = {};
    // @ts-expect-error application-flavors is not typed
    const applicationsArray = [];
    allApplications.forEach((appName: string) => {
        // @ts-ignore
        applicationsTree[appName] = {};
        const { versions, defaultVersion, build = "Default", ...appData } = getAppData(appName);
        // @ts-ignore
        applicationsTree[appName].defaultVersion = defaultVersion;
        // @ts-ignore
        versions.forEach((options) => {
            const { version } = options;
            // @ts-ignore
            if (!(version in applicationsTree[appName])) applicationsTree[appName][version] = {};
            const config = { ...appData, build, ...options };
            if (cls) {
            // @ts-ignore
                applicationsTree[appName][version][build] = new cls(config);
                applicationsArray.push(new cls(config));
            } else {
                // @ts-ignore
                applicationsTree[appName][version][build] = config;
                applicationsArray.push(config);
            }
        });
    });
    // @ts-ignore
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
// @ts-ignore applicationsTree is not typed
export function getApplication({ applicationsTree, name, version = null, build = "Default" }:{
    applicationsTree: object,
    name: string,
    version?: string,
    build?: string,
}) {
    // @ts-ignore applicationsTree is not typed
    const app = applicationsTree[name];
    const version_ = version || app.defaultVersion;
    if (!app[version_]) console.log(`Version ${version_} not available for ${name} !`);
    return app[version_]?.[build];
}

const { applicationsTree } = getAllApplications();

/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 * @param name {String}
 * @param version {String|undefined}
 * @param build {String}
 * @returns {*}
 */
export function getApplicationConfig({ name, version = undefined, build = "Default" }: {
    name: string;
    version?: string;
    build?: string;
}) {
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
 * @param execName if not provided, find the executable with isDefault === true
 * @returns an executable config
 */
export function getExecutableConfig({ appName, execName = undefined }: { appName: string, execName: string | undefined }):  ExecutableData | undefined {
    const appTree: AppTree = getAppTree(appName);
    Object.entries(appTree).forEach(([name, exec]) => {
        exec.name = name;
    });
    if (!execName) {
        console.log("No executable name provided, using default executable");
        return getOneMatchFromObject(appTree, "isDefault", true) as ExecutableData | undefined;
    }
    if (!appTree[execName]) {
        console.log(`Executable ${execName} not found for application ${appName}`);
        return undefined;
    }
    return appTree[execName];
}
