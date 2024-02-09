import { AllowedApplications, allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@exabyte-io/code.js/dist/utils";
import { ExecutableData } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";
import { ApplicationData, ApplicationTree } from "@exabyte-io/application-flavors.js/lib/js/build_templates";

type AppTree = {
    // application
    [appName in AllowedApplications]: {
        // version
        defaultVersion: string
    } & {
        [version: string]: {
            // build
            [build: string]: {
                build: string;
                name: string;
                summary: string;
                shortName: string;
            }
        }
    };
}

type AppArray = {
    version: string;
    build: string;
    name: string;
    summary: string;
    shortName: string;
    isDefault: boolean;
}[];

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
export function getAllApplications(cls?: Constructor<any>) {
    const applicationsTree: AppTree = {} as AppTree;
    const applicationsArray: AppArray = [];
    allApplications.forEach((appName) => {
        const { versions, defaultVersion, ...appData } = getAppData(appName);
        // @ts-ignore
        applicationsTree[appName] = { defaultVersion };
        versions.forEach((options) => {
            const { version, build = "Default" } = options;
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
export function getApplication({ applicationsTree, name, version, build = "Default" }:{
    applicationsTree: AppTree,
    name: AllowedApplications,
    version?: string,
    build?: string,
}) {
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
export function getApplicationConfig({ name, version, build = "Default" }: {
    name: AllowedApplications;
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
export function getExecutableConfig({ appName, execName }: { appName: AllowedApplications, execName?: string }) {
    const appTree = getAppTree(appName);
    console.log(appTree)
    if (!execName) {
        console.log("No executable name provided, using default executable");
        // iterate over properties of appTree and find the one with isDefault === true
        Object.entries(appTree).forEach(([key, value]) => {
            if (value.isDefault) {
                console.log(`Found default executable: ${key}`);
                execName = key;
            }
        });
        if (!execName) {
            console.log(`No default executable found for application ${appName}`);
            return undefined;
        }
    }
    if (!appTree[execName]) {
        console.log(`Executable ${execName} not found for application ${appName}`);
        return undefined;
    }
    return {...appTree[execName], name: execName};
}
