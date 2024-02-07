// @ts-ignore
import { allApplications, getAppData, getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@exabyte-io/code.js/dist/utils";
import { AppTree, ApplicationArray, ApplicationConfig, ApplicationData, ApplicationTree, ExecutableData, VersionData } from "./types";

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
export function getAllApplications(cls: InstanceType<any> | null) {
    const applicationsTree: ApplicationTree = {};
    const applicationsArray: ApplicationArray = [];
    allApplications.forEach((appName: string) => {
        // applicationsTree[appName] = {};
        const { versions, defaultVersion, ...appData } = getAppData(appName) as ApplicationData;
        applicationsTree[appName].defaultVersion = defaultVersion;
        versions.forEach(({ version, build, ...versionData }) => {
            // add version to applicationsTree if it doesn't exist
            if (!(version in applicationsTree[appName])) applicationsTree[appName][version] = {};
            // convert to class instance if cls is provided otherwise use the config
            const config: ApplicationConfig = { ...appData, version, build, ...versionData };
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
 * @param applicationsTree See getAllApplications applicationsTree object structure
 * @param name name of the application
 * @param version version of the application (optional, defaults to defaultVersion)
 * @param build the build to use (optional, defaults to Default)
 * @return an application
 */
export function getApplication({
    applicationsTree,
    name,
    version = undefined,
    build = "Default"
}: {
    applicationsTree: ApplicationTree,
    name: string
    version?: string
    build?: string
}): ApplicationConfig | undefined {
    const app = applicationsTree[name];
    if (!app) {
        console.log(`Application ${name} not found!`);
        return undefined;
    }
    const version_ = version || app.defaultVersion;
    const appData = app[version_][build];
    if (!appData) {
        console.log(`Version ${version_} not available for ${name}!`);
        return undefined;
    }
    return {
        ...app,
        ...appData,
        version: appData.version, // Explicitly set version and build in case they are different in app and versionData
        build: appData.build || build
    };
}

const { applicationsTree } = getAllApplications(null);

/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 * @param name name of the application
 * @param version version of the application (optional, defaults to defaultVersion)
 * @param build the build to use (optional, defaults to Default)
 * @returns an application config
 */
export function getApplicationConfig({
    name,
    version = undefined,
    build = "Default"
}: {
    name: string,
    version?: string,
    build?: string
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
