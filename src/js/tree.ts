/* eslint-disable new-cap */
import {
    type ApplicationName,
    type ApplicationTreeItem,
    allApplications,
    getAppData,
    getAppTree,
} from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@mat3ra/code/dist/js/utils";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { ApplicationSchemaBase } from "@mat3ra/esse/dist/js/types";

import type { ApplicationMixin } from "./applicationMixin";

type ApplicationVersion = {
    [build: string]: ApplicationMixin | ApplicationSchemaBase;
};

type LocalApplicationTreeItem = {
    defaultVersion: string;
    [version: string]: ApplicationVersion | string;
};

type ApplicationTreeStructure = Partial<Record<ApplicationName, LocalApplicationTreeItem>>;

/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
export function getAllApplications(cls: Constructor<ApplicationMixin> | null = null) {
    const applicationsTree: ApplicationTreeStructure = {};
    const applicationsArray: ApplicationMixin | ApplicationSchemaBase[] = [];
    allApplications.forEach((appName) => {
        const { versions, defaultVersion, build = "Default", ...appData } = getAppData(appName);
        const appTreeItem: LocalApplicationTreeItem = { defaultVersion };

        versions.forEach((options) => {
            const { version } = options;

            const appVersion =
                version in appTreeItem && typeof appTreeItem[version] === "object"
                    ? appTreeItem[version]
                    : {};

            appTreeItem[version] = appVersion;

            const config = { ...appData, build, ...options };

            if (cls) {
                appVersion[build] = new cls(config);
                applicationsArray.push(new cls(config));
            } else {
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
export function getApplication({
    applicationsTree,
    name,
    version = null,
    build = "Default",
}: {
    applicationsTree: ApplicationTreeStructure;
    name: ApplicationName;
    version?: string | null;
    build?: string;
}) {
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
    return appVersion[build] ?? null;
}

const { applicationsTree } = getAllApplications(null);

export type CreateApplicationConfig = {
    name: ApplicationName;
    version?: string | null;
    build?: string;
};

/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 */
export function getApplicationConfig({
    name,
    version = null,
    build = "Default",
}: CreateApplicationConfig) {
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
export function getExecutableConfig({
    appName,
    execName,
}: {
    appName: ApplicationName;
    execName?: string | null;
}): ApplicationTreeItem {
    const appTree = getAppTree(appName);

    Object.entries(appTree).forEach(([name, exec]) => {
        exec.name = name;
    });

    if (!execName) {
        return getOneMatchFromObject(appTree, "isDefault", true) as ApplicationTreeItem;
    }

    return appTree[execName];
}
