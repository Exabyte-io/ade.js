import {
    type ApplicationName,
    allApplications,
    allTemplates,
    getAppData,
    getAppTree,
} from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@mat3ra/code/dist/js/utils/object";
import type { ApplicationSchemaBase, ExecutableSchema } from "@mat3ra/esse/dist/js/types";

import Application from "./application";
import Executable from "./executable";
import Flavor from "./flavor";
import Template from "./template";

type ApplicationVersion = {
    [build: string]: ApplicationSchemaBase;
};

type ApplicationTreeItem = {
    defaultVersion: string;
    [version: string]: ApplicationVersion | string;
};

export type CreateApplicationConfig = {
    name: ApplicationName;
    version?: string | null;
    build?: string;
};

type ApplicationTree = Partial<Record<ApplicationName, ApplicationTreeItem>>;

export default class ApplicationRegistry {
    // applications
    static applicationsTree?: ApplicationTree;

    static applicationsArray?: ApplicationSchemaBase[];

    static createApplication({ name, version = null, build = "Default" }: CreateApplicationConfig) {
        const staticConfig = ApplicationRegistry.getApplicationConfig({ name, version, build });
        return new Application({ ...staticConfig, name, version, build });
    }

    static getUniqueAvailableApplicationNames() {
        return allApplications;
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

        const applicationsTree: ApplicationTree = {};
        const applicationsArray: ApplicationSchemaBase[] = [];

        allApplications.forEach((appName) => {
            const { versions, defaultVersion, ...appData } = getAppData(appName);
            const appTreeItem: ApplicationTreeItem = { defaultVersion };

            versions.forEach((versionInfo) => {
                const { version, build = "Default" } = versionInfo;

                const appVersion =
                    version in appTreeItem && typeof appTreeItem[version] === "object"
                        ? appTreeItem[version]
                        : {};

                appTreeItem[version] = appVersion;

                const applicationConfig: ApplicationSchemaBase = {
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
    static getApplicationConfig({
        name,
        version = null,
        build = "Default",
    }: CreateApplicationConfig) {
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

        return appVersion[build] ?? null;
    }

    static getExecutables({ name, version }: { name: ApplicationName; version?: string }) {
        const tree = getAppTree(name);

        return Object.keys(tree)
            .filter((key) => {
                const executable = tree[key];
                const { supportedApplicationVersions } = executable;
                return (
                    !supportedApplicationVersions ||
                    (version && supportedApplicationVersions.includes(version))
                );
            })
            .map((key) => new Executable({ ...tree[key], name: key }));
    }

    static getExecutableByName(appName: ApplicationName, execName?: string) {
        const appTree = getAppTree(appName);

        Object.entries(appTree).forEach(([name, exec]) => {
            exec.name = name;
        });

        const config = execName
            ? appTree[execName]
            : (getOneMatchFromObject(appTree, "isDefault", true) as ExecutableSchema);

        return new Executable(config);
    }

    // TODO: remove this method and use getApplicationExecutableByName directly
    static getExecutableByConfig(appName: ApplicationName, config?: { name: string }) {
        return this.getExecutableByName(appName, config?.name);
    }

    static getExecutableFlavors(executable: Executable) {
        const flavorsTree = executable.prop("flavors", {}) as Record<string, any>;

        return Object.keys(flavorsTree).map((key) => {
            return new Flavor({
                ...flavorsTree[key],
                name: key,
            });
        });
    }

    static getFlavorByName(executable: Executable, name?: string) {
        return this.getExecutableFlavors(executable).find((flavor) =>
            name ? flavor.name === name : flavor.isDefault,
        );
    }

    static getFlavorByConfig(executable: Executable, config?: { name: string }) {
        return this.getFlavorByName(executable, config?.name);
    }

    // flavors
    static getInputAsTemplates(flavor: Flavor) {
        const appName = flavor.prop("applicationName", "");
        const execName = flavor.prop("executableName", "");

        return flavor.input.map((input) => {
            const inputName = input.templateName || input.name;

            const filtered = allTemplates.filter(
                (temp) =>
                    temp.applicationName === appName &&
                    temp.executableName === execName &&
                    temp.name === inputName,
            );

            if (filtered.length !== 1) {
                console.log(
                    `found ${filtered.length} templates for app=${appName} exec=${execName} name=${inputName} expected 1`,
                );
            }

            return new Template({ ...filtered[0], name: input.name });
        });
    }

    static getInputAsRenderedTemplates(flavor: Flavor, context: Record<string, unknown>) {
        return this.getInputAsTemplates(flavor).map((template) => {
            return template.getRenderedJSON(context);
        });
    }

    static getAllFlavorsForApplication(appName: ApplicationName, version?: string) {
        const allExecutables = this.getExecutables({ name: appName, version });

        return allExecutables.flatMap((executable) => this.getExecutableFlavors(executable));
    }
}
