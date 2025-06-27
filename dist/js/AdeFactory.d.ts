import { type ApplicationName } from "@exabyte-io/application-flavors.js";
import type { ApplicationSchemaBase } from "@mat3ra/esse/dist/js/types";
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
export default class AdeFactory {
    static applicationsTree: ApplicationTree;
    static applicationsArray: ApplicationSchemaBase[];
    static createApplication({ name, version, build }: CreateApplicationConfig): Application;
    static getUniqueAvailableApplicationNames(): ApplicationName[];
    /**
     * @summary Return all applications as both a nested object of Applications and an array of config objects
     * @returns containing applications and applicationConfigs
     */
    static getAllApplications(): {
        applicationsTree: Partial<Record<ApplicationName, ApplicationTreeItem>>;
        applicationsArray: ApplicationSchemaBase[];
    };
    /**
     * @summary Get an application from the constructed applications
     * @param name name of the application
     * @param version version of the application (optional, defaults to defaultVersion)
     * @param build  the build to use (optional, defaults to Default)
     * @return an application
     */
    static getApplicationConfig({ name, version, build, }: CreateApplicationConfig): ApplicationSchemaBase | null;
    static getExecutables(application: Application): Executable[];
    static getExecutableByName(application: Application, name?: string): Executable;
    static getExecutableByConfig(application: Application, config?: {
        name: string;
    }): Executable;
    static getFlavorsByApplicationVersion(executable: Executable, version: string): Flavor[];
    static getExecutableFlavors(executable: Executable): Flavor[];
    static getFlavorByName(executable: Executable, name?: string): Flavor | undefined;
    static getFlavorByConfig(executable: Executable, config?: {
        name: string;
    }): Flavor | undefined;
    static getInputAsTemplates(flavor: Flavor): Template[];
    static getInputAsRenderedTemplates(flavor: Flavor, context: Record<string, unknown>): import("@mat3ra/esse/dist/js/esse/types").AnyObject[];
}
export {};
