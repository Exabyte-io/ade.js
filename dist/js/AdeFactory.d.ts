import { type ApplicationName } from "@exabyte-io/application-flavors.js";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { ApplicationSchemaBase } from "@mat3ra/esse/dist/js/types";
import Application from "./application";
import type { ApplicationMixin } from "./applicationMixin";
import Executable from "./executable";
import Flavor from "./flavor";
import Template from "./template";
type ApplicationVersion = {
    [build: string]: ApplicationMixin | ApplicationSchemaBase;
};
type LocalApplicationTreeItem = {
    defaultVersion: string;
    [version: string]: ApplicationVersion | string;
};
export type CreateApplicationConfig = {
    name: ApplicationName;
    version?: string | null;
    build?: string;
};
type ApplicationTreeStructure = Partial<Record<ApplicationName, LocalApplicationTreeItem>>;
export default class AdeFactory {
    static applicationsTree: ApplicationTreeStructure;
    static applicationsArray: (ApplicationMixin | ApplicationSchemaBase)[];
    static createApplication({ name, version, build }: CreateApplicationConfig): Application;
    static getUniqueAvailableApplicationNames(): ApplicationName[];
    /**
     * @summary Return all applications as both a nested object of Applications and an array of config objects
     * @param Cls optional class to use to create applications
     * @returns containing applications and applicationConfigs
     */
    static getAllApplications(Cls?: Constructor<ApplicationMixin> | null): {
        applicationsTree: Partial<Record<ApplicationName, LocalApplicationTreeItem>>;
        applicationsArray: (ApplicationMixin | ApplicationSchemaBase)[];
    };
    /**
     * @summary Get an application from the constructed applications
     * @param name name of the application
     * @param version version of the application (optional, defaults to defaultVersion)
     * @param build  the build to use (optional, defaults to Default)
     * @return an application
     */
    static getApplicationConfig({ name, version, build, }: CreateApplicationConfig): ApplicationMixin | ApplicationSchemaBase | null;
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
