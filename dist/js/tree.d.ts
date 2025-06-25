import { type ApplicationName } from "@exabyte-io/application-flavors.js";
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
export declare function getAllApplications(cls?: Constructor<ApplicationMixin> | null): {
    applicationsTree: Partial<Record<ApplicationName, LocalApplicationTreeItem>>;
    applicationsArray: ApplicationSchemaBase[];
};
/**
 * @summary Get an application from the constructed applications
 * @param applicationsTree See getAllApplications applicationsTree object structure
 * @param name name of the application
 * @param version version of the application (optional, defaults to defaultVersion)
 * @param build  the build to use (optional, defaults to Default)
 * @return an application
 */
export declare function getApplication({ applicationsTree, name, version, build, }: {
    applicationsTree: ApplicationTreeStructure;
    name: ApplicationName;
    version?: string | null;
    build?: string;
}): ApplicationMixin | ApplicationSchemaBase | null;
export type CreateApplicationConfig = {
    name: ApplicationName;
    version?: string | null;
    build?: string;
};
/**
 * @summary Get pre-defined application config from an already generated applicationsTree of configs
 */
export declare function getApplicationConfig({ name, version, build, }: CreateApplicationConfig): ApplicationMixin | ApplicationSchemaBase | null;
export {};
