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
/**
 * @summary Return all applications as both a nested object of Applications and an array of config objects
 * @param cls optional class to use to create applications
 * @returns containing applications and applicationConfigs
 */
export declare function getAllApplications(cls?: Constructor<ApplicationMixin> | null): {
    applicationsTree: Partial<Record<ApplicationName, LocalApplicationTreeItem>>;
    applicationsArray: ApplicationSchemaBase[];
};
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
