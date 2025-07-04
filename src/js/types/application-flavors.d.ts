/**
 * TODO: @exabyte-io/application-flavors.js package must be removed and the source code must be moved to the current repo later in the future
 */
declare module "@exabyte-io/application-flavors.js" {
    import type { ExecutableSchema, FlavorSchema } from "@mat3ra/esse/dist/js/types";

    export const allApplications: ApplicationName[];

    export interface ApplicationTreeItem {
        supportedApplicationVersions?: string[];
        name: string;
        isDefault?: boolean;
    }

    export function getAppTree(
        name: ApplicationName,
    ): Record<string, ExecutableSchema & ApplicationTreeItem & { flavors: FlavorSchema[] }>;

    export interface ApplicationVersion {
        version: string;
        isDefault?: boolean;
        build?: string;
        hasAdvancedComputeOptions?: boolean;
    }

    export interface ApplicationData {
        name: string;
        shortName: string;
        summary: string;
        defaultVersion: string;
        isLicensed?: boolean;
        build?: string;
        versions: ApplicationVersion[];
    }

    export type ApplicationName =
        | "deepmd"
        | "espresso"
        | "exabyteml"
        | "jupyterLab"
        | "nwchem"
        | "python"
        | "shell"
        | "vasp";

    export function getAppData(appName: ApplicationName): ApplicationData;

    export type Template = {
        applicationName: string;
        executableName: string;
        name: string;
        content: string;
    };

    export const allTemplates: Template[];

    export const allowedResults: string[];

    export const allowedMonitors: string[];
}
