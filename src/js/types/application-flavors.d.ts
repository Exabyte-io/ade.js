declare module "@exabyte-io/application-flavors.js" {
    export const allApplications: string[];

    export interface ApplicationTreeItem {
        supportedApplicationVersions?: string[];
        name: string;
        isDefault?: boolean;
    }

    export function getAppTree(name: string): Record<string, ApplicationTreeItem>;

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
}
