// TODO: remove these types once application-flavors is moved to TypeScript
export type VersionData = {
    version: string;
    build: string;
    hasAdvancedComputeOptions: boolean;
    isDefault: boolean;
}

export type ApplicationData = {
    name: string;
    shortName: string;
    summary: string;
    defaultVersion: string;
    versions: VersionData[];
}

export type ApplicationConfig = Omit<ApplicationData, "versions" | "defaultVersion"> & VersionData;

export type AppTree = {
    [applicationName: string]: ExecutableData;
}

export type ExecutableData = {
    name?: string;
    isDefault: boolean;
    hasAdvancedComputeOptions?: boolean;
    results: string[];
    monitors?: string[];
    preProcessors?: string[];
    postProcessors?: string[];
    flavors: {[flavorName: string]: FlavorData};
}

export type FlavorData = {
    isDefault: boolean;
    monitors?: string[];
    results?: string[];
    preProcessors?: string[];
    postProcessors?: string[];
    input: TemplateData[];
    applicationName: string;
    executableName: string;
}

export type TemplateData = {
    name: string;
    templateName: string;
    applicationName?: string;
    executableName?: string;
}
