import { getAppTree } from "@exabyte-io/application-flavors.js";
import { getOneMatchFromObject } from "@mat3ra/code/dist/js/utils/object";
import type { ExecutableSchema } from "@mat3ra/esse/dist/js/types";

import Application from "./application";
import Executable from "./executable";
import Flavor from "./flavor";
import { type CreateApplicationConfig } from "./tree";

export class AdeFactory {
    // applications
    static createApplication({ name, version = null, build = "Default" }: CreateApplicationConfig) {
        return new Application({ name, version, build });
    }

    static getApplicationExecutables(application: Application) {
        const tree = getAppTree(application.name);

        return Object.keys(tree)
            .filter((key) => {
                const { supportedApplicationVersions } = tree[key];
                return (
                    !supportedApplicationVersions ||
                    supportedApplicationVersions.includes(application.version)
                );
            })
            .map((key) => new Executable({ ...tree[key], name: key }));
    }

    static getApplicationExecutableByName(application: Application, name?: string) {
        const appTree = getAppTree(application.name);

        Object.entries(appTree).forEach(([name, exec]) => {
            exec.name = name;
        });

        const config = name
            ? appTree[name]
            : (getOneMatchFromObject(appTree, "isDefault", true) as ExecutableSchema);

        return new Executable(config);
    }

    // TODO: remove this method and use getApplicationExecutableByName directly
    static getApplicationExecutableByConfig(application: Application, config?: { name: string }) {
        return this.getApplicationExecutableByName(application, config?.name);
    }

    // executables
    static getFlavorsByApplicationVersion(executable: Executable, version: string) {
        const filteredFlavors = this.getExecutableFlavors(executable).filter((flavor) => {
            const supportedApplicationVersions = flavor.prop<string[]>(
                "supportedApplicationVersions",
            );
            return !supportedApplicationVersions || supportedApplicationVersions.includes(version);
        });

        return filteredFlavors;
    }

    static getExecutableFlavors(executable: Executable) {
        const flavorsTree = executable.prop("flavors", {}) as Record<string, any>;

        return Object.keys(flavorsTree).map((key) => {
            return new Flavor({
                ...flavorsTree[key],
                name: key,
                executable,
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
}
