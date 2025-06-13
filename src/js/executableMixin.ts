import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import { Flavor } from "./flavor";
import type { FlavorMixin } from "./flavorMixin";

type BaseFlavor = FlavorMixin & NamedInMemoryEntity & InMemoryEntity;
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;

export function executableMixin(item: Base) {
    // @ts-ignore
    const properties: ExecutableMixin & Base = {
        get defaultFlavor() {
            return this.getFlavorByName();
        },

        get flavorsTree() {
            return this.prop("flavors", {}) as Record<string, any>;
        },

        get flavors() {
            const tree = this.flavorsTree || {};

            return Object.keys(tree).map((key) => {
                return (this.constructor as unknown as ExecutableStaticProperties).constructFlavor({
                    ...tree[key],
                    name: key,
                    executable: this,
                });
            });
        },

        get flavorsFromTree() {
            return Object.keys(this.flavorsTree).map((key) => {
                return (this.constructor as unknown as ExecutableStaticProperties).constructFlavor({
                    ...this.flavorsTree[key],
                    name: key,
                });
            });
        },

        getFlavorByName(name?: string) {
            return this.getEntityByName(this.flavors, "flavor", name || "") as BaseFlavor;
        },

        getFlavorByConfig(config?: { name: string }) {
            return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        },

        getFlavorsByApplicationVersion(version: string) {
            const filteredFlavors = this.flavors.filter((flavor) => {
                const supportedApplicationVersions = flavor.prop<string[]>(
                    "supportedApplicationVersions",
                );
                return (
                    !supportedApplicationVersions || supportedApplicationVersions.includes(version)
                );
            });

            return filteredFlavors;
        },

        toJSON(exclude: string[] = []) {
            return super.toJSON(["flavors"].concat(exclude));
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return item;
}

export type CreateExecutableConfig = {
    name: string;
    flavors?: Record<string, any>;
};

export type BaseConstructor = Constructor<Base> & {
    constructCustomFlavor?: (config: object) => BaseFlavor;
};

export function executableStaticMixin<T extends BaseConstructor>(Executable: T) {
    const properties: ExecutableStaticProperties = {
        constructFlavor(this: BaseConstructor & typeof properties, config: object): BaseFlavor {
            if (this.constructCustomFlavor) {
                return this.constructCustomFlavor(config);
            }
            return new Flavor(config);
        },
    };

    Object.defineProperties(Executable, Object.getOwnPropertyDescriptors(properties));

    return properties;
}

export type ExecutableStaticProperties = {
    constructFlavor: (config: object) => BaseFlavor;
};

export type ExecutableMixin = {
    defaultFlavor: BaseFlavor;
    flavorsTree: Record<string, any>;
    flavors: BaseFlavor[];
    flavorsFromTree: BaseFlavor[];
    getFlavorByName: (name?: string) => BaseFlavor;
    getFlavorByConfig: (config?: { name: string }) => BaseFlavor;
    getFlavorsByApplicationVersion: (version: string) => BaseFlavor[];
    toJSON: (exclude?: string[]) => object;
};
