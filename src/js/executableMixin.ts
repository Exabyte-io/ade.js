import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

// import Flavor from "./flavor";
import type { FlavorMixin } from "./flavorMixin";

type BaseFlavor = FlavorMixin & NamedInMemoryEntity & InMemoryEntity;
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;

export function executableMixin(item: Base) {
    // @ts-ignore
    const properties: ExecutableMixin & Base = {
        // get defaultFlavor() {
        //     return this.getFlavorByName();
        // },

        // get flavorsTree() {
        //     return this.prop("flavors", {}) as Record<string, any>;
        // },

        //  get flavors() {
        //     return Object.keys(this.flavorsTree).map((key) => {
        //         return this.constructor.Flavor.create({
        //             ...this.flavorsTree[key],
        //             name: key,
        //             executable: this,
        //         });
        //     });
        // },

        // get flavorsFromTree() {
        //     return Object.keys(this.flavorsTree).map((key) => {
        //         return new this.constructor.Flavor({ ...this.flavorsTree[key], name: key });
        //     });
        // },

        // getFlavorByName(name?: string) {
        //     return this.getEntityByName(this.flavors, "flavor", name || "") as BaseFlavor;
        // },

        // getFlavorByConfig(config?: { name: string }) {
        //     return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        // },

        // getFlavorsByApplicationVersion(version: string) {
        //     const filteredFlavors = this.flavors.filter((flavor) => {
        //         const supportedApplicationVersions = flavor.prop<string[]>(
        //             "supportedApplicationVersions",
        //         );
        //         return (
        //             !supportedApplicationVersions || supportedApplicationVersions.includes(version)
        //         );
        //     });

        //     return filteredFlavors;
        // },

        toJSON(exclude: string[] = []) {
            const thisProto = Object.getPrototypeOf(this);
            const superProto = Object.getPrototypeOf(thisProto);
            const baseToJSON = superProto.toJSON;
            return baseToJSON.call(this, ["flavors"].concat(exclude));
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

export type ExecutableMixin = {
    toJSON: (exclude?: string[]) => object;
};
