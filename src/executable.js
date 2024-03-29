import { NamedDefaultableInMemoryEntity, RuntimeItemsMixin } from "@mat3ra/code/dist/js/entity";
import { mix } from "mixwith";

import { Flavor } from "./flavor";

export class Executable extends mix(NamedDefaultableInMemoryEntity).with(RuntimeItemsMixin) {
    static Flavor = Flavor;

    toJSON(exclude) {
        return super.toJSON(["flavors"].concat(exclude));
    }

    get flavorsTree() {
        return this.prop("flavors");
    }

    get flavors() {
        return Object.keys(this.flavorsTree).map((key) => {
            return this.constructor.Flavor.create({
                ...this.flavorsTree[key],
                name: key,
                executable: this,
            });
        });
    }

    get flavorsFromTree() {
        return Object.keys(this.flavorsTree).map((key) => {
            return new this.constructor.Flavor({ ...this.flavorsTree[key], name: key });
        });
    }

    get defaultFlavor() {
        return this.getFlavorByName();
    }

    getFlavorByName(name) {
        return this.getEntityByName(this.flavors, "flavor", name);
    }

    getFlavorByConfig(config) {
        return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
    }

    getFlavorsByApplicationVersion(version) {
        const filteredFlavors = this.flavors.filter((flavor) => {
            const supportedApplicationVersions = flavor.prop("supportedApplicationVersions");
            return !supportedApplicationVersions || supportedApplicationVersions.includes(version);
        });

        return filteredFlavors;
    }
}
