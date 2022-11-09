import { NamedDefaultableInMemoryEntity, RuntimeItemsMixin } from "@exabyte-io/code.js/dist/entity";
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
        let flavor = this.getEntityByName(this.flavors, "flavor", name);
        if (!flavor) {
            console.warn(`Could not find flavor '${name}'! Using default instead.`);
            flavor = this.getEntityByName(this.flavors, "flavor", undefined); // extracts default flavor
        }
        return flavor;
    }

    getFlavorByConfig(config) {
        return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
    }
}
