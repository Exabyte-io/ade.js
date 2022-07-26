import { mix } from "mixwith";
import { NamedDefaultableInMemoryEntity, RuntimeItemsMixin } from "@exabyte-io/code.js/dist/entity";
import { Flavor } from "./flavor";

export class Executable extends mix(NamedDefaultableInMemoryEntity).with(RuntimeItemsMixin) {
    static Flavor = Flavor;

    constructor(config) {
        super(config);
    }

    toJSON(exclude) {
        return super.toJSON(["flavors"].concat(exclude));
    }

    get flavorsTree() {
        return this.prop("flavors");
    }

    get flavors() {
        return Object.keys(this.flavorsTree).map(key => {
            return this.constructor.Flavor.create(
                Object.assign({}, this.flavorsTree[key], {name: key, executable: this})
            );
        });
    }

    get flavorsFromTree() {
        return Object.keys(this.flavorsTree).map(key => {
            return new this.constructor.Flavor(
                Object.assign({}, this.flavorsTree[key], {name: key})
            );
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

}
