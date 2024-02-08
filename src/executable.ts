import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";

import { Flavor as AdeFlavor } from "./flavor";
import { FlavorData } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";
import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";

const Base = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);
type ExecutableBaseEntity = InstanceType<typeof Base>;

export function ExecutableMixin<
    F extends AdeFlavor = AdeFlavor,
    T extends Constructor<ExecutableBaseEntity> = Constructor<ExecutableBaseEntity>,
>(superclass: T, Flavor: F) {
    return class Executable extends superclass {
        static Flavor = Flavor;

        constructor(...config: any[]) {
            super(...config);
        }

        // @ts-ignore
        toJSON(exclude?: string[]): AnyObject {
            return super.toJSON(["flavors"].concat(exclude ? exclude : []));
        }

        get flavorsTree(): Record<string, FlavorData> {
            return this.prop<Record<string, FlavorData>>("flavors");
        }

        get flavors() {
            return Object.keys(this.flavorsTree).map((key: string) => {
                return Executable.Flavor.create({
                    ...this.flavorsTree[key] as FlavorData,
                    name: key,
                    executable: this as Executable,
                });
            });
        }

        get flavorsFromTree() {
            return Object.keys(this.flavorsTree).map((key: string) => {
                return new Executable.Flavor({ ...this.flavorsTree[key], name: key });
            });
        }

        get defaultFlavor() {
            return this.getFlavorByName();
        }

        getFlavorByName(name?: string | null) {
            return name ? this.getEntityByName(this.flavors, "flavor", name) : undefined;
        }

        getFlavorByConfig(config?: {name: string}) {
            return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        }

        getFlavorsByApplicationVersion(version: string) {
            const filteredFlavors = this.flavors.filter((flavor) => {
                const supportedApplicationVersions = flavor.prop<string[]>("supportedApplicationVersions");
                return !supportedApplicationVersions || supportedApplicationVersions.includes(version);
            });

            return filteredFlavors;
        }
    }
}

export const Executable = ExecutableMixin(Base, AdeFlavor);

export type Executable = InstanceType<typeof Executable>;