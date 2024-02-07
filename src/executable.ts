import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";
import { Flavor } from "./flavor";
import { Constructor } from "@exabyte-io/code.js/dist/context";
import { FlavorData } from "./types";

const ExecutableBase = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);
type ExecutableBase = InstanceType<typeof ExecutableBase>;

export const Executable = ExecutableMixin(ExecutableBase);
export type Executable = InstanceType<typeof Executable>;

export function ExecutableMixin<
    T extends Constructor<ExecutableBase> = Constructor<ExecutableBase>,
>(superclass: T) {
    return class AdeExecutable extends superclass {
        static Flavor = Flavor;

        // @ts-ignore
        toJSON(exclude) {
            return super.toJSON(["flavors"].concat(exclude));
        }

        get flavorsTree() {
            return this.prop<Record<string, FlavorData>>("flavors");
        }

        get flavors() {
            return Object.keys(this.flavorsTree).map((key: string) => {
                return AdeExecutable.Flavor.create({
                    ...this.flavorsTree[key],
                    name: key,
                    executable: this,
                });
            });
        }

        get flavorsFromTree() {
            return Object.keys(this.flavorsTree).map((key: string) => {
                return new AdeExecutable.Flavor({ ...this.flavorsTree[key], name: key });
            });
        }

        get defaultFlavor() {
            return this.getFlavorByName();
        }

        getFlavorByName(name?: string | null) {
            return name ? this.getEntityByName(this.flavors, "flavor", name) as Flavor : undefined;
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
