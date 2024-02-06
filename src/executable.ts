import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";
import { Flavor } from "./flavor";

const RuntimeItemsEntity = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);
type ExecutableBaseEntity = InstanceType<typeof RuntimeItemsEntity>;

export type ExecutableBaseEntityConstructor<T extends ExecutableBaseEntity = ExecutableBaseEntity> = new (
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    ...args: any[]
) => T;

export function ExecutableMixin<
    T extends ExecutableBaseEntityConstructor = ExecutableBaseEntityConstructor,
>(superclass: T) {
    return class extends superclass {
        static Flavor = Flavor;

        toJSON(exclude) {
            return super.toJSON(["flavors"].concat(exclude));
        }

        get flavorsTree() {
            return this.prop<Flavor[]>("flavors");
        }

        get flavors() {
            return Object.keys(this.flavorsTree).map((key) => {
                return AdeExecutable.Flavor.create({
                    ...this.flavorsTree[key],
                    name: key,
                    executable: this,
                });
            });
        }

        get flavorsFromTree() {
            return Object.keys(this.flavorsTree).map((key) => {
                return new AdeExecutable.Flavor({ ...this.flavorsTree[key], name: key });
            });
        }

        get defaultFlavor() {
            return this.getFlavorByName();
        }

        getFlavorByName(name?: string | null) {
            return this.getEntityByName(this.flavors, "flavor", name) as Flavor;
        }

        getFlavorByConfig(config?: {name: string}) {
            return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        }

        getFlavorsByApplicationVersion(version) {
            const filteredFlavors = this.flavors.filter((flavor) => {
                const supportedApplicationVersions = flavor.prop<string[]>("supportedApplicationVersions");
                return !supportedApplicationVersions || supportedApplicationVersions.includes(version);
            });

            return filteredFlavors;
        }
    }
}

export const Executable = ExecutableMixin(
    RuntimeItemsMixin(
        NamedDefaultableHashedInMemoryEntity
    )
);

export type Executable = InstanceType<typeof Executable>;
