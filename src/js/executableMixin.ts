import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import type { FlavorMixin } from "./flavorMixin";

type BaseFlavor = FlavorMixin & NamedInMemoryEntity & InMemoryEntity;
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;

export function executableMixin(item: Base) {
    // @ts-ignore
    const properties: ExecutableMixin & Base = {
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

export type BaseConstructor = Constructor<Base> & {
    constructCustomFlavor?: (config: object) => BaseFlavor;
};

export type ExecutableMixin = {
    toJSON: (exclude?: string[]) => object;
};
