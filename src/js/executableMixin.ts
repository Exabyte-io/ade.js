import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { DefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import type { ExecutableSchema } from "@mat3ra/esse/dist/js/types";

import type { FlavorMixin } from "./flavorMixin";

type BaseFlavor = FlavorMixin & NamedInMemoryEntity & InMemoryEntity;
type Base = InMemoryEntity & NamedInMemoryEntity & DefaultableInMemoryEntity;

export function executableMixin(item: Base) {
    // @ts-expect-error
    const properties: ExecutableMixin & Base = {
        get applicationId() {
            return this.prop("applicationId", []);
        },
        set applicationId(value: string[]) {
            this.setProp("applicationId", value);
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}

export function executableStaticMixin(Executable: Constructor<Base>) {
    const properties: ExecutableStaticMixin = {
        get jsonSchema() {
            return JSONSchemasInterface.getSchemaById("software/executable") as ExecutableSchema;
        },
    };

    Object.defineProperties(Executable, Object.getOwnPropertyDescriptors(properties));
}

export type BaseConstructor = Constructor<Base> & {
    constructCustomFlavor?: (config: object) => BaseFlavor;
};

export type ExecutableMixin = {
    applicationId: string[];
    toJSON: () => ExecutableSchema & AnyObject;
};

export type ExecutableStaticMixin = {
    jsonSchema: ExecutableSchema;
};
