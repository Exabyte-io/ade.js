import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";

import { Template } from "./template";

const RuntimeItemsEntity = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);

type FlavorBaseEntity = InstanceType<typeof RuntimeItemsEntity>

export type FlavorBaseEntityConstructor<T extends FlavorBaseEntity = FlavorBaseEntity> = new (
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    ...args: any[]
) => T;

type NamedTemplate = {name: string} | {templateName: string};

export function FlavorMixin<
    T extends FlavorBaseEntityConstructor = FlavorBaseEntityConstructor,
>(superclass: T) {
    return class extends superclass {
        get input() {
            return this.prop<NamedTemplate[]>("input", []);
        }

        get inputAsTemplates() {
            return this.input.map((input) => {
                const templateName = 'name' in input ? input.name : input.templateName;
                const template = Template.fromFlavor(
                    this.prop("applicationName"),
                    this.prop("executableName"),
                    templateName,
                );
                // Override template name upon creation.
                // Used to, for example, set "KPOINTS" as name of the input file while "KPOINTS_BANDS" as the template name.
                template.name = templateName
                return template;
            });
        }

        getInputAsRenderedTemplates(context) {
            return this.inputAsTemplates.map((t) => t.getRenderedJSON(context));
        }

        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        }
    }
}

export const Flavor = FlavorMixin(
    RuntimeItemsMixin(
        NamedDefaultableHashedInMemoryEntity
    )
);

export type Flavor = InstanceType<typeof Flavor>;
