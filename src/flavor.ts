import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";

import { Template } from "./template";
import { NamedTemplate } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";

const Base = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);
type FlavorBaseEntity = InstanceType<typeof Base>;

export function FlavorMixin<
    T extends Constructor<FlavorBaseEntity> = Constructor<FlavorBaseEntity>,
>(superclass: T) {
    return class Flavor extends RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity) {
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

        getInputAsRenderedTemplates(context: object) {
            return this.inputAsTemplates.map((t) => t.getRenderedJSON(context));
        }

        get disableRenderMaterials() {
            return this.prop("isMultiMaterial", false);
        }
    }
}

export const Flavor = FlavorMixin(Base);

export type Flavor = InstanceType<typeof Flavor>;
