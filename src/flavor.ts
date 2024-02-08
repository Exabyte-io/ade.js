import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";

import { Template } from "./template";
import { NamedTemplate } from "./types";
import { Constructor } from "@exabyte-io/code.js/dist/context";
import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";

const Base = RuntimeItemsMixin(NamedDefaultableHashedInMemoryEntity);
type FlavorBaseEntity = InstanceType<typeof Base>;

export function FlavorMixin<
    T extends Constructor<FlavorBaseEntity> = Constructor<FlavorBaseEntity>,
>(superclass: T) {
    return class Flavor extends superclass {
        get input(): NamedTemplate[] {
            return this.prop<NamedTemplate[]>("input", []);
        }

        get inputAsTemplates(): Template[] {
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

        getInputAsRenderedTemplates(context: object): AnyObject[] {
            return this.inputAsTemplates.map((t) => t.getRenderedJSON(context));
        }

        get disableRenderMaterials(): boolean {
            return this.prop<boolean>("isMultiMaterial", false);
        }
    }
}

export const Flavor = FlavorMixin(Base);

export type Flavor = typeof Flavor;
