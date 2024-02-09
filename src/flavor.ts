import {
    NamedDefaultableInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";

import { Template as AdeTemplate } from "./template";
import { Constructor } from "@exabyte-io/code.js/dist/context";
import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";
import { ExecutionUnitInputItemSchemaForPhysicsBasedSimulationEngines, TemplateSchema1 } from "@exabyte-io/code.js/dist/types";

const Base = RuntimeItemsMixin(NamedDefaultableInMemoryEntity);
type FlavorBaseEntity = InstanceType<typeof Base>;

export function FlavorMixin<
    S extends AdeTemplate = AdeTemplate,
    T extends Constructor<FlavorBaseEntity> = Constructor<FlavorBaseEntity>,
>(superclass: T, Template: S) {
    return class Flavor extends superclass {
        static Template = Template;

        get input() {
            return this.prop<(ExecutionUnitInputItemSchemaForPhysicsBasedSimulationEngines & {templateName?: string})[]>("input", []);
        }

        get inputAsTemplates() {
            return this.input.map((input) => {
                const templateName = 'templateName' in input && input.templateName ? input.templateName : input.name;
                const template = Flavor.Template.fromFlavor(
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

export const Flavor = FlavorMixin(Base, AdeTemplate);

export type Flavor = typeof Flavor;
