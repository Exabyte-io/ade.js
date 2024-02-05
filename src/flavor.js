import {
    NamedDefaultableHashedInMemoryEntity,
    RuntimeItemsMixin,
} from "@exabyte-io/code.js/dist/entity";
import { mix } from "mixwith";

import { Template } from "./template";

export class Flavor extends mix(NamedDefaultableHashedInMemoryEntity).with(RuntimeItemsMixin) {
    get input() {
        return this.prop("input", []);
    }

    // TODO : prevent this from running in client
    get inputAsTemplates() {
        // get template from application-flavors rather than lookup by ID
        return this.input.map((input) => {
            const template = Template.fromFlavor(
                this.prop("applicationName"),
                this.prop("executableName"),
                input.templateName || input.name,
            );
            // Override template name upon creation.
            // Used to, for example, set "KPOINTS" as name of the input file while "KPOINTS_BANDS" as the template name.
            template.name = input.name;
            return template;
        });
    }

    getInputAsRenderedTemplates(context) {
        return this.inputAsTemplates.map((t) => t.getRenderedJSON(context));
    }

    get disableRenderMaterials() {
        return this.prop("isMultiMaterial", false);
    }

    getHashObject() {
        return {
            applicationName: this.applicationName,
            executableName: this.executableName,
            executableId: this.executableId,
            input: this.input,
            results: this.results,
            monitors: this.monitors,
            preProcessors: this.preProcessors,
            postProcessors: this.postProcessors,
        };
    }
}
