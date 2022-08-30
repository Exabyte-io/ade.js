import {
    JobContextMixin,
    MaterialContextMixin,
    MaterialsContextMixin,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
} from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";

import { ExecutableContextProvider } from "../../providers";

export class VASPContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {
    static materialCls = Made.Material;

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        // consider adjusting so that below values are read from PlanewaveDataManager
        // ECUTWFC;
        // ECUTRHO;

        return {
            // TODO: figure out whether we need two separate POSCARS, maybe one is enough
            POSCAR: this.material.getAsPOSCAR(true, true),
            POSCAR_WITH_CONSTRAINTS: this.material.getAsPOSCAR(true),
        };
    }
}

export class VASPNEBContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixin,
    MaterialsContextMixin,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {
    static materialCls = Made.Material;

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const VASPContexts = sortedMaterials.map((material) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new VASPContextProvider(config).getData();
        });

        return {
            FIRST_IMAGE: VASPContexts[0].POSCAR_WITH_CONSTRAINTS,
            LAST_IMAGE: VASPContexts[VASPContexts.length - 1].POSCAR_WITH_CONSTRAINTS,
            INTERMEDIATE_IMAGES: VASPContexts.slice(1, VASPContexts.length - 1).map(
                (data) => data.POSCAR_WITH_CONSTRAINTS,
            ),
        };
    }
}
