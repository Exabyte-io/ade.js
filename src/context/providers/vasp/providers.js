import { mix } from "mixwith";
import {
    MaterialContextMixinBuilder,
    MaterialsContextMixinBuilder,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin
} from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { ExecutableContextProvider } from "../../providers";

export class VASPContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixinBuilder(Made.Material),
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {

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
        }
    }

}

export class VASPNEBContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixinBuilder(Made.Material),
    MaterialsContextMixinBuilder(Made.Material),
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const VASPContexts = sortedMaterials.map(material => {
            const context = Object.assign({}, this.config.context, {material: material});
            const config = Object.assign({}, this.config, {context});
            return new VASPContextProvider(config).getData();
        });

        return {
            FIRST_IMAGE: VASPContexts[0].POSCAR_WITH_CONSTRAINTS,
            LAST_IMAGE: VASPContexts[VASPContexts.length - 1].POSCAR_WITH_CONSTRAINTS,
            INTERMEDIATE_IMAGES: VASPContexts.slice(1, VASPContexts.length - 1).map(data => data.POSCAR_WITH_CONSTRAINTS),
        }
    }
}
