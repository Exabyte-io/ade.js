/* eslint-disable max-classes-per-file */
import {
    Constructor,
    JobContextMixin,
    MaterialContextMixin,
    MaterialsContextMixin,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
} from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";

import { ExecutableContextProvider } from "../../providers";

const VASPContextProviderBase = JobContextMixin(
    WorkflowContextMixin(
        MethodDataContextMixin(
            MaterialsContextMixin<Constructor,Made.Material>(
                MaterialContextMixin<Constructor,Made.Material>(ExecutableContextProvider)
            )
        )
    )
);

export class VASPContextProvider extends VASPContextProviderBase {
    static Material = Made.Material;

    // eslint-disable-next-line class-methods-use-this
    buildVASPContext(material: Made.Material) {
        return {
            // TODO: figure out whether we need two separate POSCARS, maybe one is enough
            POSCAR: material.getAsPOSCAR(true, true),
            POSCAR_WITH_CONSTRAINTS: material.getAsPOSCAR(true),
        };
    }

    getDataPerMaterial() {
        if (!this.materials || this.materials.length <= 1) return {};
        return { perMaterial: this.materials.map((material: Made.Material) => this.buildVASPContext(material)) };
    }

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        // consider adjusting so that below values are read from PlanewaveDataManager
        // ECUTWFC;
        // ECUTRHO;

        return {
            ...this.buildVASPContext(this.material),
            ...this.getDataPerMaterial(),
        };
    }
}

const VASPNEBContextProviderBase = JobContextMixin(
    WorkflowContextMixin(
        MethodDataContextMixin(
            MaterialsSetContextMixin(
                MaterialsContextMixin<Constructor,Made.Material>(
                    MaterialContextMixin<Constructor,Made.Material>(ExecutableContextProvider)
                )
            )
        )
    )
);

export class VASPNEBContextProvider extends VASPNEBContextProviderBase {
    static Material = Made.Material;

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const VASPContexts = sortedMaterials.map((material: Partial<Made.Material>) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new VASPContextProvider(config).getData();
        });

        return {
            FIRST_IMAGE: VASPContexts[0].POSCAR_WITH_CONSTRAINTS,
            LAST_IMAGE: VASPContexts[VASPContexts.length - 1].POSCAR_WITH_CONSTRAINTS,
            INTERMEDIATE_IMAGES: VASPContexts.slice(1, VASPContexts.length - 1).map(
                (data: {POSCAR_WITH_CONSTRAINTS: string}) => data.POSCAR_WITH_CONSTRAINTS,
            ),
        };
    }
}
