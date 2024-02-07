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
// @ts-expect-error periodic-table.js is not typed
import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";
import lodash from "lodash";
import path from "path";
import s from "underscore.string";

import { ExecutableContextProvider } from "../../providers";

const QEPWXContextProviderBase = JobContextMixin(
    WorkflowContextMixin(
        MethodDataContextMixin(
            MaterialsContextMixin<Constructor, Made.Material>(
                MaterialContextMixin<Constructor, Made.Material>(ExecutableContextProvider)
            )
        )
    )
)
export class QEPWXContextProvider extends QEPWXContextProviderBase {
    static Material = Made.Material;

    get materials() {
        return this._materials as Made.Material[];
    }

    static atomSymbols(material: Made.Material) {
        return material.Basis.uniqueElements;
    }

    /** Returns the input text block for atomic positions WITH constraints.
     */
    static atomicPositionsWithConstraints(material: Made.Material) {
        return material.Basis.atomicPositionsWithConstraints.join("\n");
    }

    /** Returns the input text block for atomic positions
     *  Note: does NOT include constraints
     */
    static atomicPositions(material: Made.Material) {
        return material.Basis.atomicPositions.join("\n");
    }

    static NAT(material: Made.Material) {
        return material.Basis.atomicPositions.length;
    }

    static NTYP(material: Made.Material) {
        return material.Basis.uniqueElements.length;
    }

    buildQEPWXContext(material: Made.Material) {
        const IBRAV = 0; // use CELL_PARAMETERS to define Bravais lattice

        return {
            IBRAV,
            RESTART_MODE: this.RESTART_MODE,
            ATOMIC_SPECIES: this.ATOMIC_SPECIES(material),
            NAT: QEPWXContextProvider.NAT(material),
            NTYP: QEPWXContextProvider.NTYP(material),
            ATOMIC_POSITIONS: QEPWXContextProvider.atomicPositionsWithConstraints(material),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: QEPWXContextProvider.atomicPositions(material),
            CELL_PARAMETERS: QEPWXContextProvider.CELL_PARAMETERS(material),
        };
    }

    getDataPerMaterial() {
        if (!this.materials || this.materials.length <= 1) return {};
        return { perMaterial: this.materials.map((material) => this.buildQEPWXContext(material)) };
    }

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        // the below values are read from PlanewaveDataManager instead
        // ECUTWFC = 40;
        // ECUTRHO = 200;

        return {
            ...this.buildQEPWXContext(this.material),
            ...this.getDataPerMaterial(),
        };
    }

    get RESTART_MODE() {
        return this.job.parentJob || this.workflow.hasRelaxation ? "restart" : "from_scratch";
    }

    getPseudoBySymbol(symbol: string) {
        return (this.methodData.pseudo || []).find((p: {element: string}) => p.element === symbol);
    }

    /** Builds ATOMIC SPECIES block of pw.x input in the format
     *  X   Mass_X   PseudoPot_X
     *  where X            is the atom label
     *        Mass_X       is the mass of element X [amu]
     *        PseudoPot_X  is the pseudopotential filename associated with element X
     *
     *  Note: assumes this.methodData is defined
     */
    ATOMIC_SPECIES(material: Made.Material) {
        return QEPWXContextProvider.atomSymbols(material)
            .map((symbol) => {
                const pseudo = this.getPseudoBySymbol(symbol);
                return QEPWXContextProvider.symbolToAtomicSpecie(symbol, pseudo);
            })
            .join("\n");
    }

    static CELL_PARAMETERS(material: Made.Material) {
        return material.Lattice.vectorArrays
            .map((x) => {
                return x
                    .map((y) => {
                        return s.sprintf("%14.9f", y).trim();
                    })
                    .join(" ");
            })
            .join("\n");
    }

    static symbolToAtomicSpecie(symbol: string, pseudo: { filename?: string, path?: string}) {
        const el = PERIODIC_TABLE[symbol];
        const filename = pseudo?.filename || path.basename(pseudo?.path || "");
        return el ? s.sprintf("%s %f %s", symbol, el.atomic_mass, filename) : undefined;
    }
}

const QENEBContextProviderBase = JobContextMixin(
    WorkflowContextMixin(
        MethodDataContextMixin(
            MaterialsSetContextMixin(
                MaterialsContextMixin<Constructor, Made.Material>(
                    MaterialContextMixin<Constructor, Made.Material>(ExecutableContextProvider)
                )
            )
        )
    )
)

export class QENEBContextProvider extends QENEBContextProviderBase {
    static Material = Made.Material;

    get material() {
        return this._material as Made.Material;
    }

    get materials() {
        return this._materials as Made.Material[];
    }

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const PWXContexts = sortedMaterials.map((material: Partial<Made.Material>) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new QEPWXContextProvider(config).getData();
        });

        return {
            ...lodash.omit(PWXContexts[0], [
                "ATOMIC_POSITIONS",
                "ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS",
            ]),
            FIRST_IMAGE: PWXContexts[0].ATOMIC_POSITIONS,
            LAST_IMAGE: PWXContexts[PWXContexts.length - 1].ATOMIC_POSITIONS,
            INTERMEDIATE_IMAGES: PWXContexts.slice(1, PWXContexts.length - 1).map(
                (data: {ATOMIC_POSITIONS: string}) => data.ATOMIC_POSITIONS,
            ),
        };
    }
}
