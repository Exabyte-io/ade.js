import {
    JobContextMixin,
    MaterialContextMixin,
    MaterialsContextMixin,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
} from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";
import _ from "lodash";
import { mix } from "mixwith";
import s from "underscore.string";

import { ExecutableContextProvider } from "../../providers";

export class QEPWXContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixin,
    MaterialsContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {
    static Material = Made.Material;

    get atomSymbols() {
        return this.material.Basis.uniqueElements;
    }

    static atomSymbols(material) {
        return material.Basis.uniqueElements;
    }

    /** Returns the input text block for atomic positions WITH constraints.
     */
    static atomicPositionsWithConstraints(material) {
        return material.Basis.atomicPositionsWithConstraints.join("\n");
    }

    /** Returns the input text block for atomic positions
     *  Note: does NOT include constraints
     */
    static atomicPositions(material) {
        return material.Basis.atomicPositions.join("\n");
    }

    static getMaterialContext(material) {
        return {
            NAT: material.Basis.atomicPositions.length,
            NTYP: material.Basis.uniqueElements.length,
            ATOMIC_POSITIONS: QEPWXContextProvider.atomicPositionsWithConstraints(material),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: QEPWXContextProvider.atomicPositions(material),
            CELL_PARAMETERS: QEPWXContextProvider.CELL_PARAMETERS(material),
        };
    }

    static getMaterialsContext(materials) {
        if (!materials) return {};
        return { perMaterial: materials.map((material) => this.getMaterialContext(material)) };
    }

    static NAT(material) {
        return material.Basis.atomicPositions.length;
    }

    static NTYP(material) {
        return material.Basis.uniqueElements.length;
    }

    buildQEPWXContext(material) {
        const IBRAV = 0; // use CELL_PARAMETERS to define Bravais lattice

        return {
            IBRAV,
            RESTART_MODE: this.RESTART_MODE,
            ATOMIC_SPECIES: this._ATOMIC_SPECIES(material),
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

    getPseudoBySymbol(symbol) {
        return (this.methodData.pseudo || []).find((p) => p.element === symbol);
    }

    get ATOMIC_SPECIES() {
        // atomic species with pseudopotentials
        return _.map(this.atomSymbols, (symbol) => {
            const pseudo = this.getPseudoBySymbol(symbol);
            return QEPWXContextProvider.symbolToAtomicSpecie(symbol, pseudo);
        }).join("\n");
    }

    /** Builds ATOMIC SPECIES block of pw.x input in the format
     *  X   Mass_X   PseudoPot_X
     *  where X            is the atom label
     *        Mass_X       is the mass of element X [amu]
     *        PseudoPot_X  is the pseudopotential filename associated with element X
     *
     *  Note: assumes this.methodData is defined
     */
    _ATOMIC_SPECIES(material) {
        return QEPWXContextProvider.atomSymbols(material)
            .map((symbol) => {
                const pseudo = this.getPseudoBySymbol(symbol);
                return QEPWXContextProvider.symbolToAtomicSpecie(symbol, pseudo);
            })
            .join("\n");
    }

    static CELL_PARAMETERS(material) {
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

    static symbolToAtomicSpecie(symbol, pseudo) {
        const el = PERIODIC_TABLE[symbol];
        const filename = pseudo ? _.get(pseudo, "filename", s.strRightBack(pseudo.path, "/")) : "";
        return el ? s.sprintf("%s %f %s", symbol, el.atomic_mass, filename) : undefined;
    }
}

export class QENEBContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixin,
    MaterialsContextMixin,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {
    static Material = Made.Material;

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const PWXContexts = sortedMaterials.map((material) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new QEPWXContextProvider(config).getData();
        });

        return {
            ..._.omit(PWXContexts[0], ["ATOMIC_POSITIONS", "ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS"]),
            FIRST_IMAGE: PWXContexts[0].ATOMIC_POSITIONS,
            LAST_IMAGE: PWXContexts[PWXContexts.length - 1].ATOMIC_POSITIONS,
            INTERMEDIATE_IMAGES: PWXContexts.slice(1, PWXContexts.length - 1).map(
                (data) => data.ATOMIC_POSITIONS,
            ),
        };
    }
}
