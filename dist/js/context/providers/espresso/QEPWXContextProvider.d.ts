import { type MethodDataContextMixinType, type Pseudo } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity/in_memory";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import { type JobContextMixinType } from "../../mixins/JobContextMixin";
import { type MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import { type MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import { type WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
import { ExecutableContextProvider } from "../../providers";
export type Material = MaterialMixin & InMemoryEntity;
export type Base = typeof ExecutableContextProvider & Constructor<MaterialContextMixinType> & Constructor<MaterialsContextMixinType> & Constructor<MethodDataContextMixinType> & Constructor<WorkflowContextMixinType> & Constructor<JobContextMixinType>;
declare const QEPWXContextProvider_base: Base;
export default class QEPWXContextProvider extends QEPWXContextProvider_base {
    static Material: typeof import("@mat3ra/made/dist/js/material").Material;
    _material?: Material;
    _materials: Material[];
    constructor(config: ContextProviderConfig);
    static atomSymbols(material: Material): string[];
    static uniqueElementsWithLabels(material: Material): string[];
    /** Returns the input text block for atomic positions WITH constraints.
     */
    static atomicPositionsWithConstraints(material: Material): string;
    /** Returns the input text block for atomic positions
     *  Note: does NOT include constraints
     */
    static atomicPositions(material: Material): string;
    static NAT(material: Material): number;
    static NTYP(material: Material): number;
    static NTYP_WITH_LABELS(material: Material): number;
    buildQEPWXContext(material: Material): {
        IBRAV: number;
        RESTART_MODE: string;
        ATOMIC_SPECIES: string;
        ATOMIC_SPECIES_WITH_LABELS: string;
        NAT: number;
        NTYP: number;
        NTYP_WITH_LABELS: number;
        ATOMIC_POSITIONS: string;
        ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: string;
        CELL_PARAMETERS: string;
    };
    getDataPerMaterial(): Record<string, unknown>;
    getData(): {
        IBRAV: number;
        RESTART_MODE: string;
        ATOMIC_SPECIES: string;
        ATOMIC_SPECIES_WITH_LABELS: string;
        NAT: number;
        NTYP: number;
        NTYP_WITH_LABELS: number;
        ATOMIC_POSITIONS: string;
        ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: string;
        CELL_PARAMETERS: string;
    };
    get RESTART_MODE(): "restart" | "from_scratch";
    getPseudoBySymbol(symbol: string): Pseudo | undefined;
    /** Builds ATOMIC SPECIES block of pw.x input in the format
     *  X   Mass_X   PseudoPot_X
     *  where X            is the atom label
     *        Mass_X       is the mass of element X [amu]
     *        PseudoPot_X  is the pseudopotential filename associated with element X
     *
     *  Note: assumes this.methodData is defined
     */
    ATOMIC_SPECIES(material: Material): string;
    ATOMIC_SPECIES_WITH_LABELS(material: Material): string;
    static CELL_PARAMETERS(material: Material): string;
    static symbolToAtomicSpecie(symbol: string, pseudo?: Pseudo): string;
    static elementAndPseudoToAtomicSpecieWithLabels(symbol: string, pseudo?: Pseudo, label?: string): string;
}
export {};
