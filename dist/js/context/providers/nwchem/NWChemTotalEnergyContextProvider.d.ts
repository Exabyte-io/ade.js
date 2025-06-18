import { type MethodDataContextMixinType, type Pseudo } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity/in_memory";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import { type JobContextMixinType } from "../../mixins/JobContextMixin";
import { type MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import { type WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
import { ExecutableContextProvider } from "../../providers";
export type Material = MaterialMixin & InMemoryEntity;
export type Base = typeof ExecutableContextProvider & Constructor<MaterialContextMixinType> & Constructor<MethodDataContextMixinType> & Constructor<WorkflowContextMixinType> & Constructor<JobContextMixinType>;
interface NWChemContextConfig extends ContextProviderConfig {
    context?: {
        material?: Material;
    };
}
declare const NWChemTotalEnergyContextProvider_base: Base;
export default class NWChemTotalEnergyContextProvider extends NWChemTotalEnergyContextProvider_base {
    static Material: typeof import("@mat3ra/made/dist/js/material").Material;
    _material?: Material;
    constructor(config: NWChemContextConfig);
    get atomicPositionsWithoutConstraints(): string[];
    get atomicPositions(): string[];
    get atomSymbols(): string[];
    get cartesianAtomicPositions(): boolean;
    get ATOMIC_SPECIES(): string;
    getData(): {
        CHARGE: number;
        MULT: number;
        BASIS: string;
        NAT: number;
        NTYP: number;
        ATOMIC_POSITIONS: string;
        ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: string;
        ATOMIC_SPECIES: string;
        FUNCTIONAL: string;
        CARTESIAN: boolean;
    };
    static symbolToAtomicSpecies(symbol: string, pseudo?: Pseudo): string | undefined;
}
export {};
