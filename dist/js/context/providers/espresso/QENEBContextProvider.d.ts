import { type MethodDataContextMixinType } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import type { JobContextMixinType } from "../../mixins/JobContextMixin";
import type { MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import { type MaterialsSetContextMixinType } from "../../mixins/MaterialsSetContextMixin";
import { type WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
import { ExecutableContextProvider } from "../../providers";
export type Base = typeof ExecutableContextProvider & Constructor<MaterialContextMixinType> & Constructor<MaterialsContextMixinType> & Constructor<MethodDataContextMixinType> & Constructor<WorkflowContextMixinType> & Constructor<JobContextMixinType> & Constructor<MaterialsSetContextMixinType>;
export type Material = MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet;
declare const QENEBContextProvider_base: Base;
export default class QENEBContextProvider extends QENEBContextProvider_base {
    static Material: typeof import("@mat3ra/made/dist/js/material").Material;
    _material?: Material;
    _materials: Material[];
    _materialsSet: undefined;
    constructor(config: ContextProviderConfig);
    getData(): {
        FIRST_IMAGE: string;
        LAST_IMAGE: string;
        INTERMEDIATE_IMAGES: string[];
        IBRAV: number;
        RESTART_MODE: string;
        ATOMIC_SPECIES: string;
        ATOMIC_SPECIES_WITH_LABELS: string;
        NAT: number;
        NTYP: number;
        NTYP_WITH_LABELS: number;
        CELL_PARAMETERS: string;
    };
}
export {};
