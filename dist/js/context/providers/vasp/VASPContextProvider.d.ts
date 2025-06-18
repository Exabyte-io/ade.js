import type { MethodDataContextMixinType } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { JobContextMixinType } from "../../mixins/JobContextMixin";
import { type MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import type { WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
import { ExecutableContextProvider } from "../../providers";
import type { Material } from "../espresso/QENEBContextProvider";
export type Base = typeof ExecutableContextProvider & Constructor<MaterialContextMixinType> & Constructor<MaterialsContextMixinType> & Constructor<MethodDataContextMixinType> & Constructor<WorkflowContextMixinType> & Constructor<JobContextMixinType>;
declare const VASPContextProvider_base: Base;
export default class VASPContextProvider extends VASPContextProvider_base {
    static Material: typeof import("@mat3ra/made/dist/js/material").Material;
    _material?: Material;
    _materials: Material[];
    constructor(config: ContextProviderConfig);
    buildVASPContext(material: Material): {
        POSCAR: string;
        POSCAR_WITH_CONSTRAINTS: string;
    };
    getDataPerMaterial(): {
        perMaterial?: undefined;
    } | {
        perMaterial: {
            POSCAR: string;
            POSCAR_WITH_CONSTRAINTS: string;
        }[];
    };
    getData(): {
        perMaterial?: undefined;
        POSCAR: string;
        POSCAR_WITH_CONSTRAINTS: string;
    } | {
        perMaterial: {
            POSCAR: string;
            POSCAR_WITH_CONSTRAINTS: string;
        }[];
        POSCAR: string;
        POSCAR_WITH_CONSTRAINTS: string;
    };
}
export {};
