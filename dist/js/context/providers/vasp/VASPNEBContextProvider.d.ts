import type { MethodDataContextMixinType } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { JobContextMixinType } from "../../mixins/JobContextMixin";
import { type MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import { type MaterialsSetContextMixinType } from "../../mixins/MaterialsSetContextMixin";
import type { WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
import { ExecutableContextProvider } from "../../providers";
import type { Material } from "../espresso/QENEBContextProvider";
type Base = typeof ExecutableContextProvider & Constructor<MaterialContextMixinType> & Constructor<MaterialsContextMixinType> & Constructor<MaterialsSetContextMixinType> & Constructor<MethodDataContextMixinType> & Constructor<WorkflowContextMixinType> & Constructor<JobContextMixinType>;
declare const VASPNEBContextProvider_base: Base;
export default class VASPNEBContextProvider extends VASPNEBContextProvider_base {
    _materials: Material[];
    constructor(config: ContextProviderConfig);
    getData(): {
        FIRST_IMAGE: string;
        LAST_IMAGE: string;
        INTERMEDIATE_IMAGES: string[];
    };
}
export {};
