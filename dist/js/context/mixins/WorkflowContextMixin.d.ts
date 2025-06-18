import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { WorkflowSchema } from "@mat3ra/esse/dist/js/types";
type Workflow = WorkflowSchema & {
    hasRelaxation: boolean;
};
export type WorkflowContextMixinType = {
    isEdited?: boolean;
    workflow: Workflow;
    _workflow: Workflow;
    initWorkflowContextMixin: () => void;
};
export declare function workflowContextMixin(item: ContextProvider): void;
export declare function WorkflowContextMixin<T extends Constructor<ContextProvider>>(superclass: T): T & Constructor<WorkflowContextMixinType>;
export {};
