import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { WorkflowSchema } from "@mat3ra/esse/dist/js/types";

type Workflow = WorkflowSchema & {
    hasRelaxation: boolean;
};

type WorkflowConfig = {
    context?: {
        workflow?: Workflow;
    };
};

const defaultWorkflow: Workflow = {
    subworkflows: [],
    units: [],
    hasRelaxation: false,
};

export type WorkflowContextMixinType = {
    isEdited?: boolean;
    workflow: Workflow;
    _workflow: Workflow;
    initWorkflowContextMixin: () => void;
};

export function workflowContextMixin(item: ContextProvider) {
    const properties = {
        isEdited: false,

        _workflow: defaultWorkflow,

        get workflow() {
            return this._workflow;
        },

        initWorkflowContextMixin() {
            const config = this.config as WorkflowConfig;
            this._workflow = (config.context && config.context.workflow) || defaultWorkflow;
            this.isEdited = false; // we always get the `defaultData` (recalculated from scratch, not persistent)
        },
    } as WorkflowContextMixinType & typeof item;

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}

export function WorkflowContextMixin<T extends Constructor<ContextProvider>>(superclass: T) {
    workflowContextMixin(superclass.prototype);
    return superclass as T & Constructor<WorkflowContextMixinType>;
}
