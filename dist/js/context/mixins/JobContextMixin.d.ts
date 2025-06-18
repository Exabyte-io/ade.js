import type { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { JobSchema } from "@mat3ra/esse/dist/js/types";
type Job = JobSchema & {
    parentJob?: Job;
};
export type JobContextMixinType = {
    isEdited?: boolean;
    job: Job;
    _job: Job;
    initJobContextMixin: () => void;
};
export declare function jobContextMixin(item: ContextProvider): void;
export declare function JobContextMixin<T extends Constructor<ContextProvider>>(superclass: T): T & Constructor<JobContextMixinType>;
export {};
