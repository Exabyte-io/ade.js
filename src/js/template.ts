import { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import {
    type TemplateMixin,
    type TemplateStaticMixin,
    templateMixin,
    templateStaticMixin,
} from "./templateMixin";

type Base = typeof NamedInMemoryEntity & Constructor<TemplateMixin> & TemplateStaticMixin;

export class Template extends (NamedInMemoryEntity as Base) {
    context: Record<string, unknown> = {}; // TODO: this is a Typescript fix, but it didn't exist in the original code
}

// Apply mixins
templateMixin(Template.prototype);
templateStaticMixin(Template);
