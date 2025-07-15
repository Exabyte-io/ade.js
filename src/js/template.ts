import { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";

import {
    type TemplateMixin,
    type TemplateStaticMixin,
    templateMixin,
    templateStaticMixin,
} from "./templateMixin";

type Base = typeof NamedInMemoryEntity & Constructor<TemplateMixin> & TemplateStaticMixin;

export default class Template extends (NamedInMemoryEntity as Base) {}

// Apply mixins
templateMixin(Template.prototype);
templateStaticMixin(Template);
