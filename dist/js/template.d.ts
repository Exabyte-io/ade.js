import { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type TemplateMixin, type TemplateStaticMixin } from "./templateMixin";
type Base = typeof NamedInMemoryEntity & Constructor<TemplateMixin> & TemplateStaticMixin;
declare const Template_base: Base;
export default class Template extends Template_base {
}
export {};
