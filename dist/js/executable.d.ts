import { NamedDefaultableInMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { type ExecutableMixin } from "./executableMixin";
type Base = Constructor<ExecutableMixin> & typeof NamedDefaultableInMemoryEntity;
declare const Executable_base: Base;
export default class Executable extends Executable_base {
}
export {};
