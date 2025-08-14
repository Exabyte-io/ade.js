"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const entity_1 = require("@mat3ra/code/dist/js/entity");
const RuntimeItemsMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin");
const executableMixin_1 = require("./executableMixin");
class Executable extends entity_1.NamedDefaultableInMemoryEntity {
}
exports.default = Executable;
// Apply mixins
(0, RuntimeItemsMixin_1.runtimeItemsMixin)(Executable.prototype);
(0, executableMixin_1.executableMixin)(Executable.prototype);
(0, executableMixin_1.executableStaticMixin)(Executable);
